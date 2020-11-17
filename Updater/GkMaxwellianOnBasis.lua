-- Gkyl ------------------------------------------------------------------------
--
-- Updater to create the Maxwellian distribution from the coserved
-- moments and project it on basis functions. Uses Gaussian
-- quadrature.
--
--------------------------------------------------------------------------------

-- Gkyl libraries
local GaussQuadRules = require "Lib.GaussQuadRules"
local LinearDecomp = require "Lib.LinearDecomp"
local Lin = require "Lib.Linalg"
local Proto = require "Proto"
local Range = require "Lib.Range"
local Time = require "Lib.Time"
local UpdaterBase = require "Updater.Base"

-- System libraries
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"

ffi.cdef [[
  void GkMaxwellianInnerLoop(double * n, double * u, double * vtSq,
                           double * bmag, double m_,
			   double * fItr,
			   double * weights, double * dz, double * zc,
			   double * ordinates,
                           double * basisAtOrdinates,
                           double * phaseToConfOrdMap,
			   int numPhaseBasis, 
                           int numConfOrds, int numPhaseOrds,
			   int numConfDims, int numPhaseDims);

]]

-- Inherit the base Updater from UpdaterBase updater object
local GkMaxwellianOnBasis = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function GkMaxwellianOnBasis:init(tbl)
   GkMaxwellianOnBasis.super.init(self, tbl) -- setup base object

   self.confGrid = assert(tbl.confGrid,
			  "Updater.GkMaxwellianOnBasis: Must provide configuration space grid object 'confGrid'")
   self.confBasis = assert(tbl.confBasis,
			   "Updater.GkMaxwellianOnBasis: Must provide configuration space basis object 'confBasis'")
   self.phaseGrid = assert(tbl.onGrid,
			   "Updater.GkMaxwellianOnBasis: Must provide phase space grid object 'onGrid'")
   self.phaseBasis = assert(tbl.phaseBasis,
			    "Updater.GkMaxwellianOnBasis: Must provide phase space basis object 'phaseBasis'")
   self.mass = assert(tbl.gkfacs[1],
		      "Updater.GkMaxwellianOnBasis: Must provide mass object in 'gkfacs'")
   self.bmag = assert(tbl.gkfacs[2],
		      "Updater.GkMaxwellianOnBasis: Must provide bmag object in 'gkfacs'")
   self.bmagItr = self.bmag:get(1)

   -- Number of quadrature points in each direction
   local N = tbl.numConfQuad and tbl.numConfQuad or self.confBasis:polyOrder() + 1
   assert(N<=8, "Updater.GkMaxwellianOnBasis: Gaussian quadrature only implemented for numQuad<=8 in each dimension")

   self.onGhosts = xsys.pickBool(tbl.projectOnGhosts, true)

   -- 1D weights and ordinates
   local ordinates = GaussQuadRules.ordinates[N]
   local weights = GaussQuadRules.weights[N]

   -- Configuration space ordinates ----------------------------------
   self._cDim = self.confBasis:ndim()
   local l, u = {}, {}
   for d = 1, self._cDim do l[d], u[d] = 1, N end
   self.confQuadRange = Range.Range(l, u)
   self.confQuadIndexer = 
      Range.makeColMajorGenIndexer(self.confQuadRange)
   self.numConfOrds = self.confQuadRange:volume()
   self.numConfBasis = self.confBasis:numBasis()
   self.confOrdinates = Lin.Mat(self.numConfOrds, self._cDim)
   self.confBasisAtOrds = Lin.Mat(self.numConfOrds, self.numConfBasis)
   for ordIndexes in self.confQuadRange:colMajorIter() do
      local ordIdx = self.confQuadIndexer(ordIndexes)
      for d = 1, self._cDim do
	 self.confOrdinates[ordIdx][d] = ordinates[ordIndexes[d]]
      end
      self.confBasis:evalBasis(self.confOrdinates[ordIdx],
			       self.confBasisAtOrds[ordIdx])
   end

   -- Phase space ordinates and weights ------------------------------
   self._pDim = self.phaseBasis:ndim()
   for d = 1, self._pDim do l[d], u[d] = 1, N end
   self.phaseQuadRange = Range.Range(l, u)
   self.phaseQuadIndexer = 
      Range.makeColMajorGenIndexer(self.phaseQuadRange)
   self.numPhaseOrds = self.phaseQuadRange:volume()
   self.numPhaseBasis = self.phaseBasis:numBasis()
   self.phaseBasisAtOrds = Lin.Mat(self.numPhaseOrds,
				   self.numPhaseBasis)
   self.phaseOrdinates = Lin.Mat(self.numPhaseOrds, self._pDim)
   self.phaseWeights = Lin.Vec(self.numPhaseOrds) -- Needed for integration
   for ordIndexes in self.phaseQuadRange:colMajorIter() do
      local ordIdx = self.phaseQuadIndexer(ordIndexes)
      self.phaseWeights[ordIdx] = 1.0
      for d = 1, self._pDim do
	 self.phaseWeights[ordIdx] =
	    self.phaseWeights[ordIdx]*weights[ordIndexes[d]]
	 self.phaseOrdinates[ordIdx][d] = ordinates[ordIndexes[d]]
      end
      self.phaseBasis:evalBasis(self.phaseOrdinates[ordIdx],
				self.phaseBasisAtOrds[ordIdx])
   end

   -- Construct the phase space to conf space ordinate map
   self.phaseToConfOrdMap = Lin.Vec(self.numPhaseOrds)
   for ordIndexes in self.phaseQuadRange:colMajorIter() do
      local confOrdIdx = self.confQuadIndexer(ordIndexes)
      local phaseOrdIdx = self.phaseQuadIndexer(ordIndexes)
      self.phaseToConfOrdMap[phaseOrdIdx] = confOrdIdx
   end

   self._vDim = self._pDim - self._cDim

   -- Cell index, center, and dx.
   self.idxP = Lin.IntVec(self._pDim)
   self.xcP  = Lin.Vec(self._pDim)
   self.dxP  = Lin.Vec(self._pDim)
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function GkMaxwellianOnBasis:_advance(tCurr, inFld, outFld)

   -- Get the inputs and outputs
   local nIn    = assert(inFld[1], "GkMaxwellianOnBasis.advance: Must specify density 'inFld[1]'")
   local uIn = assert(inFld[2], "GkMaxwellianOnBasis.advance: Must specify drift speed 'inFld[2]'")
   local vtSqIn = assert(inFld[3], "GkMaxwellianOnBasis.advance: Must specify thermal velocity squared 'inFld[3]'")
   local fOut   = assert(outFld[1], "GkMaxwellianOnBasis.advance: Must specify an output field 'outFld[1]'")

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim
   local uDim = uIn:numComponents()/self.numConfBasis -- number of dimensions in u
   
   -- OLD CODE before projecting on ghosts --
   -- Move these declarations below self.onGhosts 
   -- local bmagItr, bmagOrd = self.bmag:get(1), Lin.Vec(self.numConfOrds)
   -- local nItr, nOrd       = nIn:get(1), Lin.Vec(self.numConfOrds)
   -- local uItr, uOrd = uIn:get(1), Lin.Vec(self.numConfOrds)
   -- local vtSqItr, vtSqOrd = vtSqIn:get(1), Lin.Vec(self.numConfOrds)
   -- local fItr             = fOut:get(1)
 
   -- Get the Ranges to loop over the domain (always project on ghosts)
   -- local confRange    = nIn:localExtRange()
   -- local confIndexer  = nIn:genIndexer()
   -- local phaseRange   = fOut:localExtRange()
   -- local phaseIndexer = fOut:genIndexer()

   -- -- construct ranges for nested loops
   -- local confRangeDecomp = LinearDecomp.LinearDecompRange {
   --    range = phaseRange:selectFirst(cDim), numSplit = self.phaseGrid:numSharedProcs() }
   --    -- range = nIn:localExtRange(), numSplit = self.phaseGrid:numSharedProcs() }
   -- local velRange = phaseRange:selectLast(vDim)
   -- local tId = self.phaseGrid:subGridSharedId() -- local thread ID

   local ordIdx = nil
   local phaseRange = fOut:localRange()
   if self.onGhosts then -- extend range to config-space ghosts
      local cdirs = {}
      for dir = 1, cDim do 
         phaseRange = phaseRange:extendDir(dir, fOut:lowerGhost(), fOut:upperGhost())
      end
   end

   --local bmagItr, bmagOrd = self.bmag:get(1), Lin.Vec(self.numConfOrds)
   local bmagOrd = Lin.Vec(self.numConfOrds)
   local nItr, nOrd = nIn:get(1), Lin.Vec(self.numConfOrds)
   local uItr, uOrd = uIn:get(1), Lin.Vec(self.numConfOrds)
   local vtSqItr, vtSqOrd = vtSqIn:get(1), Lin.Vec(self.numConfOrds)
   local fItr             = fOut:get(1)

   -- Construct ranges for nested loops.
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phaseRange:selectFirst(cDim), numSplit = self.phaseGrid:numSharedProcs() }
   local velRange = phaseRange:selectLast(vDim)
   local tId      = self.phaseGrid:subGridSharedId()    -- Local thread ID.

   local confIndexer = nIn:genIndexer()
   local phaseIndexer = fOut:genIndexer()
   
   -- The configuration space loop
   for cIdx in confRangeDecomp:rowMajorIter(tId) do
      self.bmag:fill(confIndexer(cIdx), self.bmagItr)
      nIn:fill(confIndexer(cIdx), nItr)
      uIn:fill(confIndexer(cIdx), uItr)
      vtSqIn:fill(confIndexer(cIdx), vtSqItr)
      
      -- Evaluate the primitive variables (given as expansion
      -- coefficients) on the ordinates
      for ordIndexes in self.confQuadRange:rowMajorIter() do
	 ordIdx = self.confQuadIndexer(ordIndexes)
	 bmagOrd[ordIdx],nOrd[ordIdx],uOrd[ordIdx],vtSqOrd[ordIdx] = 0.0, 0.0, 0.0, 0.0
	 
	 for k = 1, self.numConfBasis do
	    bmagOrd[ordIdx] = bmagOrd[ordIdx] + self.bmagItr[k]*self.confBasisAtOrds[ordIdx][k]
	    nOrd[ordIdx] = nOrd[ordIdx] + nItr[k]*self.confBasisAtOrds[ordIdx][k]
	    vtSqOrd[ordIdx] = vtSqOrd[ordIdx] + vtSqItr[k]*self.confBasisAtOrds[ordIdx][k]

	    if uDim > 1 and cDim == 3 then -- get z-component of fluid velocity
	       uOrd[ordIdx] = uOrd[ordIdx] + uItr[self.numConfBasis*2+k]*self.confBasisAtOrds[ordIdx][k]
	    else 
	       uOrd[ordIdx] = uOrd[ordIdx] + uItr[k]*self.confBasisAtOrds[ordIdx][k]
	    end

	 end
      end

      -- The velocity space loop
      for vIdx in velRange:rowMajorIter() do
	 -- Construct the phase space index out of the configuration
	 -- space a velocity space indices
         cIdx:copyInto(self.idxP)
         for d = 1, vDim do self.idxP[cDim+d] = vIdx[d] end
	 fOut:fill(phaseIndexer(self.idxP), fItr)

	 -- Get cell shape, cell center coordinates
	 self.phaseGrid:setIndex(self.idxP)
         self.phaseGrid:getDx(self.dxP)
         self.phaseGrid:cellCenter(self.xcP)

	 ffiC.GkMaxwellianInnerLoop(nOrd:data(), uOrd:data(), vtSqOrd:data(),
				    bmagOrd:data(), self.mass,
				    fItr:data(),
				    self.phaseWeights:data(), self.dxP:data(), self.xcP:data(),
				    self.phaseOrdinates:data(),
				    self.phaseBasisAtOrds:data(),
				    self.phaseToConfOrdMap:data(),
				    self.numPhaseBasis,
				    self.numConfOrds, self.numPhaseOrds,
				    cDim, pDim)
      end

   end

   -- set id of output to id of projection basis
   fOut:setBasisId(self.phaseBasis:id())
end

return GkMaxwellianOnBasis
