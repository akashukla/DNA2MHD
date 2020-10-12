-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate a cell-averaged Maxwellian
--
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase = require "Updater.Base"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto = require "Lib.Proto"
local CellAveMaxwellianDecl = require "Updater.cellAveMaxwellianCalcData.MaxwellianCellAvModDecl"
local xsys = require "xsys"
local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"

-- Charge exchange collisions updater object.
local CellAveMaxwellian = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function CellAveMaxwellian:init(tbl)
   CellAveMaxwellian.super.init(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid,
			     "Updater.CellAvMaxw: Must provide grid object using 'onGrid'")
   self._confBasis = assert(tbl.confBasis,
			     "Updater.CellAvMax: Must provide configuration space basis object using 'confBasis'")
   self._phaseBasis = assert(tbl.phaseBasis,
			     "Updater.CellAvMaxw: Must provide phase space basis object using 'phaseBasis'")
   self._kineticSpecies = assert(tbl.kineticSpecies,
				 "Updater.CellAvMaxwellian: Must provide solver type (Vm or Gk) using 'kineticSpecies'")
   if self._kineticSpecies == 'Gk' then
      self.mass = assert(tbl.gkfacs[1],
			 "Updater.CellAvMaxwellian: Must provide mass object in 'gkfacs'")
      self.bmag = assert(tbl.gkfacs[2],
			 "Updater.CellAvMaxwellian: Must provide bmag object in 'gkfacs'")
      self.bmagItr = self.bmag:get(1)
   end
   
   -- Dimension of spaces.
   self._cDim = self._confBasis:ndim()
   self._pDim = self._phaseBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Basis name and polynomial order.
   self._basisID = self._phaseBasis:id()
   self._polyOrder = self._phaseBasis:polyOrder()

   -- Number of basis functions.
   self._numBasis = self._confBasis:numBasis()

   -- Define cell center
   self.xc  = Lin.Vec(self._pDim)

   -- Define functions
   if self._kineticSpecies == 'Vm' then
      self._calcMax = CellAveMaxwellianDecl.CellAvMax(self._basisID, self._cDim, self._vDim, self._polyOrder)
   else
      self._calcGkMax = CellAveMaxwellianDecl.GkCellAvMax(self._basisID, self._cDim, self._vDim, self._polyOrder)
   end
   
   self.onGhosts = xsys.pickBool(false, tbl.onGhosts)

   self._tmEvalMom = 0.0
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function CellAveMaxwellian:vlasov(m0, u, vtSq, fMax)
   local tmEvalMomStart = Time.clock()
   local grid = self._onGrid
   local pDim = self._pDim
   
   local confIndexer = m0:genIndexer()
   local phaseIndexer = fMax:genIndexer()
   
   local m0Itr = m0:get(1)
   local uItr = u:get(1)
   local vtSqItr = vtSq:get(1)
   local fMaxItr = fMax:get(1)
   
   local phaseRange = fMax:localRange()

   local phaseRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phaseRange:selectFirst(pDim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId()    -- Local thread ID.
   
   -- Phase space loop
   for pIdx in phaseRangeDecomp:rowMajorIter(tId) do
      grid:setIndex(pIdx)
      grid:cellCenter(self.xc)

      m0:fill(confIndexer(pIdx), m0Itr)
      u:fill(confIndexer(pIdx), uItr)
      vtSq:fill(confIndexer(pIdx), vtSqItr)
      fMax:fill(phaseIndexer(pIdx), fMaxItr)

      self._calcMax(self.xc:data(), m0Itr:data(), uItr:data(), vtSqItr:data(), fMaxItr:data())     
   end
   
   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
end

function CellAveMaxwellian:gyrokinetic(m0, u, vtSq, bmag, fMax)
   local tmEvalMomStart = Time.clock()
   local grid = self._onGrid
   local pDim = self._pDim
   
   local confIndexer = m0:genIndexer()
   local phaseIndexer = fMax:genIndexer()
   
   local m0Itr = m0:get(1)
   local uItr = u:get(1)
   local vtSqItr = vtSq:get(1)
   local fMaxItr = fMax:get(1)
   
   local phaseRange = fMax:localRange()

   local phaseRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phaseRange:selectFirst(pDim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId()    -- Local thread ID.
   
   -- Phase space loop
   for pIdx in phaseRangeDecomp:rowMajorIter(tId) do
      grid:setIndex(pIdx)
      grid:cellCenter(self.xc)

      m0:fill(confIndexer(pIdx), m0Itr)
      u:fill(confIndexer(pIdx), uItr)
      vtSq:fill(confIndexer(pIdx), vtSqItr)
      self.bmag:fill(confIndexer(pIdx), self.bmagItr)
      fMax:fill(phaseIndexer(pIdx), fMaxItr)

      self._calcGkMax(self.mass, self.xc:data(), m0Itr:data(), uItr:data(), vtSqItr:data(), self.bmagItr:data(), fMaxItr:data())     
   end
   
   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
end

function CellAveMaxwellian:_advance(tCurr, inFld, outFld)

   if self._kineticSpecies == 'Vm' then
      local m0 = assert(inFld[1], "cellAveMaxwellian: Must specify particle density as input[1]")
      local u = assert(inFld[2], "cellAveMaxwellian: Must specify fluid velocity as input[2]")
      local vtSq = assert(inFld[3], "cellAveMaxwellian: Must specify squared thermal velocity as input[3]")
      local fMax = assert(outFld[1], "cellAveMaxwellian: Must specify an output field")
      self:vlasov(m0, u, vtSq, fMax)
   elseif self._kineticSpecies == 'Gk' then
      local m0 = assert(inFld[1], "GkCellAveMaxwellian: Must specify particle density as input[1]")
      local u = assert(inFld[2], "GkCellAveMaxwellian: Must specify fluid velocity as input[2]")
      local vtSq = assert(inFld[3], "GkCellAveMaxwellian: Must specify squared thermal velocity as input[3]")
      local fMax = assert(outFld[1], "GkCellAveMaxwellian: Must specify an output field")
      self:gyrokinetic(m0, u, vtSq, bmag, fMax)
   else
      print("Updater.CellAveMaxwellian: kineticSpecies must be 'Vm' or 'Gk'")
   end
   
end
   
function CellAveMaxwellian:evalMomTime() return self._tmEvalMom end

return CellAveMaxwellian
