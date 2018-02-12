-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute integrated quantities on a Cartesian grid.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local Lin = require "Lib.Linalg"
local MomDecl = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"
local ffi = require "ffi"

ffi.cdef [[
    void gkylCartFieldIntQuantV(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out);
    void gkylCartFieldIntQuantV2(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out);
]]

-- Integrated quantities calculator
local CartFieldIntegratedQuantCalc = Proto(UpdaterBase)

function CartFieldIntegratedQuantCalc:init(tbl)
   CartFieldIntegratedQuantCalc.super.init(self, tbl) -- setup base object

   -- grid and basis
   self.onGrid = assert(
      tbl.onGrid, "Updater.CartFieldIntegratedQuantCalc: Must provide grid object using 'onGrid'")
   self.basis = assert(
      tbl.basis,
      "Updater.CartFieldIntegratedQuantCalc: Must provide phase-space basis object using 'basis'")

   -- number of components to set
   self.numComponents = tbl.numComponents and tbl.numComponents or 1

   self.updateFunc = ffi.C.gkylCartFieldIntQuantV
   if tbl.quantity == "V2" then
      self.updateFunc = ffi.C.gkylCartFieldIntQuantV2
   end

   -- for use in advance method
   self.dxv = Lin.Vec(self.basis:ndim()) -- cell shape
   self.localVals = Lin.Vec(self.numComponents)
   self.globalVals = Lin.Vec(self.numComponents)
end   

-- advance method
function CartFieldIntegratedQuantCalc:_advance(tCurr, dt, inFld, outFld)
   local grid = self.onGrid
   local field, vals = inFld[1], outFld[1]

   local ndim = self.basis:ndim()
   local nvals = #self.localVals

   local fieldIndexer = field:genIndexer()
   local fieldItr = field:get(1)

   -- clear local values
   for i = 1, #self.localVals do
      self.localVals[i] = 0.0
   end

   local localRange = field:localRange()   
   -- loop, computing integrated moments in each cell
   for idx in localRange:colMajorIter() do
      grid:setIndex(idx)
      for d = 1, ndim do
	 self.dxv[d] = grid:dx(d)
      end

      field:fill(fieldIndexer(idx), fieldItr)
      -- compute integrated quantities
      self.updateFunc(
	 ndim, self.numComponents, self.basis:numBasis(), self.dxv:data(), fieldItr:data(), self.localVals:data())
   end

   -- all-reduce across processors and push result into dyn-vector
   local nodeComm = self:getNodeComm()
   Mpi.Allreduce(self.localVals:data(), self.globalVals:data(), nvals, Mpi.DOUBLE, Mpi.SUM, nodeComm)
   vals:appendData(tCurr, self.globalVals)
   
   return true, GKYL_MAX_DOUBLE
end

return CartFieldIntegratedQuantCalc