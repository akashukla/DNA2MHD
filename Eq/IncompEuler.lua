-- Gkyl ------------------------------------------------------------------------
--
-- 2D incompressible Euler equations
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local IncompEulerModDecl = require "Eq.poissonBracketData.CanonicalModDecl"
local Proto = require "Lib.Proto"
local HamiltonianBase = require "Eq.HamiltonianBase"
local ffi = require "ffi"
local xsys = require "xsys"

-- start from HamiltonianBase class
local IncompEuler = Proto(HamiltonianBase)

-- ctor
function IncompEuler:init(tbl)
   -- initialize generic Hamiltonian equation
   IncompEuler.super.init(self, tbl)
   self:initHamilTimeIndep()  -- equation dependent

   -- store pointers to C kernels implementing volume and surface terms
   local nm, ndim, p = self._basis:id(), self._basis:ndim(), self._basis:polyOrder()
   assert(ndim == 2, "Incompressible Euler equations only implemented in 2D")
   self._volTerm = IncompEulerModDecl.selectVol(nm, ndim, p)
   self._surfTerms = IncompEulerModDecl.selectSurf(nm, ndim, p)
end

function IncompEuler:initHamilTimeIndep()
   -- no time-independent part of hamiltonian for incomp euler
   self.hamilTimeIndep:clear(0.0)
end

function IncompEuler:setAuxFields(auxFields)
   -- get streamfunction, psi
   self.psi = auxFields[1]

   -- for incomp euler system, time-dependent part of hamiltonian is just psi
   self:setHamiltonian(self.psi)
end

return IncompEuler