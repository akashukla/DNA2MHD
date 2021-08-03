-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Ionization operator using Voronov
-- reaction rate.
-- See:
-- Voronov, G.S., 1997. A practical fit formula for ionization rate
-- coefficients of atoms and ions by electron impact: Z= 1–28. Atomic
-- Data and Nuclear Data Tables, 65(1), pp.1-35.
--
--------------------------------------------------------------------------------

local CollisionsBase = require "App.Collisions.CollisionsBase"
local Constants      = require "Lib.Constants"
local DataStruct     = require "DataStruct"
local DiagsImplBase  = require "App.Diagnostics.DiagnosticsImplBase"
local DiagsApp       = require "App.Diagnostics.SpeciesDiagnostics"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"
local Mpi            = require "Comm.Mpi"
local lume           = require "Lib.lume"
local xsys           = require "xsys"

-- GkIonization -----------------------------------------------------------
--
-- Voronov ionization operator.
---------------------------------------------------------------------------

-- ............... IMPLEMENTATION OF DIAGNOSTICS ................. --
-- Diagnostics could be placed in a separate file if they balloon in
-- number. But if we only have one or two we can just place it here.

-- ~~~~ Source integrated over the domain ~~~~~~~~~~~~~~~~~~~~~~
local gkIzDiagImpl = function()
   local _intSrcIzM0 = Proto(DiagsImplBase)
   function _intSrcIzM0:fullInit(diagApp, mySpecies, fieldIn, owner)
      self.fieldAux = mySpecies:allocMoment()
      self.updatersAux = mySpecies.numDensityCalc
      self.field    = DataStruct.DynVector { numComponents = 1 }
      self.updater  = mySpecies.volIntegral.scalar
      self.owner    = owner
      self.done     = false
   end
   function _intSrcIzM0:getType() return "integrated" end
   function _intSrcIzM0:advance(tm, inFlds, outFlds)
      self.updatersAux:advance(tm, {self.owner.ionizSrc}, {self.fieldAux})
      self.updater:advance(tm, {self.fieldAux}, {self.field})
   end

   return {intSrcIzM0 = _intSrcIzM0}
end

-- .................... END OF DIAGNOSTICS ...................... --

local GkIonization = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkIonization:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkIonization:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously store table.

   self.cfl = 0.1
   self.collKind = "Ionization"

   self.collidingSpecies = assert(tbl.collideWith, "App.GkIonization: Must specify names of species to collide with in 'collideWith'.")

   -- Set these values to be consistent with other collision apps
   self.selfCollisions  = false
   self.crossCollisions = true              
   self.varNu           = false
   self.timeDepNu       = false
   self.collFreqs       = {1}
   
   self.collideNm   = tbl.collideWith[1]
   
   self.elcNm       = assert(tbl.electrons, "App.GkIonization: Must specify electron species name in 'electrons'.")
   self.neutNm      = assert(tbl.neutrals, "App.GkIonization: Must specify electron species name in 'neutrals'.")
   
   self.plasma      = tbl.plasma
   self.mass        = tbl.elcMass
   self.charge      = tbl.elemCharge

   if self.plasma == "H" then
      self._E = 13.6
      self._P = 0
      self._A = 0.291e-7
      self._K = 0.39
      self._X = 0.232
   end

   if self.plasma == "Ar" then
      self._E = 15.8
      self._P = 1
      self._A = 0.599e-7
      self._K = 0.26
      self._X = 0.136
   end

   self.timers = {nonSlvr = 0.}
end

function GkIonization:createDiagnostics(mySpecies, field)
   -- Create source diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = gkIzDiagImpl()}
      self.diagnostics:fullInit(mySpecies, field, self)
   end
   return self.diagnostics
end

function GkIonization:setName(nm) self.name = self.speciesName.."_"..nm end
function GkIonization:setSpeciesName(nm) self.speciesName = nm end
function GkIonization:setCfl(cfl) self.cfl = cfl end
function GkIonization:setConfBasis(basis) self.confBasis = basis end
function GkIonization:setConfGrid(grid) self.confGrid = grid end
function GkIonization:setPhaseBasis(basis) self.phaseBasis = basis end
function GkIonization:setPhaseGrid(grid) self.phaseGrid = grid end

function GkIonization:createSolver(funcField)
   self.collisionSlvr = Updater.Ionization {
      onGrid     = self.confGrid,
      confBasis  = self.confBasis,
      phaseGrid  = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      elcMass    = self.mass,
      elemCharge = self.charge,
      reactRate  = true,
	 
      -- Voronov parameters
      A = self._A,
      E = self._E,
      K = self._K,
      P = self._P,
      X = self._X,
   }
   if (self.speciesName == self.elcNm) then
      self.calcIonizationTemp = Updater.Ionization {
	 onGrid     = self.confGrid,
	 confBasis  = self.confBasis,
	 phaseGrid  = self.phaseGrid,
	 phaseBasis = self.phaseBasis,
      	 elcMass    = self.mass,
      	 elemCharge = self.charge,
	 reactRate  = false, 
      	 E          = self._E,
      }
      self.sumDistF =  DataStruct.Field {
	 onGrid        = self.phaseGrid,
	 numComponents = self.phaseBasis:numBasis(),
	 ghost         = {1, 1},
	 metaData = {
	    polyOrder = self.phaseBasis:polyOrder(),
	    basisType = self.phaseBasis:id()
	 },
      }
   end
   self.confMult = Updater.CartFieldBinOp {
      onGrid     = self.confGrid,
      weakBasis  = self.confBasis,
      operation  = "Multiply",
   }
   self.confPhaseMult = Updater.CartFieldBinOp {
      onGrid     = self.phaseGrid,
      weakBasis  = self.phaseBasis,
      fieldBasis = self.confBasis,
      operation  = "Multiply",
   }
   self.m0elc = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.confBasis:polyOrder(),
	 basisType = self.confBasis:id()
      },
   }
   self.coefM0 = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.confBasis:polyOrder(),
	 basisType = self.confBasis:id()
      },
   }
   self.fMaxNeut = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
   self.ionizSrc = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
      metaData = {
	 polyOrder = self.phaseBasis:polyOrder(),
	 basisType = self.phaseBasis:id()
      },
   }
end

function GkIonization:advance(tCurr, fIn, species, fRhsOut)
   local tmNonSlvrStart = Time.clock()
   local coefIz = species[self.elcNm]:getVoronovReactRate()
   local elcM0  = species[self.elcNm]:fluidMoments()[1]

   -- Check whether particle is electron, neutral or ion species.
   if (self.speciesName == self.elcNm) then
      -- Electrons.
      local neutM0     = species[self.neutNm]:fluidMoments()[1]
      local elcDistF   = species[self.elcNm]:getDistF()
      local fMaxwellIz = species[self.elcNm]:getFMaxwellIz()
      
      self.sumDistF:combine(2.0,fMaxwellIz,-1.0,elcDistF)
      
      self.confMult:advance(tCurr, {coefIz, neutM0}, {self.coefM0})
      self.confPhaseMult:advance(tCurr, {self.coefM0, self.sumDistF}, {self.ionizSrc})
      -- Uncomment to test without fMaxwellian(Tiz)
      --self.confPhaseMult:advance(tCurr, {self.coefM0, elcDistF}, {self.ionizSrc})      

      fRhsOut:accumulate(1.0,self.ionizSrc)
   elseif (species[self.speciesName].charge == 0) then
      -- Neutrals, on Vlasov grid.
      local neutDistF = species[self.neutNm]:getDistF()
      self.m0elc:copy(elcM0)
     
      self.confMult:advance(tCurr, {coefIz, self.m0elc}, {self.coefM0})
      self.confPhaseMult:advance(tCurr, {self.coefM0, neutDistF}, {self.ionizSrc})

      fRhsOut:accumulate(-1.0,self.ionizSrc)  
   else
      -- Ions. 
      self.m0elc:copy(elcM0)
      local neutM0   = species[self.neutNm]:fluidMoments()[1]
      local neutU    = species[self.neutNm]:selfPrimitiveMoments()[1] 
      local neutVtSq = species[self.neutNm]:selfPrimitiveMoments()[2]
      -- Include only z-component of neutU.

      self.confMult:advance(tCurr, {coefIz, self.m0elc}, {self.coefM0})
      species[self.speciesName].calcMaxwell:advance(tCurr,
         {neutM0, neutU, neutVtSq, species[self.speciesName].bmag}, {self.fMaxNeut})
      self.confPhaseMult:advance(tCurr, {self.coefM0, self.fMaxNeut}, {self.ionizSrc})

      fRhsOut:accumulate(1.0,self.ionizSrc)
   end
   self.timers.nonSlvr = self.timers.nonSlvr + Time.clock() - tmNonSlvrStart
end
   
function GkIonization:write(tm, frame)
end

function GkIonization:setCfl(cfl)
   self.cfl = cfl
end

function GkIonization:getIonizSrc()
   return self.ionizSrc
end

function GkIonization:slvrTime()
   return 0.
end

function GkIonization:nonSlvrTime()
   return self.timers.nonSlvr
end

function GkIonization:projectMaxwellTime()
   local time = self.confMult.projectMaxwellTime()
   time = time + self.confPhaseMult.projectMaxwellTime
   return time
end

return GkIonization

