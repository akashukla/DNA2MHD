-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: BGK Collision operator or Vlasov-Maxwell.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase = require "App.Collisions.CollisionsBase"
local Constants      = require "Lib.Constants"
local DataStruct     = require "DataStruct"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"
local Mpi            = require "Comm.Mpi"
local lume           = require "Lib.lume"
local xsys           = require "xsys"


-- VmBGKCollisions ---------------------------------------------------------------
--
-- Bhatnagar-Gross-Krook Collision operator.
--------------------------------------------------------------------------------

local VmBGKCollisions = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmBGKCollisions:init(tbl)
   self.tbl = tbl
end

-- Function to find the index of an element in table.
local function findInd(tbl, el)
   for i, v in ipairs(tbl) do
      if v == el then
         return i
      end
   end
   return #tbl+1    -- If not found return a number larger than the length of the table.
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function VmBGKCollisions:fullInit(speciesTbl)
   local tbl = self.tbl  -- Previously stored table.

   self.collKind = "VmBGK"  -- Type of collisions model. Useful at the species app level.

   self.cfl = 0.1  -- Some default value.

   self.collidingSpecies = assert(
      tbl.collideWith,
      "App.VmLBOCollisions: Must specify names of species to collide with in 'collideWith'.")

   -- First determine if self-species and/or cross-species collisions take place,
   -- and (if cross-collisions=true) put the names of the other colliding species in a list.
   local selfSpecInd = findInd(self.collidingSpecies, self.speciesName)
   if selfSpecInd < (#self.collidingSpecies+1) then
      self.selfCollisions = true  -- Apply self-species collisions.
      if #self.collidingSpecies > 1 then
         self.crossCollisions = true  -- Apply cross-species collisions.
         self.crossSpecies = lume.clone(self.collidingSpecies)
         table.remove(self.crossSpecies, selfSpecInd)
      else
         self.crossCollisions = false  -- Don't apply cross-species collisions.
      end
   else
      self.selfCollisions = false  -- Don't apply self-species collisions.
      self.crossCollisions = true  -- Apply cross-species collisions.
      self.crossSpecies = lume.clone(self.collidingSpecies)  -- All species in collidingSpecies must be cross-species.
   end

   -- Now establish if user wants constant or spatially varying collisionality.
   -- For constant nu, separate self and cross collision frequencies.
   self.collFreqs = tbl.frequencies -- List of collision frequencies, if using spatially constant nu.
   if self.collFreqs then
      self.varNu = false  -- Not spatially varying nu.
      self.cellConstNu = true  -- Cell-wise constant nu?
      if self.selfCollisions then
         self.collFreqSelf = self.collFreqs[selfSpecInd]
      end
      if self.crossCollisions then
         self.collFreqCross = lume.clone(self.collFreqs)
         table.remove(self.collFreqCross, selfSpecInd)
      end
   else
      self.varNu  = true                 -- Spatially varying nu.
      self.mass   = speciesTbl.mass      -- Mass of this species.
      self.charge = speciesTbl.charge    -- Charge of this species.
      -- For now only cell-wise constant nu is implemented.
      -- self.cellConstNu = assert(tbl.cellAvFrequencies, "App.GkLBOCollisions: Must specify 'useCellAverageNu=true/false' for using cellwise constant/expanded spatially varying collisionality.")
      self.cellConstNu = true
      -- If no constant collision frequencies provided ('frequencies'), user can specify 'normNu'
      -- list of collisionalities normalized by (T_0^(3/2)/n_0) evaluated somewhere in the
      -- simulation. Otherwise code compute Spitzer collisionality from scratch.
      self.normNuIn = tbl.normNu
      -- normNuSelf, epsilon0 and elemCharge may not used, but are
      -- initialized below to avoid if-statements in advance method.
      if self.normNuIn then
         self.userInputNormNu = true
         if self.selfCollisions then
            self.normNuSelf = self.normNuIn[selfSpecInd]
         end
         if self.crossCollisions then
            self.normNuCross = lume.clone(self.normNuIn)
            table.remove(self.normNuCross, selfSpecInd)
         end
         self.epsilon0   = Constants.EPSILON0
         self.elemCharge = Constants.ELEMENTARY_CHARGE
      else
         self.userInputNormNu = false
         if self.selfCollisions then
            self.normNuSelf = 0.0
         end
         if self.crossCollisions then
            self.normNuCross = lume.clone(self.collidingSpecies)
            table.remove(self.normNuCross, selfSpecInd)
            for i, _ in ipairs(self.normNuCross) do self.normNuCross[i] = 0.0 end
         end
         self.epsilon0 = assert(
            tbl.epsilon0, "Updater.VmLBOCollisions: Must specify vacuum permittivity 'epsilon0' ('frequencies' and 'normNu' are not specified, so nu is calculated via Spitzer).")
         self.elemCharge = assert(
            tbl.elemCharge, "Updater.VmLBOCollisions: Must specify elementary charge with 'elemCharge' ('frequencies' and 'normNu' are not specified, so nu is calculated via Spitzer).")
--         self.coulombLog = assert(
--      	    tbl.coulombLog, "Updater.VmBGKCollisions: Must specify Coulomb logaritm with 'coulombLog' ('frequencies' and 'normNu' are not specified, so nu is calculated via Spitzer).")
      end
   end

   if self.crossCollisions then
      self.mass   = speciesTbl.mass          -- Mass of this species.
      self.charge = speciesTbl.charge        -- Charge of this species.
      local betaGreeneIn = tbl.betaGreene    -- Can specify 'betaGreene' free parameter in Grene cross-species collisions.
      if betaGreeneIn then
        self.betaGreene = betaGreeneIn
      else
        self.betaGreene = 0.0                -- Default value.
      end
   end

   if tbl.nuFrac then
      self.nuFrac = tbl.nuFrac
   else
      self.nuFrac = 1.0
   end

   self.exactLagFixM012 = xsys.pickBool(tbl.exactLagFixM012, true) 

   self.tmEvalMom = 0.0
end

function VmBGKCollisions:setName(nm)
   self.name = nm
end
function VmBGKCollisions:setSpeciesName(nm)
   self.speciesName = nm
end

function VmBGKCollisions:setCfl(cfl)
   self.cfl = cfl
end
function VmBGKCollisions:setConfBasis(basis)
   self.confBasis = basis
end
function VmBGKCollisions:setConfGrid(grid)
   self.confGrid = grid
end
function VmBGKCollisions:setPhaseBasis(basis)
   self.phaseBasis = basis
end
function VmBGKCollisions:setPhaseGrid(grid)
   self.phaseGrid = grid
end

function VmBGKCollisions:createSolver()
   self.numVelDims = self.phaseGrid:ndim() - self.confGrid:ndim()

   local function createConfFieldCompV()
      return DataStruct.Field {
	 onGrid        = self.confGrid,
	 numComponents = self.confBasis:numBasis()*self.numVelDims,
	 ghost         = {1, 1},
      }
   end
   local function createConfFieldComp1()
      return DataStruct.Field {
	 onGrid        = self.confGrid,
	 numComponents = self.confBasis:numBasis(),
	 ghost         = {1, 1},
      }
   end

   if self.exactLagFixM012 then
      -- Intermediate moments used in Lagrange fixing.
      self.dM1 = createConfFieldCompV()
      self.dM2 = createConfFieldComp1()
      self.dM0 = createConfFieldComp1()
      -- Create updater to compute M0, M1i, M2 moments sequentially.
      self.fiveMomentsCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.phaseGrid,
         phaseBasis = self.phaseBasis,
         confBasis  = self.confBasis,
         moment     = "FiveMoments",
      }
      self.lagFix = Updater.LagrangeFix {
         onGrid     = self.phaseGrid,
         phaseGrid  = self.phaseGrid,
         phaseBasis = self.phaseBasis,
         confGrid   = self.confGrid,
         confBasis  = self.confBasis,
         mode       = 'Vlasov',
      }
   end

   -- Sum of Maxwellians multiplied by respective collisionalities.
   self.nufMaxwellSum = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
   }

   if self.varNu then
      -- Collisionality, nu, summed over all species pairs
      self.nuSum = createConfFieldComp1()
      -- Updater to compute spatially varying (Spitzer) nu
      self.spitzerNu = Updater.SpitzerCollisionality {
         onGrid = self.confGrid,
         confBasis        = self.confBasis,
         useCellAverageNu = self.cellConstNu,
         willInputNormNu  = self.userInputNormNu,
         elemCharge       = self.elemCharge,
         epsilon0         = self.epsilon0,
      }
      -- Weak multiplication to multiply nu(x) with fMaxwell
      self.phaseMul = Updater.CartFieldBinOp {
         onGrid     = self.phaseGrid,
         weakBasis  = self.phaseBasis,
         fieldBasis = self.confBasis,
         operation  = "Multiply",
      }
   else
      self.nuSum = 0.0  -- Assigned in advance method
   end

   if self.crossCollisions then
      -- Cross-collision Maxwellian multiplied by collisionality
      self.nufMaxwellCross = DataStruct.Field {
         onGrid        = self.phaseGrid,
         numComponents = self.phaseBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Dummy fields for the primitive moment calculator
      self.uCrossSq = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      if self.varNu then
         -- Temporary collisionality field
         self.nuCrossSelf  = createConfFieldComp1()
         self.nuCrossOther = createConfFieldComp1()
      else
         self.nuCrossSelf  = 0.0
         self.nuCrossOther = 0.0
      end
      -- Updater to compute cross-species primitive moments.
      self.primMomCross = Updater.CrossPrimMoments {
         onGrid     = self.confGrid,
         phaseBasis = self.phaseBasis,
         confBasis  = self.confBasis,
         operator   = "VmBGK",
         betaGreene = self.betaGreene,
      }
      if self.exactLagFixM012 then
         -- Will need the 1st and 2nd moment of the cross-species Maxwellian.
         self.crossMaxwellianM1 = createConfFieldCompV()
         self.crossMaxwellianM2 = createConfFieldComp1()
         self.confMultiply = Updater.CartFieldBinOp {
            onGrid    = self.confGrid,
            weakBasis = self.confBasis,
            operation = "Multiply",
         }
         self.confDotProduct = Updater.CartFieldBinOp {
            onGrid    = self.confGrid,
            weakBasis = self.confBasis,
            operation = "DotProduct",
         }
      end
   end

   -- Maxwellian solver
   self.maxwellian = Updater.MaxwellianOnBasis {
      onGrid     = self.phaseGrid,
      confGrid   = self.confGrid,
      confBasis  = self.confBasis,
      phaseGrid  = self.phaseGrid,
      phaseBasis = self.phaseBasis,
   }
   -- BGK Collision solver itself
   self.collisionSlvr = Updater.BgkCollisions {
      onGrid           = self.phaseGrid,
      confGrid         = self.confGrid,
      confBasis        = self.confBasis,
      phaseBasis       = self.phaseBasis,
      varyingNu        = self.varNu,
      useCellAverageNu = self.cellConstNu,
   }

end

function VmBGKCollisions:advance(tCurr, fIn, species, fRhsOut)

   -- Fetch coupling moments and primitive moments of this species
   local selfMom     = species[self.speciesName]:fluidMoments()
   local primMomSelf = species[self.speciesName]:selfPrimitiveMoments()

   if self.varNu then
      self.nuSum:clear(0.0)
   else
      self.nuSum = 0.0
   end
   self.nufMaxwellSum:clear(0.0)

   if self.selfCollisions then
      self.maxwellian:advance(tCurr, {selfMom[1], primMomSelf[1], primMomSelf[2]},
                              {self.nufMaxwellSum})
      if self.exactLagFixM012 then
         self.fiveMomentsCalc:advance(tCurr, {self.nufMaxwellSum}, { self.dM0, self.dM1, self.dM2 })
         -- Barrier before manipulations to moments before passing them to Lagrange Fix updater
         Mpi.Barrier(self.phaseGrid:commSet().sharedComm)
	 self.dM0:scale(-1)
	 self.dM0:accumulate(1, selfMom[1])
	 self.dM1:scale(-1)
	 self.dM1:accumulate(1, selfMom[2])
	 self.dM2:scale(-1)
	 self.dM2:accumulate(1, selfMom[3])
	 self.lagFix:advance(tCurr, {self.dM0, self.dM1, self.dM2}, {self.nufMaxwellSum})
      end

      if self.varNu then
         -- Compute the collisionality
         self.spitzerNu:advance(tCurr, {self.mass, self.charge, selfMom[1], primMomSelf[2], self.normNuSelf},
				{self.nuSum})
         self.phaseMul:advance(tCurr, {self.nuSum, self.nufMaxwellSum},
			       {self.nufMaxwellSum})
      else
         self.nuSum = self.collFreqSelf
         -- Barrier before scaling by the collisionality
         -- Note that this is not necessary when nu varies because the scaling is done via advance methods
         Mpi.Barrier(self.phaseGrid:commSet().sharedComm)
         self.nufMaxwellSum:scale(self.collFreqSelf)
      end
   end    -- end if self.selfCollisions

   if self.crossSpecies then
      for sInd, otherNm in ipairs(self.crossSpecies) do
	 local mOther       = species[otherNm]:getMass()
	 local otherMom     = species[otherNm]:fluidMoments()
         local primMomOther = species[otherNm]:selfPrimitiveMoments()

         local tmEvalMomStart = Time.clock()

	 -- Collision frequency established before computing
         -- crossPrimMom in case we want to generalize Greene without
         -- m_s*n_s*nu_sr=m_r*n_r*nu_rs.
         if self.varNu then
            -- Compute the collisionality if another species hasn't
            -- already done so.
            if (not species[self.speciesName].momentFlags[6][otherNm]) then
               self.spitzerNu:advance(tCurr, {self.mass, self.charge, otherMom[1], primMomSelf[2], self.normNuCross[sInd]},
                                      {species[self.speciesName].nuVarXCross[otherNm]})
               species[self.speciesName].momentFlags[6][otherNm] = true
            end
            if (not species[otherNm].momentFlags[6][self.speciesName]) then
               local chargeOther = species[otherNm]:getCharge()
               self.spitzerNu:advance(tCurr, {mOther, chargeOther, selfMom[1], primMomOther[2],
                                              species[otherNm].collPairs[otherNm][self.speciesName].normNu},
                                      {species[otherNm].nuVarXCross[self.speciesName]})
               species[otherNm].momentFlags[6][self.speciesName] = true
            end
            self.nuCrossSelf:copy(species[self.speciesName].nuVarXCross[otherNm])
            self.nuCrossOther:copy(species[otherNm].nuVarXCross[self.speciesName])
         else
            self.nuCrossSelf  = self.collFreqCross[sInd]
            self.nuCrossOther = species[otherNm].collPairs[otherNm][self.speciesName].nu
         end

         if (not (species[self.speciesName].momentFlags[5][otherNm] and
                  species[otherNm].momentFlags[5][self.speciesName])) then
            -- Cross-primitive moments for the collision of these two species has not been computed.
            self.primMomCross:advance(tCurr, {self.mass, self.nuCrossSelf, selfMom, primMomSelf,
                                              mOther, self.nuCrossOther, otherMom, primMomOther},
                                      {species[self.speciesName].uCross[otherNm], species[self.speciesName].vtSqCross[otherNm],
                                       species[otherNm].uCross[self.speciesName], species[otherNm].vtSqCross[self.speciesName]})

            species[self.speciesName].momentFlags[5][otherNm] = true
            species[otherNm].momentFlags[5][self.speciesName] = true
         end

         self.tmEvalMom = self.tmEvalMom + Time.clock() - tmEvalMomStart

	 self.maxwellian:advance(tCurr, {selfMom[1], species[self.speciesName].uCross[otherNm],
                                         species[self.speciesName].vtSqCross[otherNm]}, {self.nufMaxwellCross})
         if self.exactLagFixM012 then
            -- Need to compute the first and second moment of the cross-species Maxwellian (needed by Lagrange fix).
            self.confMultiply:advance(tCurr, {species[self.speciesName].uCross[otherNm], selfMom[1]}, {self.crossMaxwellianM1})
            self.confDotProduct:advance(tCurr, {species[self.speciesName].uCross[otherNm],
                                                species[self.speciesName].uCross[otherNm]}, {self.crossMaxwellianM2})
            self.crossMaxwellianM2:accumulate(1.0, species[self.speciesName].vtSqCross[otherNm])
            self.confMultiply:advance(tCurr, {selfMom[1],self.crossMaxwellianM2}, {self.crossMaxwellianM2})
            -- Now compute current moments of Maxwellians and perform Lagrange fix. 
            self.fiveMomentsCalc:advance(tCurr, {self.nufMaxwellCross}, { self.dM0, self.dM1, self.dM2 })
            -- Barrier before manipulations to moments before passing them to Lagrange Fix updater.
            Mpi.Barrier(self.phaseGrid:commSet().sharedComm)
            self.dM0:scale(-1)
            self.dM0:accumulate(1, selfMom[1])
            self.dM1:scale(-1)
            self.dM1:accumulate(1, self.crossMaxwellianM1)
            self.dM2:scale(-1)
            self.dM2:accumulate(1, self.crossMaxwellianM2)
            self.lagFix:advance(tCurr, {self.dM0, self.dM1, self.dM2}, {self.nufMaxwellSum})
         end

         if self.varNu then
            self.phaseMul:advance(tCurr, {self.nuCrossSelf, self.nufMaxwellCross}, {self.nufMaxwellCross})

            -- Barrier over shared communicator before accumulate
            Mpi.Barrier(self.phaseGrid:commSet().sharedComm)

            self.nuSum:accumulate(1.0, self.nuCrossSelf)
            self.nufMaxwellSum:accumulate(1.0, self.nufMaxwellCross)
         else
            self.nuSum = self.nuSum+self.nuCrossSelf

            -- Barrier over shared communicator before accumulate
            Mpi.Barrier(self.phaseGrid:commSet().sharedComm)

            self.nufMaxwellSum:accumulate(self.nuCrossSelf, self.nufMaxwellCross)
         end
      end    -- end loop over other species that this species collides with
   end    -- end if self.crossCollisions

   self.collisionSlvr:advance(tCurr, {fIn, self.nufMaxwellSum, self.nuSum}, {fRhsOut})

end

function VmBGKCollisions:write(tm, frame)
end

function VmBGKCollisions:totalTime()
   return self.collisionSlvr.totalTime + self.maxwellian.totalTime + self.tmEvalMom
end

function VmBGKCollisions:slvrTime()
   return self.collisionSlvr.totalTime + self.maxwellian.totalTime
end

function VmBGKCollisions:momTime()
   return self.tmEvalMom
end

return VmBGKCollisions