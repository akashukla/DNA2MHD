-- Gkyl ------------------------------------------------------------------------
--
-- Test for cartesian grid objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Grid = require "Grid"
local Lin  = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats        = Unit.stats

function test_1()
   local grid = Grid.RectCart {
      cells = {10, 20}
   }

   -- Just make sure setIndex() method works. For RectCart object
   -- setting index is not needed.
   idx = Lin.IntVec(grid:ndim())
   idx[1], idx[2] = 1, 1
   grid:setIndex(idx)
   
   assert_equal(2, grid:ndim(), "Checking NDIM")
   assert_equal("uniform", grid:id(), "Checking ID")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")

   local localRange = grid:localRange()
   assert_equal(1, localRange:lower(1), "Checking region bounds")
   assert_equal(10, localRange:upper(1), "Checking region bounds")

   assert_equal(1, localRange:lower(2), "Checking region bounds")
   assert_equal(20, localRange:upper(2), "Checking region bounds")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(0.0, grid:lower(2), "Checking lower")

   assert_equal(1.0, grid:upper(1), "Checking upper")
   assert_equal(1.0, grid:upper(2), "Checking upper")

   assert_equal(0.1, grid:dx(1), "Checking dx")
   assert_equal(0.05, grid:dx(2), "Checking dx")

   -- Test cell-center coordinates.
   local lox, dx = grid:lower(1), grid:dx(1)
   local loy, dy = grid:lower(2), grid:dx(2)

   
   local xc = Lin.Vec(2)
   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 grid:setIndex( idx:setValues {i, j} )
	 grid:cellCenter(xc)

	 assert_equal(xc[1], lox+(i-0.5)*dx, "Testing cell-center coordinate")
	 assert_equal(xc[2], loy+(j-0.5)*dy, "Testing cell-center coordinate")
      end
   end
   
end

function test_2()
   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {2.0, 5.0},
      cells = {10, 20}
   }

   assert_equal(2, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")

   assert_equal(0.2*0.2, grid:cellVolume(), "Checking volume")

   local localRange = grid:localRange()
   -- Test cell-center coordinates.
   local lox, dx = grid:lower(1), grid:dx(1)
   local loy, dy = grid:lower(2), grid:dx(2)

   local idx = Lin.IntVec(grid:ndim())
   local xc = Lin.Vec(2)
   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 grid:setIndex( idx:setValues {i, j} )
	 grid:cellCenter(xc)

	 assert_equal(xc[1], lox+(i-0.5)*dx, "Testing cell-center coordinate")
	 assert_equal(xc[2], loy+(j-0.5)*dy, "Testing cell-center coordinate")
      end
   end
end

function test_3()
   local grid = Grid.RectCart {
      lower = {0.0, 1.0, 2.0},
      upper = {2.0, 5.0, 10.0},
      cells = {10, 20, 40}
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")
   assert_equal(40, grid:numCells(3), "Checking numCells")

   local localRange = grid:localRange()
   assert_equal(1, localRange:lower(1), "Checking region bounds")
   assert_equal(10, localRange:upper(1), "Checking region bounds")

   assert_equal(1, localRange:lower(2), "Checking region bounds")
   assert_equal(20, localRange:upper(2), "Checking region bounds")

   assert_equal(1, localRange:lower(3), "Checking region bounds")
   assert_equal(40, localRange:upper(3), "Checking region bounds")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")
   assert_equal(2.0, grid:lower(3), "Checking lower")   

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")
   assert_equal(10.0, grid:upper(3), "Checking upper")   

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")
   assert_equal(0.2, grid:dx(3), "Checking dx")

   assert_equal(0.2*0.2*0.2, grid:cellVolume(), "Checking volume")

   -- Test cell-center coordinates.
   local lox, dx = grid:lower(1), grid:dx(1)
   local loy, dy = grid:lower(2), grid:dx(2)
   local loz, dz = grid:lower(3), grid:dx(3)

   local idx = Lin.IntVec(grid:ndim())
   local xc = Lin.Vec(3)
   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 for k = localRange:lower(3), localRange:upper(3) do
	    grid:setIndex( idx:setValues {i, j, k} )
	    grid:cellCenter(xc)
	    
	    assert_equal(xc[1], lox+(i-0.5)*dx, "Testing cell-center coordinate")
	    assert_equal(xc[2], loy+(j-0.5)*dy, "Testing cell-center coordinate")
	    assert_equal(xc[3], loz+(k-0.5)*dz, "Testing cell-center coordinate")
	 end
      end
   end
end

function test_4()
   local grid = Grid.NonUniformRectCart {
      cells = {10, 10}
   }

   assert_equal(2, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(10, grid:numCells(2), "Checking numCells")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(0.0, grid:lower(2), "Checking lower")

   assert_equal(1.0, grid:upper(1), "Checking upper")
   assert_equal(1.0, grid:upper(2), "Checking upper")

   assert_equal(0.1, grid:dx(1), "Checking dx 1")
   assert_equal(0.1, grid:dx(2), "Checking dx 2")

   assert_equal(0.1*0.1, grid:cellVolume(), "Checking volume")

end

function test_5a()
   local grid = Grid.NonUniformRectCart {   
      lower = {0.0, 1.0, 2.0},
      upper = {2.0, 5.0, 10.0},
      cells = {10, 20, 40},
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")
   assert_equal(40, grid:numCells(3), "Checking numCells")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")
   assert_equal(2.0, grid:lower(3), "Checking lower")   

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")
   assert_equal(10.0, grid:upper(3), "Checking upper")   

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")
   assert_equal(0.2, grid:dx(3), "Checking dx")

   assert_equal(0.2*0.2*0.2, grid:cellVolume(), "Checking volume")

   -- Test cell-center coordinates (THIS WORKS AS MESH IS ACTUALLY UNIFORM).
   local lox, dx = grid:lower(1), grid:dx(1)
   local loy, dy = grid:lower(2), grid:dx(2)
   local loz, dz = grid:lower(3), grid:dx(3)

   local idx        = Lin.IntVec(grid:ndim())
   local localRange = grid:localRange()
   local xc         = Lin.Vec(3)
   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 for k = localRange:lower(3), localRange:upper(3) do
            grid:setIndex( idx:setValues {i, j, k} )
            grid:cellCenter(xc)
        
            assert_equal(lox+(i-0.5)*dx, xc[1], "Testing cell-center x coordinate")
            assert_equal(loy+(j-0.5)*dy, xc[2], "Testing cell-center y coordinate")
            assert_equal(loz+(k-0.5)*dz, xc[3], "Testing cell-center z coordinate")
        
            local xcDir = {grid:cellCenterInDir(1), grid:cellCenterInDir(2), grid:cellCenterInDir(3)} 
            assert_equal(xc[1], xcDir[1], "Testing cell center in x direction")
            assert_equal(xc[2], xcDir[2], "Testing cell center in y direction")
            assert_equal(xc[3], xcDir[3], "Testing cell center in z direction")
            
            local loDir = {grid:cellLowerInDir(1), grid:cellLowerInDir(2), grid:cellLowerInDir(3)} 
            assert_equal(lox+(i-1)*dx, loDir[1], "Testing cell lower in x direction")
            assert_equal(loy+(j-1)*dy, loDir[2], "Testing cell lower in y direction")
            assert_equal(loz+(k-1)*dz, loDir[3], "Testing cell lower in z direction")
            
            local upDir = {grid:cellUpperInDir(1), grid:cellUpperInDir(2), grid:cellUpperInDir(3)} 
            assert_equal(lox+i*dx, upDir[1], "Testing cell upper in x direction")
            assert_equal(loy+j*dy, upDir[2], "Testing cell upper in y direction")
            assert_equal(loz+k*dz, upDir[3], "Testing cell upper in z direction")
	 end
      end
   end
end

function test_5()
   local grid = Grid.NonUniformRectCart {   
      lower    = {0.0, 1.0, 2.0},
      upper    = {2.0, 5.0, 10.0},
      cells    = {10, 20, 40},
      -- Functions mapping computational space to physical space.
      mappings = {
	 function (zeta)
	    return zeta
	 end,
	 function (zeta)
	    return zeta	    
	 end,
	 function (zeta)
	    return zeta
	 end,
      }
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")
   assert_equal(40, grid:numCells(3), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:lower(2), "Checking lower")
   assert_equal(2.0, grid:lower(3), "Checking lower")   

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(5.0, grid:upper(2), "Checking upper")
   assert_equal(10.0, grid:upper(3), "Checking upper")   

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.2, grid:dx(2), "Checking dx")
   assert_equal(0.2, grid:dx(3), "Checking dx")

   assert_equal(0.2*0.2*0.2, grid:cellVolume(), "Checking volume")

   -- Test cell-center coordinates (THIS WORKS AS MESH IS ACTUALLY UNIFORM).
   local lox, dx = grid:lower(1), grid:dx(1)
   local loy, dy = grid:lower(2), grid:dx(2)
   local loz, dz = grid:lower(3), grid:dx(3)

   local idx        = Lin.IntVec(grid:ndim())
   local localRange = grid:localRange()
   local xc         = Lin.Vec(3)
   for i = localRange:lower(1), localRange:upper(1) do
      for j = localRange:lower(2), localRange:upper(2) do
	 for k = localRange:lower(3), localRange:upper(3) do
	    grid:setIndex( idx:setValues {i, j, k} )
	    grid:cellCenter(xc)

	    assert_equal(lox+(i-0.5)*dx, xc[1], "Testing cell-center x coordinate")
	    assert_equal(loy+(j-0.5)*dy, xc[2], "Testing cell-center y coordinate")
	    assert_equal(loz+(k-0.5)*dz, xc[3], "Testing cell-center z coordinate")

            local xcDir = {grid:cellCenterInDir(1), grid:cellCenterInDir(2), grid:cellCenterInDir(3)} 
	    assert_equal(xc[1], xcDir[1], "Testing cell center in x direction")
	    assert_equal(xc[2], xcDir[2], "Testing cell center in y direction")
	    assert_equal(xc[3], xcDir[3], "Testing cell center in z direction")
            
            local loDir = {grid:cellLowerInDir(1), grid:cellLowerInDir(2), grid:cellLowerInDir(3)} 
            assert_equal(lox+(i-1)*dx, loDir[1], "Testing cell lower in x direction")
            assert_equal(loy+(j-1)*dy, loDir[2], "Testing cell lower in y direction")
            assert_equal(loz+(k-1)*dz, loDir[3], "Testing cell lower in z direction")
            
            local upDir = {grid:cellUpperInDir(1), grid:cellUpperInDir(2), grid:cellUpperInDir(3)} 
            assert_equal(lox+i*dx, upDir[1], "Testing cell upper in x direction")
            assert_equal(loy+j*dy, upDir[2], "Testing cell upper in y direction")
            assert_equal(loz+k*dz, upDir[3], "Testing cell upper in z direction")
	 end
      end
   end   
end

function test_6()
   local grid = Grid.NonUniformRectCart {   
      cells    = {10, 20, 40},
      -- Functions mapping computational space to physical space.
      mappings = {
	 function (zeta)
	    return 2*zeta
	 end,
	 function (zeta)
	    return 2*zeta	    
	 end,
	 function (zeta)
	    return 2*zeta
	 end,
      }
   }

   assert_equal(3, grid:ndim(), "Checking NDIM")

   assert_equal(10, grid:numCells(1), "Checking numCells")
   assert_equal(20, grid:numCells(2), "Checking numCells")
   assert_equal(40, grid:numCells(3), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(0.0, grid:lower(2), "Checking lower")
   assert_equal(0.0, grid:lower(3), "Checking lower")   

   assert_equal(2.0, grid:upper(1), "Checking upper")
   assert_equal(2.0, grid:upper(2), "Checking upper")
   assert_equal(2.0, grid:upper(3), "Checking upper")   

   assert_equal(0.2, grid:dx(1), "Checking dx")
   assert_equal(0.1, grid:dx(2), "Checking dx")
   assert_equal(0.05, grid:dx(3), "Checking dx")

   assert_equal(0.2*0.1*0.05, grid:cellVolume(), "Checking volume")
end

function test_7()
   local grid = Grid.NonUniformRectCart {   
      cells    = {3},
      -- Functions mapping computational space to physical space.
      mappings = {
	 function (zeta)
	    return zeta*zeta
	 end,
      }
   }

   local idx = Lin.IntVec(grid:ndim())
   
   assert_equal(1, grid:ndim(), "Checking NDIM")

   assert_equal(3, grid:numCells(1), "Checking numCells")

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:upper(1), "Checking upper")

   grid:setIndex( idx:setValues {1} )
   assert_equal(1/9, grid:dx(1), "Checking dx")
   assert_equal(1/9, grid:cellVolume(), "Checking volume")

   grid:setIndex( idx:setValues {2} )
   assert_equal(1/3, grid:dx(1), "Checking dx")
   assert_equal(1/3, grid:cellVolume(), "Checking volume")

   grid:setIndex(idx:setValues {3} )
   assert_equal(5/9, grid:dx(1), "Checking dx")
   assert_equal(5/9, grid:cellVolume(), "Checking volume")
end

function test_8()
   local grid = Grid.NonUniformRectCart { cells = {3} }
   local xn   = grid:nodeCoords(1)
   -- Set nodes manually (this is the mapping zeta^2).
   xn[1] = 0.0
   xn[2] = 1/3*1/3
   xn[3] = 2/3*2/3
   xn[4] = 1*1
   -- Done.

   local idx = Lin.IntVec(grid:ndim())

   assert_equal(1, grid:ndim(), "Checking NDIM")

   assert_equal(3, grid:numCells(1), "Checking numCells")   

   assert_equal(0.0, grid:lower(1), "Checking lower")
   assert_equal(1.0, grid:upper(1), "Checking upper")

   grid:setIndex( idx:setValues {1} )
   assert_equal(1/9, grid:dx(1), "Checking dx")
   assert_equal(1/9, grid:cellVolume(), "Checking volume")

   grid:setIndex( idx:setValues {2} )
   assert_equal(1/3, grid:dx(1), "Checking dx")
   assert_equal(1/3, grid:cellVolume(), "Checking volume")

   grid:setIndex( idx:setValues {3} )
   assert_equal(5/9, grid:dx(1), "Checking dx")
   assert_equal(5/9, grid:cellVolume(), "Checking volume")   
end

function test_9()
   local grid = Grid.RectCart {
      lower        = {0.0, 1.0, 1.0},
      upper        = {2.0, 5.0, 10.0},
      cells        = {10, 20, 30},
      periodicDirs = {1, 3},
   }

   -- Check periodicity.
   assert_equal(true, grid:isDirPeriodic(1), "Checking periodicity")
   assert_equal(false, grid:isDirPeriodic(2), "Checking periodicity")
   assert_equal(true, grid:isDirPeriodic(3), "Checking periodicity")
end

function test_10()
   local grid = Grid.NonUniformRectCart {
      lower        = {0.0, 1.0, 1.0},
      upper        = {2.0, 5.0, 10.0},
      cells        = {10, 20, 30},
      periodicDirs = {2, 3},
   }

   -- Check periodicity.
   assert_equal(false, grid:isDirPeriodic(1), "Checking periodicity (NU)")
   assert_equal(true, grid:isDirPeriodic(2), "Checking periodicity (NU)")
   assert_equal(true, grid:isDirPeriodic(3), "Checking periodicity (NU)")
end

function test_12()
   -- Test the grid's findCell method.
   local grid = Grid.RectCart {
      lower = {0.0, 1.0},
      upper = {2.0, 5.0},
      cells = {10, 20}
   }

   -- Find a point inside a cell.
   local fIdx, np = {0, 0}, {0.35, 2.1}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking fIdx[1] for interior point")
   assert_equal(6, fIdx[2], "Checking fIdx[2] for interior point")

   -- Find a point on boundary of two cells.
   local fIdx, np = {0, 0}, {0.35, 2.2}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking fIdx[1] for cell boundary point (lower)")
   assert_equal(6, fIdx[2], "Checking fIdx[2] for cell boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(2, fIdx[1], "Checking fIdx[1] for cell boundary point (upper)")
   assert_equal(7, fIdx[2], "Checking fIdx[2] for cell boundary point (upper)")

   -- Find a point on corner of 4 cells.
   local fIdx, np = {0, 0}, {0.4, 3.0}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking fIdx[1] for corner point (lower)")
   assert_equal(10, fIdx[2], "Checking fIdx[2] for corner point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking fIdx[1] for corner point (upper)")
   assert_equal(11, fIdx[2], "Checking fIdx[2] for corner point (upper)")

   -- Find a point on lower x domain boundary.
   local fIdx, np = {0, 0}, {0.0, 3.1}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking fIdx[1] for lower-x point")
   assert_equal(11, fIdx[2], "Checking fIdx[2] for lower-x point")
   local fIdx, np = {0, 0}, {0.0, 3.}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking fIdx[1] for lower-x boundary point (lower)")
   assert_equal(10, fIdx[2], "Checking fIdx[2] for lower-x boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(1, fIdx[1], "Checking fIdx[1] for lower-x boundary point (upper)")
   assert_equal(11, fIdx[2], "Checking fIdx[2] for lower-x boundary point (upper)")

   -- Find a point on lower y domain boundary.
   local fIdx, np = {0, 0}, {0.5, 1.}
   grid:findCell(np, fIdx)
   assert_equal(3, fIdx[1], "Checking fIdx[1] for lower-y point")
   assert_equal(1, fIdx[2], "Checking fIdx[2] for lower-y point")
   local fIdx, np = {0, 0}, {0.4, 1.}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking fIdx[1] for lower-y boundary point (lower)")
   assert_equal(1, fIdx[2], "Checking fIdx[2] for lower-y boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking fIdx[1] for lower-y boundary point (upper)")
   assert_equal(1, fIdx[2], "Checking fIdx[2] for lower-y boundary point (upper)")

   -- Find a point on upper x domain boundary.
   local fIdx, np = {0, 0}, {2.0, 3.1}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking fIdx[1] for upper-x point")
   assert_equal(11, fIdx[2], "Checking fIdx[2] for upper-x point")
   local fIdx, np = {0, 0}, {2.0, 3.}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking fIdx[1] for upper-x boundary point (lower)")
   assert_equal(10, fIdx[2], "Checking fIdx[2] for upper-x boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(10, fIdx[1], "Checking fIdx[1] for upper-x boundary point (upper)")
   assert_equal(11, fIdx[2], "Checking fIdx[2] for upper-x boundary point (upper)")

   -- Find a point on upper y domain boundary.
   local fIdx, np = {0, 0}, {0.5, 5.}
   grid:findCell(np, fIdx)
   assert_equal(3, fIdx[1], "Checking fIdx[1] for upper-y point")
   assert_equal(20, fIdx[2], "Checking fIdx[2] for upper-y point")
   local fIdx, np = {0, 0}, {0.4, 5.}
   grid:findCell(np, fIdx)
   assert_equal(2, fIdx[1], "Checking fIdx[1] for upper-y boundary point (lower)")
   assert_equal(20, fIdx[2], "Checking fIdx[2] for upper-y boundary point (lower)")
   grid:findCell(np, fIdx, false)
   assert_equal(3, fIdx[1], "Checking fIdx[1] for upper-y boundary point (upper)")
   assert_equal(20, fIdx[2], "Checking fIdx[2] for upper-y boundary point (upper)")

   -- Find domain corner points.
   local fIdx, np = {0, 0}, {0.0, 1.}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking fIdx[1] for lower left corner point")
   assert_equal(1, fIdx[2], "Checking fIdx[2] for lower left corner point")
   local fIdx, np = {0, 0}, {2.0, 1.}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking fIdx[1] for lower right corner point")
   assert_equal(1, fIdx[2], "Checking fIdx[2] for lower right corner point")
   local fIdx, np = {0, 0}, {0.0, 5.}
   grid:findCell(np, fIdx)
   assert_equal(1, fIdx[1], "Checking fIdx[1] for upper left corner point")
   assert_equal(20, fIdx[2], "Checking fIdx[2] for upper left corner point")
   local fIdx, np = {0, 0}, {2.0, 5.}
   grid:findCell(np, fIdx)
   assert_equal(10, fIdx[1], "Checking fIdx[1] for upper right corner point")
   assert_equal(20, fIdx[2], "Checking fIdx[2] for upper right corner point")
end

-- Run tests.
test_1()
test_2()
test_3()
test_4()
test_5a()
test_5()
test_6()
test_7()
test_8()
test_9()
test_10()
test_12()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
