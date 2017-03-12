local _M = {} 
_M[1] = function (z, b) 
   local z1, z2, z3 = z[1], z[2], z[3] 
   b[1] = 0.3535533905932734 
   b[2] = 0.6123724356957931*z1 
   b[3] = 0.6123724356957931*z2 
   b[4] = 0.6123724356957931*z3 
end 
_M[2] = function (z, b) 
   local z1, z2, z3 = z[1], z[2], z[3] 
   b[1] = 0.3535533905932734 
   b[2] = 0.6123724356957931*z1 
   b[3] = 0.6123724356957931*z2 
   b[4] = 0.6123724356957931*z3 
   b[5] = 1.060660171779822*z1*z2 
   b[6] = 1.060660171779822*z1*z3 
   b[7] = 1.060660171779822*z2*z3 
   b[8] = 1.185854122563141*z1^2-0.3952847075210471 
   b[9] = 1.185854122563141*z2^2-0.3952847075210471 
   b[10] = 1.185854122563141*z3^2-0.3952847075210471 
end 
_M[3] = function (z, b) 
   local z1, z2, z3 = z[1], z[2], z[3] 
   b[1] = 0.3535533905932734 
   b[2] = 0.6123724356957931*z1 
   b[3] = 0.6123724356957931*z2 
   b[4] = 0.6123724356957931*z3 
   b[5] = 1.060660171779822*z1*z2 
   b[6] = 1.060660171779822*z1*z3 
   b[7] = 1.060660171779822*z2*z3 
   b[8] = 1.185854122563141*z1^2-0.3952847075210471 
   b[9] = 1.185854122563141*z2^2-0.3952847075210471 
   b[10] = 1.185854122563141*z3^2-0.3952847075210471 
   b[11] = 1.837117307087383*z1*z2*z3 
   b[12] = 2.053959590644373*z1^2*z2-0.6846531968814578*z2 
   b[13] = 2.053959590644373*z1*z2^2-0.6846531968814578*z1 
   b[14] = 2.053959590644373*z1^2*z3-0.6846531968814578*z3 
   b[15] = 2.053959590644373*z2^2*z3-0.6846531968814578*z3 
   b[16] = 2.053959590644373*z1*z3^2-0.6846531968814578*z1 
   b[17] = 2.053959590644373*z2*z3^2-0.6846531968814578*z2 
   b[18] = 2.338535866733713*z1^3-1.403121520040228*z1 
   b[19] = 2.338535866733713*z2^3-1.403121520040228*z2 
   b[20] = 2.338535866733713*z3^3-1.403121520040228*z3 
end 
_M[4] = function (z, b) 
   local z1, z2, z3 = z[1], z[2], z[3] 
   b[1] = 0.3535533905932734 
   b[2] = 0.6123724356957931*z1 
   b[3] = 0.6123724356957931*z2 
   b[4] = 0.6123724356957931*z3 
   b[5] = 1.060660171779822*z1*z2 
   b[6] = 1.060660171779822*z1*z3 
   b[7] = 1.060660171779822*z2*z3 
   b[8] = 1.185854122563141*z1^2-0.3952847075210471 
   b[9] = 1.185854122563141*z2^2-0.3952847075210471 
   b[10] = 1.185854122563141*z3^2-0.3952847075210471 
   b[11] = 1.837117307087383*z1*z2*z3 
   b[12] = 2.053959590644373*z1^2*z2-0.6846531968814578*z2 
   b[13] = 2.053959590644373*z1*z2^2-0.6846531968814578*z1 
   b[14] = 2.053959590644373*z1^2*z3-0.6846531968814578*z3 
   b[15] = 2.053959590644373*z2^2*z3-0.6846531968814578*z3 
   b[16] = 2.053959590644373*z1*z3^2-0.6846531968814578*z1 
   b[17] = 2.053959590644373*z2*z3^2-0.6846531968814578*z2 
   b[18] = 2.338535866733713*z1^3-1.403121520040228*z1 
   b[19] = 2.338535866733713*z2^3-1.403121520040228*z2 
   b[20] = 2.338535866733713*z3^3-1.403121520040228*z3 
   b[21] = 3.557562367689424*z1^2*z2*z3-1.185854122563141*z2*z3 
   b[22] = 3.557562367689424*z1*z2^2*z3-1.185854122563141*z1*z3 
   b[23] = 3.557562367689424*z1*z2*z3^2-1.185854122563141*z1*z2 
   b[24] = 3.977475644174331*z1^2*z2^2-1.325825214724777*z2^2-1.325825214724777*z1^2+0.4419417382415923 
   b[25] = 3.977475644174331*z1^2*z3^2-1.325825214724777*z3^2-1.325825214724777*z1^2+0.4419417382415923 
   b[26] = 3.977475644174331*z2^2*z3^2-1.325825214724777*z3^2-1.325825214724777*z2^2+0.4419417382415923 
   b[27] = 4.050462936504911*z1^3*z2-2.430277761902947*z1*z2 
   b[28] = 4.050462936504911*z1*z2^3-2.430277761902947*z1*z2 
   b[29] = 4.050462936504911*z1^3*z3-2.430277761902947*z1*z3 
   b[30] = 4.050462936504911*z2^3*z3-2.430277761902947*z2*z3 
   b[31] = 4.050462936504911*z1*z3^3-2.430277761902947*z1*z3 
   b[32] = 4.050462936504911*z2*z3^3-2.430277761902947*z2*z3 
   b[33] = 4.640388251536713*z1^4-3.977475644174326*z1^2+0.3977475644174325 
   b[34] = 4.640388251536713*z2^4-3.977475644174326*z2^2+0.3977475644174325 
   b[35] = 4.640388251536713*z3^4-3.977475644174326*z3^2+0.3977475644174325 
end 
return _M 