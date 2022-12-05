local _M = {} 
_M[1] = function (z, b) 
   local z1, z2, z3, z4 = z[1], z[2], z[3], z[4] 
   b[1] = 0.25 
   b[2] = 0.4330127018922193*z1 
   b[3] = 0.4330127018922193*z2 
   b[4] = 0.4330127018922193*z3 
   b[5] = 0.4330127018922193*z4 
end 
_M[2] = function (z, b) 
   local z1, z2, z3, z4 = z[1], z[2], z[3], z[4] 
   b[1] = 0.25 
   b[2] = 0.4330127018922193*z1 
   b[3] = 0.4330127018922193*z2 
   b[4] = 0.4330127018922193*z3 
   b[5] = 0.4330127018922193*z4 
   b[6] = 0.75*z1*z2 
   b[7] = 0.75*z1*z3 
   b[8] = 0.75*z2*z3 
   b[9] = 0.75*z1*z4 
   b[10] = 0.75*z2*z4 
   b[11] = 0.75*z3*z4 
   b[12] = 0.8385254915624212*z1^2-0.2795084971874737 
   b[13] = 0.8385254915624212*z2^2-0.2795084971874737 
   b[14] = 0.8385254915624212*z3^2-0.2795084971874737 
   b[15] = 0.8385254915624212*z4^2-0.2795084971874737 
end 
_M[3] = function (z, b) 
   local z1, z2, z3, z4 = z[1], z[2], z[3], z[4] 
   b[1] = 0.25 
   b[2] = 0.4330127018922193*z1 
   b[3] = 0.4330127018922193*z2 
   b[4] = 0.4330127018922193*z3 
   b[5] = 0.4330127018922193*z4 
   b[6] = 0.75*z1*z2 
   b[7] = 0.75*z1*z3 
   b[8] = 0.75*z2*z3 
   b[9] = 0.75*z1*z4 
   b[10] = 0.75*z2*z4 
   b[11] = 0.75*z3*z4 
   b[12] = 0.8385254915624212*z1^2-0.2795084971874737 
   b[13] = 0.8385254915624212*z2^2-0.2795084971874737 
   b[14] = 0.8385254915624212*z3^2-0.2795084971874737 
   b[15] = 0.8385254915624212*z4^2-0.2795084971874737 
   b[16] = 1.299038105676658*z1*z2*z3 
   b[17] = 1.299038105676658*z1*z2*z4 
   b[18] = 1.299038105676658*z1*z3*z4 
   b[19] = 1.299038105676658*z2*z3*z4 
   b[20] = 1.452368754827781*z1^2*z2-0.4841229182759271*z2 
   b[21] = 1.452368754827781*z1*z2^2-0.4841229182759271*z1 
   b[22] = 1.452368754827781*z1^2*z3-0.4841229182759271*z3 
   b[23] = 1.452368754827781*z2^2*z3-0.4841229182759271*z3 
   b[24] = 1.452368754827781*z1*z3^2-0.4841229182759271*z1 
   b[25] = 1.452368754827781*z2*z3^2-0.4841229182759271*z2 
   b[26] = 1.452368754827781*z1^2*z4-0.4841229182759271*z4 
   b[27] = 1.452368754827781*z2^2*z4-0.4841229182759271*z4 
   b[28] = 1.452368754827781*z3^2*z4-0.4841229182759271*z4 
   b[29] = 1.452368754827781*z1*z4^2-0.4841229182759271*z1 
   b[30] = 1.452368754827781*z2*z4^2-0.4841229182759271*z2 
   b[31] = 1.452368754827781*z3*z4^2-0.4841229182759271*z3 
   b[32] = 1.653594569415369*z1^3-0.9921567416492214*z1 
   b[33] = 1.653594569415369*z2^3-0.9921567416492214*z2 
   b[34] = 1.653594569415369*z3^3-0.9921567416492214*z3 
   b[35] = 1.653594569415369*z4^3-0.9921567416492214*z4 
end
return _M 