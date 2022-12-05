local _M = {} 
_M[1] = function (z, b) 
   local z1, z2, z3, z4, z5 = z[1], z[2], z[3], z[4], z[5] 
   b[1] = 0.1767766952966368 
   b[2] = 0.3061862178478971*z1 
   b[3] = 0.3061862178478971*z2 
   b[4] = 0.3061862178478971*z3 
   b[5] = 0.3061862178478971*z4 
   b[6] = 0.3061862178478971*z5 
end 
_M[2] = function (z, b) 
   local z1, z2, z3, z4, z5 = z[1], z[2], z[3], z[4], z[5] 
   b[1] = 0.1767766952966368 
   b[2] = 0.3061862178478971*z1 
   b[3] = 0.3061862178478971*z2 
   b[4] = 0.3061862178478971*z3 
   b[5] = 0.3061862178478971*z4 
   b[6] = 0.3061862178478971*z5 
   b[7] = 0.5303300858899105*z1*z2 
   b[8] = 0.5303300858899105*z1*z3 
   b[9] = 0.5303300858899105*z2*z3 
   b[10] = 0.5303300858899105*z1*z4 
   b[11] = 0.5303300858899105*z2*z4 
   b[12] = 0.5303300858899105*z3*z4 
   b[13] = 0.5303300858899105*z1*z5 
   b[14] = 0.5303300858899105*z2*z5 
   b[15] = 0.5303300858899105*z3*z5 
   b[16] = 0.5303300858899105*z4*z5 
   b[17] = 0.592927061281571*z1^2-0.1976423537605237 
   b[18] = 0.592927061281571*z2^2-0.1976423537605237 
   b[19] = 0.592927061281571*z3^2-0.1976423537605237 
   b[20] = 0.592927061281571*z4^2-0.1976423537605237 
   b[21] = 0.592927061281571*z5^2-0.1976423537605237 
end 
_M[3] = function (z, b) 
   local z1, z2, z3, z4, z5 = z[1], z[2], z[3], z[4], z[5] 
   b[1] = 0.1767766952966368 
   b[2] = 0.3061862178478971*z1 
   b[3] = 0.3061862178478971*z2 
   b[4] = 0.3061862178478971*z3 
   b[5] = 0.3061862178478971*z4 
   b[6] = 0.3061862178478971*z5 
   b[7] = 0.5303300858899105*z1*z2 
   b[8] = 0.5303300858899105*z1*z3 
   b[9] = 0.5303300858899105*z2*z3 
   b[10] = 0.5303300858899105*z1*z4 
   b[11] = 0.5303300858899105*z2*z4 
   b[12] = 0.5303300858899105*z3*z4 
   b[13] = 0.5303300858899105*z1*z5 
   b[14] = 0.5303300858899105*z2*z5 
   b[15] = 0.5303300858899105*z3*z5 
   b[16] = 0.5303300858899105*z4*z5 
   b[17] = 0.592927061281571*z1^2-0.1976423537605237 
   b[18] = 0.592927061281571*z2^2-0.1976423537605237 
   b[19] = 0.592927061281571*z3^2-0.1976423537605237 
   b[20] = 0.592927061281571*z4^2-0.1976423537605237 
   b[21] = 0.592927061281571*z5^2-0.1976423537605237 
   b[22] = 0.9185586535436913*z1*z2*z3 
   b[23] = 0.9185586535436913*z1*z2*z4 
   b[24] = 0.9185586535436913*z1*z3*z4 
   b[25] = 0.9185586535436913*z2*z3*z4 
   b[26] = 0.9185586535436913*z1*z2*z5 
   b[27] = 0.9185586535436913*z1*z3*z5 
   b[28] = 0.9185586535436913*z2*z3*z5 
   b[29] = 0.9185586535436913*z1*z4*z5 
   b[30] = 0.9185586535436913*z2*z4*z5 
   b[31] = 0.9185586535436913*z3*z4*z5 
   b[32] = 1.026979795322186*z1^2*z2-0.3423265984407287*z2 
   b[33] = 1.026979795322186*z1*z2^2-0.3423265984407287*z1 
   b[34] = 1.026979795322186*z1^2*z3-0.3423265984407287*z3 
   b[35] = 1.026979795322186*z2^2*z3-0.3423265984407287*z3 
   b[36] = 1.026979795322186*z1*z3^2-0.3423265984407287*z1 
   b[37] = 1.026979795322186*z2*z3^2-0.3423265984407287*z2 
   b[38] = 1.026979795322186*z1^2*z4-0.3423265984407287*z4 
   b[39] = 1.026979795322186*z2^2*z4-0.3423265984407287*z4 
   b[40] = 1.026979795322186*z3^2*z4-0.3423265984407287*z4 
   b[41] = 1.026979795322186*z1*z4^2-0.3423265984407287*z1 
   b[42] = 1.026979795322186*z2*z4^2-0.3423265984407287*z2 
   b[43] = 1.026979795322186*z3*z4^2-0.3423265984407287*z3 
   b[44] = 1.026979795322186*z1^2*z5-0.3423265984407287*z5 
   b[45] = 1.026979795322186*z2^2*z5-0.3423265984407287*z5 
   b[46] = 1.026979795322186*z3^2*z5-0.3423265984407287*z5 
   b[47] = 1.026979795322186*z4^2*z5-0.3423265984407287*z5 
   b[48] = 1.026979795322186*z1*z5^2-0.3423265984407287*z1 
   b[49] = 1.026979795322186*z2*z5^2-0.3423265984407287*z2 
   b[50] = 1.026979795322186*z3*z5^2-0.3423265984407287*z3 
   b[51] = 1.026979795322186*z4*z5^2-0.3423265984407287*z4 
   b[52] = 1.169267933366856*z1^3-0.7015607600201138*z1 
   b[53] = 1.169267933366856*z2^3-0.7015607600201138*z2 
   b[54] = 1.169267933366856*z3^3-0.7015607600201138*z3 
   b[55] = 1.169267933366856*z4^3-0.7015607600201138*z4 
   b[56] = 1.169267933366856*z5^3-0.7015607600201138*z5 
end
return _M 