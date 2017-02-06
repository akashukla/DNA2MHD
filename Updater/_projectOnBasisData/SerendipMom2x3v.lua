local _M = { numDensity = {}, momentum = {} } 
 
_M.numDensity[1] = function (f, out, dv, w) 
   out[1] = out[1] + 2.828427124746191*f[1] 
   out[2] = out[2] + 2.828427124746191*f[2] 
   out[3] = out[3] + 2.828427124746191*f[3] 
   out[4] = out[4] + 2.828427124746191*f[7] 
end 
_M.numDensity[2] = function (f, out, dv, w) 
   out[1] = out[1] + 2.828427124746191*f[1] 
   out[2] = out[2] + 2.828427124746191*f[2] 
   out[3] = out[3] + 2.828427124746191*f[3] 
   out[4] = out[4] + 2.828427124746191*f[7] 
   out[5] = out[5] + 2.828427124746191*f[17] 
   out[6] = out[6] + 2.828427124746191*f[18] 
   out[7] = out[7] + 2.828427124746191*f[32] 
   out[8] = out[8] + 2.828427124746191*f[33] 
end 
_M.numDensity[3] = function (f, out, dv, w) 
   out[1] = out[1] + 2.828427124746191*f[1] 
   out[2] = out[2] + 2.828427124746191*f[2] 
   out[3] = out[3] + 2.828427124746191*f[3] 
   out[4] = out[4] + 2.828427124746191*f[7] 
   out[5] = out[5] + 2.828427124746191*f[17] 
   out[6] = out[6] + 2.828427124746191*f[18] 
   out[7] = out[7] + 2.828427124746191*f[32] 
   out[8] = out[8] + 2.828427124746191*f[33] 
   out[9] = out[9] + 2.828427124746191*f[52] 
   out[10] = out[10] + 2.828427124746191*f[53] 
   out[11] = out[11] + 2.828427124746191*f[92] 
   out[12] = out[12] + 2.828427124746191*f[93] 
end 
_M.numDensity[4] = function (f, out, dv, w) 
   out[1] = out[1] + 2.828427124746191*f[1] 
   out[2] = out[2] + 2.828427124746191*f[2] 
   out[3] = out[3] + 2.828427124746191*f[3] 
   out[4] = out[4] + 2.828427124746191*f[7] 
   out[5] = out[5] + 2.828427124746191*f[17] 
   out[6] = out[6] + 2.828427124746191*f[18] 
   out[7] = out[7] + 2.828427124746191*f[32] 
   out[8] = out[8] + 2.828427124746191*f[33] 
   out[9] = out[9] + 2.828427124746191*f[52] 
   out[10] = out[10] + 2.828427124746191*f[53] 
   out[11] = out[11] + 2.828427124746191*f[92] 
   out[12] = out[12] + 2.828427124746191*f[102] 
   out[13] = out[13] + 2.828427124746191*f[103] 
   out[14] = out[14] + 2.828427124746191*f[122] 
   out[15] = out[15] + 2.828427124746191*f[123] 
   out[16] = out[16] + 2.828427124746191*f[208] 
   out[17] = out[17] + 2.828427124746191*f[209] 
end 

 
_M.momentum[1] = function (f, out, dv, w) 
   local dv1, dv2, dv3 = dv[1], dv[2], dv[3] 
   local w1, w2, w3 = w[1], w[2], w[3] 
   out[1] = out[1] + 2.828427124746191*f[1]*w1+0.8164965809277261*f[4]*dv1 
   out[2] = out[2] + 2.828427124746191*f[2]*w1+0.8164965809277261*f[8]*dv1 
   out[3] = out[3] + 2.828427124746191*f[3]*w1+0.8164965809277261*f[9]*dv1 
   out[4] = out[4] + 2.828427124746191*f[7]*w1+0.8164965809277261*f[17]*dv1 
   out[5] = out[5] + 2.828427124746191*f[1]*w2+0.8164965809277261*f[5]*dv2 
   out[6] = out[6] + 2.828427124746191*f[2]*w2+0.8164965809277261*f[10]*dv2 
   out[7] = out[7] + 2.828427124746191*f[3]*w2+0.8164965809277261*f[11]*dv2 
   out[8] = out[8] + 2.828427124746191*f[7]*w2+0.8164965809277261*f[18]*dv2 
   out[9] = out[9] + 2.828427124746191*f[1]*w3+0.8164965809277261*f[6]*dv3 
   out[10] = out[10] + 2.828427124746191*f[2]*w3+0.8164965809277261*f[13]*dv3 
   out[11] = out[11] + 2.828427124746191*f[3]*w3+0.8164965809277261*f[14]*dv3 
   out[12] = out[12] + 2.828427124746191*f[7]*w3+0.8164965809277261*f[21]*dv3 
end 
_M.momentum[2] = function (f, out, dv, w) 
   local dv1, dv2, dv3 = dv[1], dv[2], dv[3] 
   local w1, w2, w3 = w[1], w[2], w[3] 
   out[1] = out[1] + 2.828427124746191*f[1]*w1+0.8164965809277261*f[4]*dv1 
   out[2] = out[2] + 2.828427124746191*f[2]*w1+0.8164965809277261*f[8]*dv1 
   out[3] = out[3] + 2.828427124746191*f[3]*w1+0.8164965809277261*f[9]*dv1 
   out[4] = out[4] + 2.828427124746191*f[7]*w1+0.8164965809277261*f[22]*dv1 
   out[5] = out[5] + 2.828427124746191*f[17]*w1+0.8164965809277261*f[34]*dv1 
   out[6] = out[6] + 2.828427124746191*f[18]*w1+0.8164965809277261*f[35]*dv1 
   out[7] = out[7] + 2.828427124746191*f[32]*w1+0.8164965809277261*f[57]*dv1 
   out[8] = out[8] + 2.828427124746191*f[33]*w1+0.8164965809277261*f[58]*dv1 
   out[9] = out[9] + 2.828427124746191*f[1]*w2+0.8164965809277261*f[5]*dv2 
   out[10] = out[10] + 2.828427124746191*f[2]*w2+0.8164965809277261*f[10]*dv2 
   out[11] = out[11] + 2.828427124746191*f[3]*w2+0.8164965809277261*f[11]*dv2 
   out[12] = out[12] + 2.828427124746191*f[7]*w2+0.8164965809277261*f[23]*dv2 
   out[13] = out[13] + 2.828427124746191*f[17]*w2+0.8164965809277261*f[38]*dv2 
   out[14] = out[14] + 2.828427124746191*f[18]*w2+0.8164965809277261*f[39]*dv2 
   out[15] = out[15] + 2.828427124746191*f[32]*w2+0.8164965809277261*f[60]*dv2 
   out[16] = out[16] + 2.828427124746191*f[33]*w2+0.8164965809277261*f[61]*dv2 
   out[17] = out[17] + 2.828427124746191*f[1]*w3+0.8164965809277261*f[6]*dv3 
   out[18] = out[18] + 2.828427124746191*f[2]*w3+0.8164965809277261*f[13]*dv3 
   out[19] = out[19] + 2.828427124746191*f[3]*w3+0.8164965809277261*f[14]*dv3 
   out[20] = out[20] + 2.828427124746191*f[7]*w3+0.8164965809277261*f[26]*dv3 
   out[21] = out[21] + 2.828427124746191*f[17]*w3+0.8164965809277261*f[44]*dv3 
   out[22] = out[22] + 2.828427124746191*f[18]*w3+0.8164965809277261*f[45]*dv3 
   out[23] = out[23] + 2.828427124746191*f[32]*w3+0.8164965809277261*f[69]*dv3 
   out[24] = out[24] + 2.828427124746191*f[33]*w3+0.8164965809277261*f[70]*dv3 
end 
_M.momentum[3] = function (f, out, dv, w) 
   local dv1, dv2, dv3 = dv[1], dv[2], dv[3] 
   local w1, w2, w3 = w[1], w[2], w[3] 
   out[1] = out[1] + 2.828427124746191*f[1]*w1+0.8164965809277261*f[4]*dv1 
   out[2] = out[2] + 2.828427124746191*f[2]*w1+0.8164965809277261*f[8]*dv1 
   out[3] = out[3] + 2.828427124746191*f[3]*w1+0.8164965809277261*f[9]*dv1 
   out[4] = out[4] + 2.828427124746191*f[7]*w1+0.8164965809277261*f[22]*dv1 
   out[5] = out[5] + 2.828427124746191*f[17]*w1+0.8164965809277261*f[34]*dv1 
   out[6] = out[6] + 2.828427124746191*f[18]*w1+0.8164965809277261*f[35]*dv1 
   out[7] = out[7] + 2.828427124746191*f[32]*w1+0.8164965809277261*f[62]*dv1 
   out[8] = out[8] + 2.828427124746191*f[33]*w1+0.8164965809277261*f[63]*dv1 
   out[9] = out[9] + 2.828427124746191*f[52]*w1+0.8164965809277261*f[94]*dv1 
   out[10] = out[10] + 2.828427124746191*f[53]*w1+0.8164965809277261*f[95]*dv1 
   out[11] = out[11] + 2.828427124746191*f[92]*w1+0.8164965809277261*f[133]*dv1 
   out[12] = out[12] + 2.828427124746191*f[93]*w1+0.8164965809277261*f[134]*dv1 
   out[13] = out[13] + 2.828427124746191*f[1]*w2+0.8164965809277261*f[5]*dv2 
   out[14] = out[14] + 2.828427124746191*f[2]*w2+0.8164965809277261*f[10]*dv2 
   out[15] = out[15] + 2.828427124746191*f[3]*w2+0.8164965809277261*f[11]*dv2 
   out[16] = out[16] + 2.828427124746191*f[7]*w2+0.8164965809277261*f[23]*dv2 
   out[17] = out[17] + 2.828427124746191*f[17]*w2+0.8164965809277261*f[38]*dv2 
   out[18] = out[18] + 2.828427124746191*f[18]*w2+0.8164965809277261*f[39]*dv2 
   out[19] = out[19] + 2.828427124746191*f[32]*w2+0.8164965809277261*f[65]*dv2 
   out[20] = out[20] + 2.828427124746191*f[33]*w2+0.8164965809277261*f[66]*dv2 
   out[21] = out[21] + 2.828427124746191*f[52]*w2+0.8164965809277261*f[98]*dv2 
   out[22] = out[22] + 2.828427124746191*f[53]*w2+0.8164965809277261*f[99]*dv2 
   out[23] = out[23] + 2.828427124746191*f[92]*w2+0.8164965809277261*f[136]*dv2 
   out[24] = out[24] + 2.828427124746191*f[93]*w2+0.8164965809277261*f[137]*dv2 
   out[25] = out[25] + 2.828427124746191*f[1]*w3+0.8164965809277261*f[6]*dv3 
   out[26] = out[26] + 2.828427124746191*f[2]*w3+0.8164965809277261*f[13]*dv3 
   out[27] = out[27] + 2.828427124746191*f[3]*w3+0.8164965809277261*f[14]*dv3 
   out[28] = out[28] + 2.828427124746191*f[7]*w3+0.8164965809277261*f[26]*dv3 
   out[29] = out[29] + 2.828427124746191*f[17]*w3+0.8164965809277261*f[44]*dv3 
   out[30] = out[30] + 2.828427124746191*f[18]*w3+0.8164965809277261*f[45]*dv3 
   out[31] = out[31] + 2.828427124746191*f[32]*w3+0.8164965809277261*f[74]*dv3 
   out[32] = out[32] + 2.828427124746191*f[33]*w3+0.8164965809277261*f[75]*dv3 
   out[33] = out[33] + 2.828427124746191*f[52]*w3+0.8164965809277261*f[104]*dv3 
   out[34] = out[34] + 2.828427124746191*f[53]*w3+0.8164965809277261*f[105]*dv3 
   out[35] = out[35] + 2.828427124746191*f[92]*w3+0.8164965809277261*f[145]*dv3 
   out[36] = out[36] + 2.828427124746191*f[93]*w3+0.8164965809277261*f[146]*dv3 
end 
_M.momentum[4] = function (f, out, dv, w) 
   local dv1, dv2, dv3 = dv[1], dv[2], dv[3] 
   local w1, w2, w3 = w[1], w[2], w[3] 
   out[1] = out[1] + 2.828427124746191*f[1]*w1+0.8164965809277261*f[4]*dv1 
   out[2] = out[2] + 2.828427124746191*f[2]*w1+0.8164965809277261*f[8]*dv1 
   out[3] = out[3] + 2.828427124746191*f[3]*w1+0.8164965809277261*f[9]*dv1 
   out[4] = out[4] + 2.828427124746191*f[7]*w1+0.8164965809277261*f[22]*dv1 
   out[5] = out[5] + 2.828427124746191*f[17]*w1+0.8164965809277261*f[34]*dv1 
   out[6] = out[6] + 2.828427124746191*f[18]*w1+0.8164965809277261*f[35]*dv1 
   out[7] = out[7] + 2.828427124746191*f[32]*w1+0.8164965809277261*f[62]*dv1 
   out[8] = out[8] + 2.828427124746191*f[33]*w1+0.8164965809277261*f[63]*dv1 
   out[9] = out[9] + 2.828427124746191*f[52]*w1+0.8164965809277261*f[104]*dv1 
   out[10] = out[10] + 2.828427124746191*f[53]*w1+0.8164965809277261*f[105]*dv1 
   out[11] = out[11] + 2.828427124746191*f[92]*w1+0.8164965809277261*f[148]*dv1 
   out[12] = out[12] + 2.828427124746191*f[102]*w1+0.8164965809277261*f[178]*dv1 
   out[13] = out[13] + 2.828427124746191*f[103]*w1+0.8164965809277261*f[179]*dv1 
   out[14] = out[14] + 2.828427124746191*f[122]*w1+0.8164965809277261*f[210]*dv1 
   out[15] = out[15] + 2.828427124746191*f[123]*w1+0.8164965809277261*f[211]*dv1 
   out[16] = out[16] + 2.828427124746191*f[208]*w1+0.8164965809277261*f[283]*dv1 
   out[17] = out[17] + 2.828427124746191*f[209]*w1+0.8164965809277261*f[284]*dv1 
   out[18] = out[18] + 2.828427124746191*f[1]*w2+0.8164965809277261*f[5]*dv2 
   out[19] = out[19] + 2.828427124746191*f[2]*w2+0.8164965809277261*f[10]*dv2 
   out[20] = out[20] + 2.828427124746191*f[3]*w2+0.8164965809277261*f[11]*dv2 
   out[21] = out[21] + 2.828427124746191*f[7]*w2+0.8164965809277261*f[23]*dv2 
   out[22] = out[22] + 2.828427124746191*f[17]*w2+0.8164965809277261*f[38]*dv2 
   out[23] = out[23] + 2.828427124746191*f[18]*w2+0.8164965809277261*f[39]*dv2 
   out[24] = out[24] + 2.828427124746191*f[32]*w2+0.8164965809277261*f[65]*dv2 
   out[25] = out[25] + 2.828427124746191*f[33]*w2+0.8164965809277261*f[66]*dv2 
   out[26] = out[26] + 2.828427124746191*f[52]*w2+0.8164965809277261*f[108]*dv2 
   out[27] = out[27] + 2.828427124746191*f[53]*w2+0.8164965809277261*f[109]*dv2 
   out[28] = out[28] + 2.828427124746191*f[92]*w2+0.8164965809277261*f[151]*dv2 
   out[29] = out[29] + 2.828427124746191*f[102]*w2+0.8164965809277261*f[181]*dv2 
   out[30] = out[30] + 2.828427124746191*f[103]*w2+0.8164965809277261*f[182]*dv2 
   out[31] = out[31] + 2.828427124746191*f[122]*w2+0.8164965809277261*f[214]*dv2 
   out[32] = out[32] + 2.828427124746191*f[123]*w2+0.8164965809277261*f[215]*dv2 
   out[33] = out[33] + 2.828427124746191*f[208]*w2+0.8164965809277261*f[286]*dv2 
   out[34] = out[34] + 2.828427124746191*f[209]*w2+0.8164965809277261*f[287]*dv2 
   out[35] = out[35] + 2.828427124746191*f[1]*w3+0.8164965809277261*f[6]*dv3 
   out[36] = out[36] + 2.828427124746191*f[2]*w3+0.8164965809277261*f[13]*dv3 
   out[37] = out[37] + 2.828427124746191*f[3]*w3+0.8164965809277261*f[14]*dv3 
   out[38] = out[38] + 2.828427124746191*f[7]*w3+0.8164965809277261*f[26]*dv3 
   out[39] = out[39] + 2.828427124746191*f[17]*w3+0.8164965809277261*f[44]*dv3 
   out[40] = out[40] + 2.828427124746191*f[18]*w3+0.8164965809277261*f[45]*dv3 
   out[41] = out[41] + 2.828427124746191*f[32]*w3+0.8164965809277261*f[74]*dv3 
   out[42] = out[42] + 2.828427124746191*f[33]*w3+0.8164965809277261*f[75]*dv3 
   out[43] = out[43] + 2.828427124746191*f[52]*w3+0.8164965809277261*f[114]*dv3 
   out[44] = out[44] + 2.828427124746191*f[53]*w3+0.8164965809277261*f[115]*dv3 
   out[45] = out[45] + 2.828427124746191*f[92]*w3+0.8164965809277261*f[160]*dv3 
   out[46] = out[46] + 2.828427124746191*f[102]*w3+0.8164965809277261*f[190]*dv3 
   out[47] = out[47] + 2.828427124746191*f[103]*w3+0.8164965809277261*f[191]*dv3 
   out[48] = out[48] + 2.828427124746191*f[122]*w3+0.8164965809277261*f[220]*dv3 
   out[49] = out[49] + 2.828427124746191*f[123]*w3+0.8164965809277261*f[221]*dv3 
   out[50] = out[50] + 2.828427124746191*f[208]*w3+0.8164965809277261*f[295]*dv3 
   out[51] = out[51] + 2.828427124746191*f[209]*w3+0.8164965809277261*f[296]*dv3 
end 
return _M 
