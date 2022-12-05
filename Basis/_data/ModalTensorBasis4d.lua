local _M = {} 
_M[1] = function (z, b) 
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
   b[12] = 1.299038105676658*z1*z2*z3 
   b[13] = 1.299038105676658*z1*z2*z4 
   b[14] = 1.299038105676658*z1*z3*z4 
   b[15] = 1.299038105676658*z2*z3*z4 
   b[16] = 2.25*z1*z2*z3*z4 
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
   b[32] = 2.25*z1*z2*z3*z4 
   b[33] = 2.515576474687264*z1^2*z2*z3-0.8385254915624212*z2*z3 
   b[34] = 2.515576474687264*z1*z2^2*z3-0.8385254915624212*z1*z3 
   b[35] = 2.515576474687264*z1*z2*z3^2-0.8385254915624212*z1*z2 
   b[36] = 2.515576474687264*z1^2*z2*z4-0.8385254915624212*z2*z4 
   b[37] = 2.515576474687264*z1*z2^2*z4-0.8385254915624212*z1*z4 
   b[38] = 2.515576474687264*z1^2*z3*z4-0.8385254915624212*z3*z4 
   b[39] = 2.515576474687264*z2^2*z3*z4-0.8385254915624212*z3*z4 
   b[40] = 2.515576474687264*z1*z3^2*z4-0.8385254915624212*z1*z4 
   b[41] = 2.515576474687264*z2*z3^2*z4-0.8385254915624212*z2*z4 
   b[42] = 2.515576474687264*z1*z2*z4^2-0.8385254915624212*z1*z2 
   b[43] = 2.515576474687264*z1*z3*z4^2-0.8385254915624212*z1*z3 
   b[44] = 2.515576474687264*z2*z3*z4^2-0.8385254915624212*z2*z3 
   b[45] = 2.8125*z1^2*z2^2-0.9375*z2^2-0.9375*z1^2+0.3125 
   b[46] = 2.8125*z1^2*z3^2-0.9375*z3^2-0.9375*z1^2+0.3125 
   b[47] = 2.8125*z2^2*z3^2-0.9375*z3^2-0.9375*z2^2+0.3125 
   b[48] = 2.8125*z1^2*z4^2-0.9375*z4^2-0.9375*z1^2+0.3125 
   b[49] = 2.8125*z2^2*z4^2-0.9375*z4^2-0.9375*z2^2+0.3125 
   b[50] = 2.8125*z3^2*z4^2-0.9375*z4^2-0.9375*z3^2+0.3125 
   b[51] = 4.357106264483344*z1^2*z2*z3*z4-1.452368754827781*z2*z3*z4 
   b[52] = 4.357106264483344*z1*z2^2*z3*z4-1.452368754827781*z1*z3*z4 
   b[53] = 4.357106264483344*z1*z2*z3^2*z4-1.452368754827781*z1*z2*z4 
   b[54] = 4.357106264483344*z1*z2*z3*z4^2-1.452368754827781*z1*z2*z3 
   b[55] = 4.871392896287466*z1^2*z2^2*z3-1.623797632095822*z2^2*z3-1.623797632095822*z1^2*z3+0.541265877365274*z3 
   b[56] = 4.871392896287466*z1^2*z2*z3^2-1.623797632095822*z2*z3^2-1.623797632095822*z1^2*z2+0.541265877365274*z2 
   b[57] = 4.871392896287466*z1*z2^2*z3^2-1.623797632095822*z1*z3^2-1.623797632095822*z1*z2^2+0.541265877365274*z1 
   b[58] = 4.871392896287466*z1^2*z2^2*z4-1.623797632095822*z2^2*z4-1.623797632095822*z1^2*z4+0.541265877365274*z4 
   b[59] = 4.871392896287466*z1^2*z3^2*z4-1.623797632095822*z3^2*z4-1.623797632095822*z1^2*z4+0.541265877365274*z4 
   b[60] = 4.871392896287466*z2^2*z3^2*z4-1.623797632095822*z3^2*z4-1.623797632095822*z2^2*z4+0.541265877365274*z4 
   b[61] = 4.871392896287466*z1^2*z2*z4^2-1.623797632095822*z2*z4^2-1.623797632095822*z1^2*z2+0.541265877365274*z2 
   b[62] = 4.871392896287466*z1*z2^2*z4^2-1.623797632095822*z1*z4^2-1.623797632095822*z1*z2^2+0.541265877365274*z1 
   b[63] = 4.871392896287466*z1^2*z3*z4^2-1.623797632095822*z3*z4^2-1.623797632095822*z1^2*z3+0.541265877365274*z3 
   b[64] = 4.871392896287466*z2^2*z3*z4^2-1.623797632095822*z3*z4^2-1.623797632095822*z2^2*z3+0.541265877365274*z3 
   b[65] = 4.871392896287466*z1*z3^2*z4^2-1.623797632095822*z1*z4^2-1.623797632095822*z1*z3^2+0.541265877365274*z1 
   b[66] = 4.871392896287466*z2*z3^2*z4^2-1.623797632095822*z2*z4^2-1.623797632095822*z2*z3^2+0.541265877365274*z2 
   b[67] = 8.4375*z1^2*z2^2*z3*z4-2.8125*z2^2*z3*z4-2.8125*z1^2*z3*z4+0.9375*z3*z4 
   b[68] = 8.4375*z1^2*z2*z3^2*z4-2.8125*z2*z3^2*z4-2.8125*z1^2*z2*z4+0.9375*z2*z4 
   b[69] = 8.4375*z1*z2^2*z3^2*z4-2.8125*z1*z3^2*z4-2.8125*z1*z2^2*z4+0.9375*z1*z4 
   b[70] = 8.4375*z1^2*z2*z3*z4^2-2.8125*z2*z3*z4^2-2.8125*z1^2*z2*z3+0.9375*z2*z3 
   b[71] = 8.4375*z1*z2^2*z3*z4^2-2.8125*z1*z3*z4^2-2.8125*z1*z2^2*z3+0.9375*z1*z3 
   b[72] = 8.4375*z1*z2*z3^2*z4^2-2.8125*z1*z2*z4^2-2.8125*z1*z2*z3^2+0.9375*z1*z2 
   b[73] = 9.43341178007724*z1^2*z2^2*z3^2-3.14447059335908*z2^2*z3^2-3.14447059335908*z1^2*z3^2+1.048156864453027*z3^2-3.14447059335908*z1^2*z2^2+1.048156864453027*z2^2+1.048156864453027*z1^2-0.3493856214843422 
   b[74] = 9.43341178007724*z1^2*z2^2*z4^2-3.14447059335908*z2^2*z4^2-3.14447059335908*z1^2*z4^2+1.048156864453027*z4^2-3.14447059335908*z1^2*z2^2+1.048156864453027*z2^2+1.048156864453027*z1^2-0.3493856214843422 
   b[75] = 9.43341178007724*z1^2*z3^2*z4^2-3.14447059335908*z3^2*z4^2-3.14447059335908*z1^2*z4^2+1.048156864453027*z4^2-3.14447059335908*z1^2*z3^2+1.048156864453027*z3^2+1.048156864453027*z1^2-0.3493856214843422 
   b[76] = 9.43341178007724*z2^2*z3^2*z4^2-3.14447059335908*z3^2*z4^2-3.14447059335908*z2^2*z4^2+1.048156864453027*z4^2-3.14447059335908*z2^2*z3^2+1.048156864453027*z3^2+1.048156864453027*z2^2-0.3493856214843422 
   b[77] = 16.33914849181254*z1^2*z2^2*z3^2*z4-5.44638283060418*z2^2*z3^2*z4-5.44638283060418*z1^2*z3^2*z4+1.815460943534727*z3^2*z4-5.44638283060418*z1^2*z2^2*z4+1.815460943534727*z2^2*z4+1.815460943534727*z1^2*z4-0.6051536478449089*z4 
   b[78] = 16.33914849181254*z1^2*z2^2*z3*z4^2-5.44638283060418*z2^2*z3*z4^2-5.44638283060418*z1^2*z3*z4^2+1.815460943534727*z3*z4^2-5.44638283060418*z1^2*z2^2*z3+1.815460943534727*z2^2*z3+1.815460943534727*z1^2*z3-0.6051536478449089*z3 
   b[79] = 16.33914849181254*z1^2*z2*z3^2*z4^2-5.44638283060418*z2*z3^2*z4^2-5.44638283060418*z1^2*z2*z4^2+1.815460943534727*z2*z4^2-5.44638283060418*z1^2*z2*z3^2+1.815460943534727*z2*z3^2+1.815460943534727*z1^2*z2-0.6051536478449089*z2 
   b[80] = 16.33914849181254*z1*z2^2*z3^2*z4^2-5.44638283060418*z1*z3^2*z4^2-5.44638283060418*z1*z2^2*z4^2+1.815460943534727*z1*z4^2-5.44638283060418*z1*z2^2*z3^2+1.815460943534727*z1*z3^2+1.815460943534727*z1*z2^2-0.6051536478449089*z1 
   b[81] = 31.640625*z1^2*z2^2*z3^2*z4^2-10.546875*z2^2*z3^2*z4^2-10.546875*z1^2*z3^2*z4^2+3.515625*z3^2*z4^2-10.546875*z1^2*z2^2*z4^2+3.515625*z2^2*z4^2+3.515625*z1^2*z4^2-1.171875*z4^2-10.546875*z1^2*z2^2*z3^2+3.515625*z2^2*z3^2+3.515625*z1^2*z3^2-1.171875*z3^2+3.515625*z1^2*z2^2-1.171875*z2^2-1.171875*z1^2+0.390625 
end 
return _M 