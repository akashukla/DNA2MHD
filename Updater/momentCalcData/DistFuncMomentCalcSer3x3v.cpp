#include <DistFuncMomentCalcModDecl.h> 
void MomentCalc3x3vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[3]*volFact; 
  out[4] += 2.828427124746191*f[7]*volFact; 
  out[5] += 2.828427124746191*f[8]*volFact; 
  out[6] += 2.828427124746191*f[9]*volFact; 
  out[7] += 2.828427124746191*f[22]*volFact; 
} 
void MomentCalc3x3vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  out[0] += (2.828427124746191*f[0]*wx1+0.8164965809277261*f[4]*dv1)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx1+0.8164965809277261*f[10]*dv1)*volFact; 
  out[2] += (2.828427124746191*f[2]*wx1+0.8164965809277261*f[11]*dv1)*volFact; 
  out[3] += (2.828427124746191*f[3]*wx1+0.8164965809277261*f[12]*dv1)*volFact; 
  out[4] += (2.828427124746191*f[7]*wx1+0.8164965809277261*f[23]*dv1)*volFact; 
  out[5] += (2.828427124746191*f[8]*wx1+0.8164965809277261*f[24]*dv1)*volFact; 
  out[6] += (2.828427124746191*f[9]*wx1+0.8164965809277261*f[25]*dv1)*volFact; 
  out[7] += (2.828427124746191*f[22]*wx1+0.8164965809277261*f[42]*dv1)*volFact; 
  out[8] += (2.828427124746191*f[0]*wx2+0.8164965809277261*f[5]*dv2)*volFact; 
  out[9] += (2.828427124746191*f[1]*wx2+0.8164965809277261*f[13]*dv2)*volFact; 
  out[10] += (2.828427124746191*f[2]*wx2+0.8164965809277261*f[14]*dv2)*volFact; 
  out[11] += (2.828427124746191*f[3]*wx2+0.8164965809277261*f[15]*dv2)*volFact; 
  out[12] += (2.828427124746191*f[7]*wx2+0.8164965809277261*f[26]*dv2)*volFact; 
  out[13] += (2.828427124746191*f[8]*wx2+0.8164965809277261*f[27]*dv2)*volFact; 
  out[14] += (2.828427124746191*f[9]*wx2+0.8164965809277261*f[28]*dv2)*volFact; 
  out[15] += (2.828427124746191*f[22]*wx2+0.8164965809277261*f[43]*dv2)*volFact; 
  out[16] += (2.828427124746191*f[0]*wx3+0.8164965809277261*f[6]*dv3)*volFact; 
  out[17] += (2.828427124746191*f[1]*wx3+0.8164965809277261*f[17]*dv3)*volFact; 
  out[18] += (2.828427124746191*f[2]*wx3+0.8164965809277261*f[18]*dv3)*volFact; 
  out[19] += (2.828427124746191*f[3]*wx3+0.8164965809277261*f[19]*dv3)*volFact; 
  out[20] += (2.828427124746191*f[7]*wx3+0.8164965809277261*f[32]*dv3)*volFact; 
  out[21] += (2.828427124746191*f[8]*wx3+0.8164965809277261*f[33]*dv3)*volFact; 
  out[22] += (2.828427124746191*f[9]*wx3+0.8164965809277261*f[34]*dv3)*volFact; 
  out[23] += (2.828427124746191*f[22]*wx3+0.8164965809277261*f[47]*dv3)*volFact; 
} 
void MomentCalc3x3vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  double zeta1[64]; 

  zeta1[0] = 8.0*wx1_sq+0.6666666666666666*dv1_sq; 
  zeta1[4] = 4.618802153517007*dv1*wx1; 
  double zeta2[64]; 

  zeta2[0] = 8.0*wx1*wx2; 
  zeta2[4] = 2.309401076758503*dv1*wx2; 
  zeta2[5] = 2.309401076758503*dv2*wx1; 
  zeta2[16] = 0.6666666666666666*dv1*dv2; 
  double zeta3[64]; 

  zeta3[0] = 8.0*wx1*wx3; 
  zeta3[4] = 2.309401076758503*dv1*wx3; 
  zeta3[6] = 2.309401076758503*dv3*wx1; 
  zeta3[20] = 0.6666666666666666*dv1*dv3; 
  double zeta4[64]; 

  zeta4[0] = 8.0*wx2_sq+0.6666666666666666*dv2_sq; 
  zeta4[5] = 4.618802153517007*dv2*wx2; 
  double zeta5[64]; 

  zeta5[0] = 8.0*wx2*wx3; 
  zeta5[5] = 2.309401076758503*dv2*wx3; 
  zeta5[6] = 2.309401076758503*dv3*wx2; 
  zeta5[21] = 0.6666666666666666*dv2*dv3; 
  double zeta6[64]; 

  zeta6[0] = 8.0*wx3_sq+0.6666666666666666*dv3_sq; 
  zeta6[6] = 4.618802153517007*dv3*wx3; 
  out[0] += (0.3535533905932737*f[4]*zeta1[4]+0.3535533905932737*f[0]*zeta1[0])*volFact; 
  out[1] += (0.3535533905932737*zeta1[4]*f[10]+0.3535533905932737*zeta1[0]*f[1])*volFact; 
  out[2] += (0.3535533905932737*zeta1[4]*f[11]+0.3535533905932737*zeta1[0]*f[2])*volFact; 
  out[3] += (0.3535533905932737*zeta1[4]*f[12]+0.3535533905932737*zeta1[0]*f[3])*volFact; 
  out[4] += (0.3535533905932737*zeta1[4]*f[23]+0.3535533905932737*zeta1[0]*f[7])*volFact; 
  out[5] += (0.3535533905932737*zeta1[4]*f[24]+0.3535533905932737*zeta1[0]*f[8])*volFact; 
  out[6] += (0.3535533905932737*zeta1[4]*f[25]+0.3535533905932737*zeta1[0]*f[9])*volFact; 
  out[7] += (0.3535533905932737*zeta1[4]*f[42]+0.3535533905932737*zeta1[0]*f[22])*volFact; 
  out[8] += (0.3535533905932737*f[16]*zeta2[16]+0.3535533905932737*f[5]*zeta2[5]+0.3535533905932737*f[4]*zeta2[4]+0.3535533905932737*f[0]*zeta2[0])*volFact; 
  out[9] += (0.3535533905932737*zeta2[16]*f[29]+0.3535533905932737*zeta2[5]*f[13]+0.3535533905932737*zeta2[4]*f[10]+0.3535533905932737*zeta2[0]*f[1])*volFact; 
  out[10] += (0.3535533905932737*zeta2[16]*f[30]+0.3535533905932737*zeta2[5]*f[14]+0.3535533905932737*zeta2[4]*f[11]+0.3535533905932737*zeta2[0]*f[2])*volFact; 
  out[11] += (0.3535533905932737*zeta2[16]*f[31]+0.3535533905932737*zeta2[5]*f[15]+0.3535533905932737*zeta2[4]*f[12]+0.3535533905932737*zeta2[0]*f[3])*volFact; 
  out[12] += (0.3535533905932737*zeta2[16]*f[44]+0.3535533905932737*zeta2[5]*f[26]+0.3535533905932737*zeta2[4]*f[23]+0.3535533905932737*zeta2[0]*f[7])*volFact; 
  out[13] += (0.3535533905932737*zeta2[16]*f[45]+0.3535533905932737*zeta2[5]*f[27]+0.3535533905932737*zeta2[4]*f[24]+0.3535533905932737*zeta2[0]*f[8])*volFact; 
  out[14] += (0.3535533905932737*zeta2[16]*f[46]+0.3535533905932737*zeta2[5]*f[28]+0.3535533905932737*zeta2[4]*f[25]+0.3535533905932737*zeta2[0]*f[9])*volFact; 
  out[15] += (0.3535533905932737*zeta2[16]*f[57]+0.3535533905932737*zeta2[5]*f[43]+0.3535533905932737*zeta2[4]*f[42]+0.3535533905932737*zeta2[0]*f[22])*volFact; 
  out[16] += (0.3535533905932737*f[20]*zeta3[20]+0.3535533905932737*f[6]*zeta3[6]+0.3535533905932737*f[4]*zeta3[4]+0.3535533905932737*f[0]*zeta3[0])*volFact; 
  out[17] += (0.3535533905932737*zeta3[20]*f[35]+0.3535533905932737*zeta3[6]*f[17]+0.3535533905932737*zeta3[4]*f[10]+0.3535533905932737*zeta3[0]*f[1])*volFact; 
  out[18] += (0.3535533905932737*zeta3[20]*f[36]+0.3535533905932737*zeta3[6]*f[18]+0.3535533905932737*zeta3[4]*f[11]+0.3535533905932737*zeta3[0]*f[2])*volFact; 
  out[19] += (0.3535533905932737*zeta3[20]*f[37]+0.3535533905932737*zeta3[6]*f[19]+0.3535533905932737*zeta3[4]*f[12]+0.3535533905932737*zeta3[0]*f[3])*volFact; 
  out[20] += (0.3535533905932737*zeta3[20]*f[48]+0.3535533905932737*zeta3[6]*f[32]+0.3535533905932737*zeta3[4]*f[23]+0.3535533905932737*zeta3[0]*f[7])*volFact; 
  out[21] += (0.3535533905932737*zeta3[20]*f[49]+0.3535533905932737*zeta3[6]*f[33]+0.3535533905932737*zeta3[4]*f[24]+0.3535533905932737*zeta3[0]*f[8])*volFact; 
  out[22] += (0.3535533905932737*zeta3[20]*f[50]+0.3535533905932737*zeta3[6]*f[34]+0.3535533905932737*zeta3[4]*f[25]+0.3535533905932737*zeta3[0]*f[9])*volFact; 
  out[23] += (0.3535533905932737*zeta3[20]*f[58]+0.3535533905932737*zeta3[6]*f[47]+0.3535533905932737*zeta3[4]*f[42]+0.3535533905932737*zeta3[0]*f[22])*volFact; 
  out[24] += (0.3535533905932737*f[5]*zeta4[5]+0.3535533905932737*f[0]*zeta4[0])*volFact; 
  out[25] += (0.3535533905932737*zeta4[5]*f[13]+0.3535533905932737*zeta4[0]*f[1])*volFact; 
  out[26] += (0.3535533905932737*zeta4[5]*f[14]+0.3535533905932737*zeta4[0]*f[2])*volFact; 
  out[27] += (0.3535533905932737*zeta4[5]*f[15]+0.3535533905932737*zeta4[0]*f[3])*volFact; 
  out[28] += (0.3535533905932737*zeta4[5]*f[26]+0.3535533905932737*zeta4[0]*f[7])*volFact; 
  out[29] += (0.3535533905932737*zeta4[5]*f[27]+0.3535533905932737*zeta4[0]*f[8])*volFact; 
  out[30] += (0.3535533905932737*zeta4[5]*f[28]+0.3535533905932737*zeta4[0]*f[9])*volFact; 
  out[31] += (0.3535533905932737*zeta4[5]*f[43]+0.3535533905932737*zeta4[0]*f[22])*volFact; 
  out[32] += (0.3535533905932737*f[21]*zeta5[21]+0.3535533905932737*f[6]*zeta5[6]+0.3535533905932737*f[5]*zeta5[5]+0.3535533905932737*f[0]*zeta5[0])*volFact; 
  out[33] += (0.3535533905932737*zeta5[21]*f[38]+0.3535533905932737*zeta5[6]*f[17]+0.3535533905932737*zeta5[5]*f[13]+0.3535533905932737*zeta5[0]*f[1])*volFact; 
  out[34] += (0.3535533905932737*zeta5[21]*f[39]+0.3535533905932737*zeta5[6]*f[18]+0.3535533905932737*zeta5[5]*f[14]+0.3535533905932737*zeta5[0]*f[2])*volFact; 
  out[35] += (0.3535533905932737*zeta5[21]*f[40]+0.3535533905932737*zeta5[6]*f[19]+0.3535533905932737*zeta5[5]*f[15]+0.3535533905932737*zeta5[0]*f[3])*volFact; 
  out[36] += (0.3535533905932737*zeta5[21]*f[51]+0.3535533905932737*zeta5[6]*f[32]+0.3535533905932737*zeta5[5]*f[26]+0.3535533905932737*zeta5[0]*f[7])*volFact; 
  out[37] += (0.3535533905932737*zeta5[21]*f[52]+0.3535533905932737*zeta5[6]*f[33]+0.3535533905932737*zeta5[5]*f[27]+0.3535533905932737*zeta5[0]*f[8])*volFact; 
  out[38] += (0.3535533905932737*zeta5[21]*f[53]+0.3535533905932737*zeta5[6]*f[34]+0.3535533905932737*zeta5[5]*f[28]+0.3535533905932737*zeta5[0]*f[9])*volFact; 
  out[39] += (0.3535533905932737*zeta5[21]*f[59]+0.3535533905932737*zeta5[6]*f[47]+0.3535533905932737*zeta5[5]*f[43]+0.3535533905932737*zeta5[0]*f[22])*volFact; 
  out[40] += (0.3535533905932737*f[6]*zeta6[6]+0.3535533905932737*f[0]*zeta6[0])*volFact; 
  out[41] += (0.3535533905932737*zeta6[6]*f[17]+0.3535533905932737*zeta6[0]*f[1])*volFact; 
  out[42] += (0.3535533905932737*zeta6[6]*f[18]+0.3535533905932737*zeta6[0]*f[2])*volFact; 
  out[43] += (0.3535533905932737*zeta6[6]*f[19]+0.3535533905932737*zeta6[0]*f[3])*volFact; 
  out[44] += (0.3535533905932737*zeta6[6]*f[32]+0.3535533905932737*zeta6[0]*f[7])*volFact; 
  out[45] += (0.3535533905932737*zeta6[6]*f[33]+0.3535533905932737*zeta6[0]*f[8])*volFact; 
  out[46] += (0.3535533905932737*zeta6[6]*f[34]+0.3535533905932737*zeta6[0]*f[9])*volFact; 
  out[47] += (0.3535533905932737*zeta6[6]*f[47]+0.3535533905932737*zeta6[0]*f[22])*volFact; 
} 
void MomentCalc3x3vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  double zeta[64]; 

  zeta[0] = 8.0*wx3_sq+8.0*wx2_sq+8.0*wx1_sq+0.6666666666666666*dv3_sq+0.6666666666666666*dv2_sq+0.6666666666666666*dv1_sq; 
  zeta[4] = 4.618802153517007*dv1*wx1; 
  zeta[5] = 4.618802153517007*dv2*wx2; 
  zeta[6] = 4.618802153517007*dv3*wx3; 
  out[0] += (0.3535533905932737*f[6]*zeta[6]+0.3535533905932737*f[5]*zeta[5]+0.3535533905932737*f[4]*zeta[4]+0.3535533905932737*f[0]*zeta[0])*volFact; 
  out[1] += (0.3535533905932737*zeta[6]*f[17]+0.3535533905932737*zeta[5]*f[13]+0.3535533905932737*zeta[4]*f[10]+0.3535533905932737*zeta[0]*f[1])*volFact; 
  out[2] += (0.3535533905932737*zeta[6]*f[18]+0.3535533905932737*zeta[5]*f[14]+0.3535533905932737*zeta[4]*f[11]+0.3535533905932737*zeta[0]*f[2])*volFact; 
  out[3] += (0.3535533905932737*zeta[6]*f[19]+0.3535533905932737*zeta[5]*f[15]+0.3535533905932737*zeta[4]*f[12]+0.3535533905932737*zeta[0]*f[3])*volFact; 
  out[4] += (0.3535533905932737*zeta[6]*f[32]+0.3535533905932737*zeta[5]*f[26]+0.3535533905932737*zeta[4]*f[23]+0.3535533905932737*zeta[0]*f[7])*volFact; 
  out[5] += (0.3535533905932737*zeta[6]*f[33]+0.3535533905932737*zeta[5]*f[27]+0.3535533905932737*zeta[4]*f[24]+0.3535533905932737*zeta[0]*f[8])*volFact; 
  out[6] += (0.3535533905932737*zeta[6]*f[34]+0.3535533905932737*zeta[5]*f[28]+0.3535533905932737*zeta[4]*f[25]+0.3535533905932737*zeta[0]*f[9])*volFact; 
  out[7] += (0.3535533905932737*zeta[6]*f[47]+0.3535533905932737*zeta[5]*f[43]+0.3535533905932737*zeta[4]*f[42]+0.3535533905932737*zeta[0]*f[22])*volFact; 
} 
void MomentCalc3x3vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  const double wx3_cu = wx3*wx3*wx3, dv3_cu = dv3*dv3*dv3; 
  double zeta1[64]; 

  zeta1[0] = 8.0*wx1*wx3_sq+8.0*wx1*wx2_sq+8.0*wx1_cu+0.6666666666666666*dv3_sq*wx1+0.6666666666666666*dv2_sq*wx1+2.0*dv1_sq*wx1; 
  zeta1[4] = 2.309401076758503*dv1*wx3_sq+2.309401076758503*dv1*wx2_sq+6.928203230275509*dv1*wx1_sq+0.1924500897298753*dv1*dv3_sq+0.1924500897298753*dv1*dv2_sq+0.3464101615137755*dv1_cu; 
  zeta1[5] = 4.618802153517007*dv2*wx1*wx2; 
  zeta1[6] = 4.618802153517007*dv3*wx1*wx3; 
  zeta1[16] = 1.333333333333333*dv1*dv2*wx2; 
  zeta1[20] = 1.333333333333333*dv1*dv3*wx3; 
  double zeta2[64]; 

  zeta2[0] = 8.0*wx2*wx3_sq+8.0*wx2_cu+8.0*wx1_sq*wx2+0.6666666666666666*dv3_sq*wx2+2.0*dv2_sq*wx2+0.6666666666666666*dv1_sq*wx2; 
  zeta2[4] = 4.618802153517007*dv1*wx1*wx2; 
  zeta2[5] = 2.309401076758503*dv2*wx3_sq+6.928203230275509*dv2*wx2_sq+2.309401076758503*dv2*wx1_sq+0.1924500897298753*dv2*dv3_sq+0.3464101615137755*dv2_cu+0.1924500897298753*dv1_sq*dv2; 
  zeta2[6] = 4.618802153517007*dv3*wx2*wx3; 
  zeta2[16] = 1.333333333333333*dv1*dv2*wx1; 
  zeta2[21] = 1.333333333333333*dv2*dv3*wx3; 
  double zeta3[64]; 

  zeta3[0] = 8.0*wx3_cu+8.0*wx2_sq*wx3+8.0*wx1_sq*wx3+2.0*dv3_sq*wx3+0.6666666666666666*dv2_sq*wx3+0.6666666666666666*dv1_sq*wx3; 
  zeta3[4] = 4.618802153517007*dv1*wx1*wx3; 
  zeta3[5] = 4.618802153517007*dv2*wx2*wx3; 
  zeta3[6] = 6.928203230275509*dv3*wx3_sq+2.309401076758503*dv3*wx2_sq+2.309401076758503*dv3*wx1_sq+0.3464101615137755*dv3_cu+0.1924500897298753*dv2_sq*dv3+0.1924500897298753*dv1_sq*dv3; 
  zeta3[20] = 1.333333333333333*dv1*dv3*wx1; 
  zeta3[21] = 1.333333333333333*dv2*dv3*wx2; 
  out[0] += (0.3535533905932737*f[20]*zeta1[20]+0.3535533905932737*f[16]*zeta1[16]+0.3535533905932737*f[6]*zeta1[6]+0.3535533905932737*f[5]*zeta1[5]+0.3535533905932737*f[4]*zeta1[4]+0.3535533905932737*f[0]*zeta1[0])*volFact; 
  out[1] += (0.3535533905932737*zeta1[20]*f[35]+0.3535533905932737*zeta1[16]*f[29]+0.3535533905932737*zeta1[6]*f[17]+0.3535533905932737*zeta1[5]*f[13]+0.3535533905932737*zeta1[4]*f[10]+0.3535533905932737*zeta1[0]*f[1])*volFact; 
  out[2] += (0.3535533905932737*zeta1[20]*f[36]+0.3535533905932737*zeta1[16]*f[30]+0.3535533905932737*zeta1[6]*f[18]+0.3535533905932737*zeta1[5]*f[14]+0.3535533905932737*zeta1[4]*f[11]+0.3535533905932737*zeta1[0]*f[2])*volFact; 
  out[3] += (0.3535533905932737*zeta1[20]*f[37]+0.3535533905932737*zeta1[16]*f[31]+0.3535533905932737*zeta1[6]*f[19]+0.3535533905932737*zeta1[5]*f[15]+0.3535533905932737*zeta1[4]*f[12]+0.3535533905932737*zeta1[0]*f[3])*volFact; 
  out[4] += (0.3535533905932737*zeta1[20]*f[48]+0.3535533905932737*zeta1[16]*f[44]+0.3535533905932737*zeta1[6]*f[32]+0.3535533905932737*zeta1[5]*f[26]+0.3535533905932737*zeta1[4]*f[23]+0.3535533905932737*zeta1[0]*f[7])*volFact; 
  out[5] += (0.3535533905932737*zeta1[20]*f[49]+0.3535533905932737*zeta1[16]*f[45]+0.3535533905932737*zeta1[6]*f[33]+0.3535533905932737*zeta1[5]*f[27]+0.3535533905932737*zeta1[4]*f[24]+0.3535533905932737*zeta1[0]*f[8])*volFact; 
  out[6] += (0.3535533905932737*zeta1[20]*f[50]+0.3535533905932737*zeta1[16]*f[46]+0.3535533905932737*zeta1[6]*f[34]+0.3535533905932737*zeta1[5]*f[28]+0.3535533905932737*zeta1[4]*f[25]+0.3535533905932737*zeta1[0]*f[9])*volFact; 
  out[7] += (0.3535533905932737*zeta1[20]*f[58]+0.3535533905932737*zeta1[16]*f[57]+0.3535533905932737*zeta1[6]*f[47]+0.3535533905932737*zeta1[5]*f[43]+0.3535533905932737*zeta1[4]*f[42]+0.3535533905932737*zeta1[0]*f[22])*volFact; 
  out[8] += (0.3535533905932737*f[21]*zeta2[21]+0.3535533905932737*f[16]*zeta2[16]+0.3535533905932737*f[6]*zeta2[6]+0.3535533905932737*f[5]*zeta2[5]+0.3535533905932737*f[4]*zeta2[4]+0.3535533905932737*f[0]*zeta2[0])*volFact; 
  out[9] += (0.3535533905932737*zeta2[21]*f[38]+0.3535533905932737*zeta2[16]*f[29]+0.3535533905932737*zeta2[6]*f[17]+0.3535533905932737*zeta2[5]*f[13]+0.3535533905932737*zeta2[4]*f[10]+0.3535533905932737*zeta2[0]*f[1])*volFact; 
  out[10] += (0.3535533905932737*zeta2[21]*f[39]+0.3535533905932737*zeta2[16]*f[30]+0.3535533905932737*zeta2[6]*f[18]+0.3535533905932737*zeta2[5]*f[14]+0.3535533905932737*zeta2[4]*f[11]+0.3535533905932737*zeta2[0]*f[2])*volFact; 
  out[11] += (0.3535533905932737*zeta2[21]*f[40]+0.3535533905932737*zeta2[16]*f[31]+0.3535533905932737*zeta2[6]*f[19]+0.3535533905932737*zeta2[5]*f[15]+0.3535533905932737*zeta2[4]*f[12]+0.3535533905932737*zeta2[0]*f[3])*volFact; 
  out[12] += (0.3535533905932737*zeta2[21]*f[51]+0.3535533905932737*zeta2[16]*f[44]+0.3535533905932737*zeta2[6]*f[32]+0.3535533905932737*zeta2[5]*f[26]+0.3535533905932737*zeta2[4]*f[23]+0.3535533905932737*zeta2[0]*f[7])*volFact; 
  out[13] += (0.3535533905932737*zeta2[21]*f[52]+0.3535533905932737*zeta2[16]*f[45]+0.3535533905932737*zeta2[6]*f[33]+0.3535533905932737*zeta2[5]*f[27]+0.3535533905932737*zeta2[4]*f[24]+0.3535533905932737*zeta2[0]*f[8])*volFact; 
  out[14] += (0.3535533905932737*zeta2[21]*f[53]+0.3535533905932737*zeta2[16]*f[46]+0.3535533905932737*zeta2[6]*f[34]+0.3535533905932737*zeta2[5]*f[28]+0.3535533905932737*zeta2[4]*f[25]+0.3535533905932737*zeta2[0]*f[9])*volFact; 
  out[15] += (0.3535533905932737*zeta2[21]*f[59]+0.3535533905932737*zeta2[16]*f[57]+0.3535533905932737*zeta2[6]*f[47]+0.3535533905932737*zeta2[5]*f[43]+0.3535533905932737*zeta2[4]*f[42]+0.3535533905932737*zeta2[0]*f[22])*volFact; 
  out[16] += (0.3535533905932737*f[21]*zeta3[21]+0.3535533905932737*f[20]*zeta3[20]+0.3535533905932737*f[6]*zeta3[6]+0.3535533905932737*f[5]*zeta3[5]+0.3535533905932737*f[4]*zeta3[4]+0.3535533905932737*f[0]*zeta3[0])*volFact; 
  out[17] += (0.3535533905932737*zeta3[21]*f[38]+0.3535533905932737*zeta3[20]*f[35]+0.3535533905932737*zeta3[6]*f[17]+0.3535533905932737*zeta3[5]*f[13]+0.3535533905932737*zeta3[4]*f[10]+0.3535533905932737*zeta3[0]*f[1])*volFact; 
  out[18] += (0.3535533905932737*zeta3[21]*f[39]+0.3535533905932737*zeta3[20]*f[36]+0.3535533905932737*zeta3[6]*f[18]+0.3535533905932737*zeta3[5]*f[14]+0.3535533905932737*zeta3[4]*f[11]+0.3535533905932737*zeta3[0]*f[2])*volFact; 
  out[19] += (0.3535533905932737*zeta3[21]*f[40]+0.3535533905932737*zeta3[20]*f[37]+0.3535533905932737*zeta3[6]*f[19]+0.3535533905932737*zeta3[5]*f[15]+0.3535533905932737*zeta3[4]*f[12]+0.3535533905932737*zeta3[0]*f[3])*volFact; 
  out[20] += (0.3535533905932737*zeta3[21]*f[51]+0.3535533905932737*zeta3[20]*f[48]+0.3535533905932737*zeta3[6]*f[32]+0.3535533905932737*zeta3[5]*f[26]+0.3535533905932737*zeta3[4]*f[23]+0.3535533905932737*zeta3[0]*f[7])*volFact; 
  out[21] += (0.3535533905932737*zeta3[21]*f[52]+0.3535533905932737*zeta3[20]*f[49]+0.3535533905932737*zeta3[6]*f[33]+0.3535533905932737*zeta3[5]*f[27]+0.3535533905932737*zeta3[4]*f[24]+0.3535533905932737*zeta3[0]*f[8])*volFact; 
  out[22] += (0.3535533905932737*zeta3[21]*f[53]+0.3535533905932737*zeta3[20]*f[50]+0.3535533905932737*zeta3[6]*f[34]+0.3535533905932737*zeta3[5]*f[28]+0.3535533905932737*zeta3[4]*f[25]+0.3535533905932737*zeta3[0]*f[9])*volFact; 
  out[23] += (0.3535533905932737*zeta3[21]*f[59]+0.3535533905932737*zeta3[20]*f[58]+0.3535533905932737*zeta3[6]*f[47]+0.3535533905932737*zeta3[5]*f[43]+0.3535533905932737*zeta3[4]*f[42]+0.3535533905932737*zeta3[0]*f[22])*volFact; 
} 
void MomentCalc3x3vSer_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  double tempM0[8], tempM1i[24]; 

  double zeta[64]; 

  zeta[0] = 8.0*wx3_sq+8.0*wx2_sq+8.0*wx1_sq+0.6666666666666666*dv3_sq+0.6666666666666666*dv2_sq+0.6666666666666666*dv1_sq; 
  zeta[4] = 4.618802153517007*dv1*wx1; 
  zeta[5] = 4.618802153517007*dv2*wx2; 
  zeta[6] = 4.618802153517007*dv3*wx3; 
  tempM0[0] = 2.828427124746191*f[0]*volFact; 
  tempM0[1] = 2.828427124746191*f[1]*volFact; 
  tempM0[2] = 2.828427124746191*f[2]*volFact; 
  tempM0[3] = 2.828427124746191*f[3]*volFact; 
  tempM0[4] = 2.828427124746191*f[7]*volFact; 
  tempM0[5] = 2.828427124746191*f[8]*volFact; 
  tempM0[6] = 2.828427124746191*f[9]*volFact; 
  tempM0[7] = 2.828427124746191*f[22]*volFact; 

  tempM1i[0] = 0.8164965809277261*f[4]*dv1*volFact+tempM0[0]*wx1; 
  tempM1i[1] = 0.8164965809277261*f[10]*dv1*volFact+tempM0[1]*wx1; 
  tempM1i[2] = 0.8164965809277261*f[11]*dv1*volFact+tempM0[2]*wx1; 
  tempM1i[3] = 0.8164965809277261*f[12]*dv1*volFact+tempM0[3]*wx1; 
  tempM1i[4] = 0.8164965809277261*f[23]*dv1*volFact+tempM0[4]*wx1; 
  tempM1i[5] = 0.8164965809277261*f[24]*dv1*volFact+tempM0[5]*wx1; 
  tempM1i[6] = 0.8164965809277261*f[25]*dv1*volFact+tempM0[6]*wx1; 
  tempM1i[7] = 0.8164965809277261*f[42]*dv1*volFact+tempM0[7]*wx1; 
  tempM1i[8] = 0.8164965809277261*f[5]*dv2*volFact+tempM0[0]*wx2; 
  tempM1i[9] = 0.8164965809277261*f[13]*dv2*volFact+tempM0[1]*wx2; 
  tempM1i[10] = 0.8164965809277261*f[14]*dv2*volFact+tempM0[2]*wx2; 
  tempM1i[11] = 0.8164965809277261*f[15]*dv2*volFact+tempM0[3]*wx2; 
  tempM1i[12] = 0.8164965809277261*f[26]*dv2*volFact+tempM0[4]*wx2; 
  tempM1i[13] = 0.8164965809277261*f[27]*dv2*volFact+tempM0[5]*wx2; 
  tempM1i[14] = 0.8164965809277261*f[28]*dv2*volFact+tempM0[6]*wx2; 
  tempM1i[15] = 0.8164965809277261*f[43]*dv2*volFact+tempM0[7]*wx2; 
  tempM1i[16] = 0.8164965809277261*f[6]*dv3*volFact+tempM0[0]*wx3; 
  tempM1i[17] = 0.8164965809277261*f[17]*dv3*volFact+tempM0[1]*wx3; 
  tempM1i[18] = 0.8164965809277261*f[18]*dv3*volFact+tempM0[2]*wx3; 
  tempM1i[19] = 0.8164965809277261*f[19]*dv3*volFact+tempM0[3]*wx3; 
  tempM1i[20] = 0.8164965809277261*f[32]*dv3*volFact+tempM0[4]*wx3; 
  tempM1i[21] = 0.8164965809277261*f[33]*dv3*volFact+tempM0[5]*wx3; 
  tempM1i[22] = 0.8164965809277261*f[34]*dv3*volFact+tempM0[6]*wx3; 
  tempM1i[23] = 0.8164965809277261*f[47]*dv3*volFact+tempM0[7]*wx3; 

  outM0[0] += tempM0[0]; 
  outM0[1] += tempM0[1]; 
  outM0[2] += tempM0[2]; 
  outM0[3] += tempM0[3]; 
  outM0[4] += tempM0[4]; 
  outM0[5] += tempM0[5]; 
  outM0[6] += tempM0[6]; 
  outM0[7] += tempM0[7]; 
  outM1i[0] += tempM1i[0]; 
  outM1i[1] += tempM1i[1]; 
  outM1i[2] += tempM1i[2]; 
  outM1i[3] += tempM1i[3]; 
  outM1i[4] += tempM1i[4]; 
  outM1i[5] += tempM1i[5]; 
  outM1i[6] += tempM1i[6]; 
  outM1i[7] += tempM1i[7]; 
  outM1i[8] += tempM1i[8]; 
  outM1i[9] += tempM1i[9]; 
  outM1i[10] += tempM1i[10]; 
  outM1i[11] += tempM1i[11]; 
  outM1i[12] += tempM1i[12]; 
  outM1i[13] += tempM1i[13]; 
  outM1i[14] += tempM1i[14]; 
  outM1i[15] += tempM1i[15]; 
  outM1i[16] += tempM1i[16]; 
  outM1i[17] += tempM1i[17]; 
  outM1i[18] += tempM1i[18]; 
  outM1i[19] += tempM1i[19]; 
  outM1i[20] += tempM1i[20]; 
  outM1i[21] += tempM1i[21]; 
  outM1i[22] += tempM1i[22]; 
  outM1i[23] += tempM1i[23]; 
  outM2[0] += (0.3535533905932737*f[6]*zeta[6]+0.3535533905932737*f[5]*zeta[5]+0.3535533905932737*f[4]*zeta[4]+0.3535533905932737*f[0]*zeta[0])*volFact; 
  outM2[1] += (0.3535533905932737*zeta[6]*f[17]+0.3535533905932737*zeta[5]*f[13]+0.3535533905932737*zeta[4]*f[10]+0.3535533905932737*zeta[0]*f[1])*volFact; 
  outM2[2] += (0.3535533905932737*zeta[6]*f[18]+0.3535533905932737*zeta[5]*f[14]+0.3535533905932737*zeta[4]*f[11]+0.3535533905932737*zeta[0]*f[2])*volFact; 
  outM2[3] += (0.3535533905932737*zeta[6]*f[19]+0.3535533905932737*zeta[5]*f[15]+0.3535533905932737*zeta[4]*f[12]+0.3535533905932737*zeta[0]*f[3])*volFact; 
  outM2[4] += (0.3535533905932737*zeta[6]*f[32]+0.3535533905932737*zeta[5]*f[26]+0.3535533905932737*zeta[4]*f[23]+0.3535533905932737*zeta[0]*f[7])*volFact; 
  outM2[5] += (0.3535533905932737*zeta[6]*f[33]+0.3535533905932737*zeta[5]*f[27]+0.3535533905932737*zeta[4]*f[24]+0.3535533905932737*zeta[0]*f[8])*volFact; 
  outM2[6] += (0.3535533905932737*zeta[6]*f[34]+0.3535533905932737*zeta[5]*f[28]+0.3535533905932737*zeta[4]*f[25]+0.3535533905932737*zeta[0]*f[9])*volFact; 
  outM2[7] += (0.3535533905932737*zeta[6]*f[47]+0.3535533905932737*zeta[5]*f[43]+0.3535533905932737*zeta[4]*f[42]+0.3535533905932737*zeta[0]*f[22])*volFact; 
} 
