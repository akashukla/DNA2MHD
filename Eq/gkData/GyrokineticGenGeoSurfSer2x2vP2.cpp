#include <GyrokineticModDecl.h>
double GyrokineticGenGeoSurf2x2vSer_x_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[11] = 2.0*phi[4]*q_; 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 

  double BstarXdBmagR[48]; 

  double alphaR[20]; 
  alphaR[0] = -(0.1767766952966368*b_z[0]*jacobTotInv[0]*(3.872983346207417*hamilR[19]-3.0*hamilR[5]+1.732050807568877*hamilR[2])*rdy2R)/q_; 
  alphaR[1] = (0.1767766952966368*b_z[0]*jacobTotInv[0]*(6.708203932499369*hamilR[20]-3.872983346207417*hamilR[12])*rdy2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(3.872983346207417*b_z[0]*jacobTotInv[0]*hamilR[19]-3.0*b_z[0]*jacobTotInv[0]*hamilR[5]+1.732050807568877*b_z[0]*jacobTotInv[0]*hamilR[2])*rdy2R)/q_; 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[11]+fL[11])+1.732050807568877*(fL[1]-1.0*fR[1])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[19]+fL[19])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[2]+fL[2])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[21]+fL[21])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[25]+fL[25])+3.0*(fL[8]-1.0*fR[8]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[32]+fL[32])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[7]+fL[7]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[35]+fL[35])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[9]+fL[9]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[37]+fL[37])+1.732050807568877*(fL[17]-1.0*fR[17])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[20]-1.0*(8.660254037844387*fL[20]+5.0*(fR[12]+fL[12]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[23]-1.0*(8.660254037844387*fL[23]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[28]-1.0*(8.660254037844387*fL[28]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[44]+fL[44])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[18]+fL[18])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[33]-1.0*(8.660254037844387*fL[33]+5.0*(fR[22]+fL[22]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[24]+fL[24]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[36]-1.0*(8.660254037844387*fL[36]+5.0*(fR[26]+fL[26]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[39]-1.0*(8.660254037844387*fL[39]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[29]+fL[29]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[42]-1.0*(8.660254037844387*fL[42]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[45]-1.0*(8.660254037844387*fL[45]+5.0*(fR[38]+fL[38]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[40]+fL[40]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[43]+fL[43]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[11]-6.324555320336761*fL[11]-4.898979485566357*(fR[1]+fL[1])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[19]-94.8683298050514*fL[19]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[2]-42.42640687119286*fL[2])*amax-18.97366596101028*alphaR[1]*fAvg[7]-21.21320343559643*(alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.04166666666666666*((18.97366596101028*fR[21]-18.97366596101028*fL[21]-14.69693845669907*(fR[6]+fL[6])+8.485281374238571*fR[3]-8.485281374238571*fL[3])*amax-4.242640687119286*(alphaR[1]*fAvg[4]+alphaR[0]*fAvg[2])); 
  Ghat[3] = -0.04166666666666666*((18.97366596101028*fR[25]-18.97366596101028*fL[25]-14.69693845669907*(fR[8]+fL[8])+8.485281374238571*fR[4]-8.485281374238571*fL[4])*amax-4.242640687119286*(alphaR[1]*fAvg[5]+alphaR[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[32]-94.86832980505142*fL[32]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[7]-42.42640687119286*fL[7])*amax-18.97366596101028*alphaR[1]*fAvg[11]-21.21320343559643*(alphaR[0]*fAvg[4]+alphaR[1]*fAvg[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[35]-94.86832980505142*fL[35]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[9]-42.42640687119286*fL[9])*amax-18.97366596101028*alphaR[1]*fAvg[13]-21.21320343559643*(alphaR[0]*fAvg[5]+alphaR[1]*fAvg[3])); 
  Ghat[6] = -0.125*((6.324555320336761*fR[37]-6.324555320336761*fL[37]-4.898979485566357*(fR[17]+fL[17])+2.828427124746191*fR[10]-2.828427124746191*fL[10])*amax-1.414213562373095*(alphaR[1]*fAvg[10]+alphaR[0]*fAvg[6])); 
  Ghat[7] = 0.025*((24.49489742783179*(fR[20]+fL[20])-14.14213562373095*fR[12]+14.14213562373095*fL[12])*amax+7.071067811865476*alphaR[0]*fAvg[7]+6.324555320336761*alphaR[1]*fAvg[1]); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[23]+fL[23])-42.42640687119286*fR[13]+42.42640687119286*fL[13])*amax+21.21320343559643*alphaR[1]*fAvg[12]+21.21320343559643*alphaR[0]*fAvg[8]); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[28]+fL[28])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[1]*fAvg[15]+21.21320343559643*alphaR[0]*fAvg[9]); 
  Ghat[10] = -0.008333333333333333*((94.8683298050514*fR[44]-94.8683298050514*fL[44]-73.48469228349535*(fR[31]+fL[31])+42.42640687119286*fR[18]-42.42640687119286*fL[18])*amax-18.97366596101028*alphaR[1]*fAvg[17]-21.21320343559643*(alphaR[0]*fAvg[10]+alphaR[1]*fAvg[6])); 
  Ghat[11] = 0.008333333333333333*((73.48469228349536*(fR[33]+fL[33])-42.42640687119286*fR[22]+42.42640687119286*fL[22])*amax+21.21320343559643*alphaR[0]*fAvg[11]+18.97366596101028*alphaR[1]*fAvg[4]); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[34]+fL[34])-42.42640687119286*fR[24]+42.42640687119286*fL[24])*amax+21.21320343559643*alphaR[0]*fAvg[12]+21.21320343559643*alphaR[1]*fAvg[8]); 
  Ghat[13] = 0.008333333333333333*((73.48469228349536*(fR[36]+fL[36])-42.42640687119286*fR[26]+42.42640687119286*fL[26])*amax+21.21320343559643*alphaR[0]*fAvg[13]+18.97366596101028*alphaR[1]*fAvg[5]); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[39]+fL[39])-42.42640687119286*fR[27]+42.42640687119286*fL[27])*amax+21.21320343559643*alphaR[1]*fAvg[18]+21.21320343559643*alphaR[0]*fAvg[14]); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[29]+42.42640687119286*fL[29])*amax+21.21320343559643*alphaR[0]*fAvg[15]+21.21320343559643*alphaR[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[42]+fL[42])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+21.21320343559643*alphaR[1]*fAvg[19]+21.21320343559643*alphaR[0]*fAvg[16]); 
  Ghat[17] = 0.025*((24.49489742783179*(fR[45]+fL[45])-14.14213562373095*fR[38]+14.14213562373095*fL[38])*amax+7.071067811865476*alphaR[0]*fAvg[17]+6.324555320336761*alphaR[1]*fAvg[10]); 
  Ghat[18] = 0.008333333333333333*((73.48469228349536*(fR[46]+fL[46])-42.42640687119286*fR[40]+42.42640687119286*fL[40])*amax+21.21320343559643*alphaR[0]*fAvg[18]+21.21320343559643*alphaR[1]*fAvg[14]); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[43]+42.42640687119286*fL[43])*amax+21.21320343559643*alphaR[0]*fAvg[19]+21.21320343559643*alphaR[1]*fAvg[16]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = -1.224744871391589*Ghat[0]; 
  incr[2] = 0.7071067811865475*Ghat[1]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = -1.224744871391589*Ghat[2]; 
  incr[7] = 0.7071067811865475*Ghat[4]; 
  incr[8] = -1.224744871391589*Ghat[3]; 
  incr[9] = 0.7071067811865475*Ghat[5]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 1.58113883008419*Ghat[0]; 
  incr[12] = 0.7071067811865475*Ghat[7]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = -1.224744871391589*Ghat[6]; 
  incr[18] = 0.7071067811865475*Ghat[10]; 
  incr[19] = 1.58113883008419*Ghat[1]; 
  incr[20] = -1.224744871391589*Ghat[7]; 
  incr[21] = 1.58113883008419*Ghat[2]; 
  incr[22] = 0.7071067811865475*Ghat[11]; 
  incr[23] = -1.224744871391589*Ghat[8]; 
  incr[24] = 0.7071067811865475*Ghat[12]; 
  incr[25] = 1.58113883008419*Ghat[3]; 
  incr[26] = 0.7071067811865475*Ghat[13]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = -1.224744871391589*Ghat[9]; 
  incr[29] = 0.7071067811865475*Ghat[15]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = 1.58113883008419*Ghat[4]; 
  incr[33] = -1.224744871391589*Ghat[11]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = 1.58113883008419*Ghat[5]; 
  incr[36] = -1.224744871391589*Ghat[13]; 
  incr[37] = 1.58113883008419*Ghat[6]; 
  incr[38] = 0.7071067811865475*Ghat[17]; 
  incr[39] = -1.224744871391589*Ghat[14]; 
  incr[40] = 0.7071067811865475*Ghat[18]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = -1.224744871391589*Ghat[16]; 
  incr[43] = 0.7071067811865475*Ghat[19]; 
  incr[44] = 1.58113883008419*Ghat[10]; 
  incr[45] = -1.224744871391589*Ghat[17]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 
  outR[8] += incr[8]*rdx2R; 
  outR[9] += incr[9]*rdx2R; 
  outR[10] += incr[10]*rdx2R; 
  outR[11] += incr[11]*rdx2R; 
  outR[12] += incr[12]*rdx2R; 
  outR[13] += incr[13]*rdx2R; 
  outR[14] += incr[14]*rdx2R; 
  outR[15] += incr[15]*rdx2R; 
  outR[16] += incr[16]*rdx2R; 
  outR[17] += incr[17]*rdx2R; 
  outR[18] += incr[18]*rdx2R; 
  outR[19] += incr[19]*rdx2R; 
  outR[20] += incr[20]*rdx2R; 
  outR[21] += incr[21]*rdx2R; 
  outR[22] += incr[22]*rdx2R; 
  outR[23] += incr[23]*rdx2R; 
  outR[24] += incr[24]*rdx2R; 
  outR[25] += incr[25]*rdx2R; 
  outR[26] += incr[26]*rdx2R; 
  outR[27] += incr[27]*rdx2R; 
  outR[28] += incr[28]*rdx2R; 
  outR[29] += incr[29]*rdx2R; 
  outR[30] += incr[30]*rdx2R; 
  outR[31] += incr[31]*rdx2R; 
  outR[32] += incr[32]*rdx2R; 
  outR[33] += incr[33]*rdx2R; 
  outR[34] += incr[34]*rdx2R; 
  outR[35] += incr[35]*rdx2R; 
  outR[36] += incr[36]*rdx2R; 
  outR[37] += incr[37]*rdx2R; 
  outR[38] += incr[38]*rdx2R; 
  outR[39] += incr[39]*rdx2R; 
  outR[40] += incr[40]*rdx2R; 
  outR[41] += incr[41]*rdx2R; 
  outR[42] += incr[42]*rdx2R; 
  outR[43] += incr[43]*rdx2R; 
  outR[44] += incr[44]*rdx2R; 
  outR[45] += incr[45]*rdx2R; 
  outR[46] += incr[46]*rdx2R; 
  outR[47] += incr[47]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += -1.0*incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += incr[6]*rdx2L; 
  outL[7] += -1.0*incr[7]*rdx2L; 
  outL[8] += incr[8]*rdx2L; 
  outL[9] += -1.0*incr[9]*rdx2L; 
  outL[10] += -1.0*incr[10]*rdx2L; 
  outL[11] += -1.0*incr[11]*rdx2L; 
  outL[12] += -1.0*incr[12]*rdx2L; 
  outL[13] += -1.0*incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 
  outL[16] += incr[16]*rdx2L; 
  outL[17] += incr[17]*rdx2L; 
  outL[18] += -1.0*incr[18]*rdx2L; 
  outL[19] += -1.0*incr[19]*rdx2L; 
  outL[20] += incr[20]*rdx2L; 
  outL[21] += -1.0*incr[21]*rdx2L; 
  outL[22] += -1.0*incr[22]*rdx2L; 
  outL[23] += incr[23]*rdx2L; 
  outL[24] += -1.0*incr[24]*rdx2L; 
  outL[25] += -1.0*incr[25]*rdx2L; 
  outL[26] += -1.0*incr[26]*rdx2L; 
  outL[27] += -1.0*incr[27]*rdx2L; 
  outL[28] += incr[28]*rdx2L; 
  outL[29] += -1.0*incr[29]*rdx2L; 
  outL[30] += -1.0*incr[30]*rdx2L; 
  outL[31] += incr[31]*rdx2L; 
  outL[32] += -1.0*incr[32]*rdx2L; 
  outL[33] += incr[33]*rdx2L; 
  outL[34] += incr[34]*rdx2L; 
  outL[35] += -1.0*incr[35]*rdx2L; 
  outL[36] += incr[36]*rdx2L; 
  outL[37] += -1.0*incr[37]*rdx2L; 
  outL[38] += -1.0*incr[38]*rdx2L; 
  outL[39] += incr[39]*rdx2L; 
  outL[40] += -1.0*incr[40]*rdx2L; 
  outL[41] += incr[41]*rdx2L; 
  outL[42] += incr[42]*rdx2L; 
  outL[43] += -1.0*incr[43]*rdx2L; 
  outL[44] += -1.0*incr[44]*rdx2L; 
  outL[45] += incr[45]*rdx2L; 
  outL[46] += incr[46]*rdx2L; 
  outL[47] += incr[47]*rdx2L; 
return std::abs(alphaSurfAvgR); 
} 
double GyrokineticGenGeoSurf2x2vSer_y_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[11] = 2.0*phi[4]*q_; 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 

  double BstarYdBmagR[48]; 

  double alphaR[20]; 
  alphaR[0] = (0.1767766952966368*b_z[0]*jacobTotInv[0]*(3.872983346207417*hamilR[20]-3.0*hamilR[5]+1.732050807568877*hamilR[1])*rdx2R)/q_; 
  alphaR[1] = -(0.1767766952966368*b_z[0]*jacobTotInv[0]*(6.708203932499369*hamilR[19]-3.872983346207417*hamilR[11])*rdx2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*(3.872983346207417*b_z[0]*jacobTotInv[0]*hamilR[20]-3.0*b_z[0]*jacobTotInv[0]*hamilR[5]+1.732050807568877*b_z[0]*jacobTotInv[0]*hamilR[1])*rdx2R)/q_; 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[12]+fL[12])+1.732050807568877*(fL[2]-1.0*fR[2])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[20]+fL[20])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[22]+fL[22])+3.0*(fL[7]-1.0*fR[7]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[26]+fL[26])+3.0*(fL[9]-1.0*fR[9]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[33]+fL[33])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[6]+fL[6]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[36]+fL[36])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[8]+fL[8]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[38]+fL[38])+1.732050807568877*(fL[18]-1.0*fR[18])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[19]-1.0*(8.660254037844387*fL[19]+5.0*(fR[11]+fL[11]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[24]-1.0*(8.660254037844387*fL[24]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[29]-1.0*(8.660254037844387*fL[29]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[45]+fL[45])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[17]+fL[17])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[32]-1.0*(8.660254037844387*fL[32]+5.0*(fR[21]+fL[21]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[23]+fL[23]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[35]-1.0*(8.660254037844387*fL[35]+5.0*(fR[25]+fL[25]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[40]-1.0*(8.660254037844387*fL[40]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[28]+fL[28]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[43]-1.0*(8.660254037844387*fL[43]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[44]-1.0*(8.660254037844387*fL[44]+5.0*(fR[37]+fL[37]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[39]+fL[39]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[42]+fL[42]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[12]-6.324555320336761*fL[12]-4.898979485566357*(fR[2]+fL[2])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[20]-94.8683298050514*fL[20]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[1]-42.42640687119286*fL[1])*amax-18.97366596101028*alphaR[1]*fAvg[7]-21.21320343559643*(alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.04166666666666666*((18.97366596101028*fR[22]-18.97366596101028*fL[22]-14.69693845669907*(fR[7]+fL[7])+8.485281374238571*fR[3]-8.485281374238571*fL[3])*amax-4.242640687119286*(alphaR[1]*fAvg[4]+alphaR[0]*fAvg[2])); 
  Ghat[3] = -0.04166666666666666*((18.97366596101028*fR[26]-18.97366596101028*fL[26]-14.69693845669907*(fR[9]+fL[9])+8.485281374238571*fR[4]-8.485281374238571*fL[4])*amax-4.242640687119286*(alphaR[1]*fAvg[5]+alphaR[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[33]-94.86832980505142*fL[33]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[6]-42.42640687119286*fL[6])*amax-18.97366596101028*alphaR[1]*fAvg[11]-21.21320343559643*(alphaR[0]*fAvg[4]+alphaR[1]*fAvg[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[36]-94.86832980505142*fL[36]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[8]-42.42640687119286*fL[8])*amax-18.97366596101028*alphaR[1]*fAvg[13]-21.21320343559643*(alphaR[0]*fAvg[5]+alphaR[1]*fAvg[3])); 
  Ghat[6] = -0.125*((6.324555320336761*fR[38]-6.324555320336761*fL[38]-4.898979485566357*(fR[18]+fL[18])+2.828427124746191*fR[10]-2.828427124746191*fL[10])*amax-1.414213562373095*(alphaR[1]*fAvg[10]+alphaR[0]*fAvg[6])); 
  Ghat[7] = 0.025*((24.49489742783179*(fR[19]+fL[19])-14.14213562373095*fR[11]+14.14213562373095*fL[11])*amax+7.071067811865476*alphaR[0]*fAvg[7]+6.324555320336761*alphaR[1]*fAvg[1]); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[24]+fL[24])-42.42640687119286*fR[13]+42.42640687119286*fL[13])*amax+21.21320343559643*alphaR[1]*fAvg[12]+21.21320343559643*alphaR[0]*fAvg[8]); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[29]+fL[29])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[1]*fAvg[15]+21.21320343559643*alphaR[0]*fAvg[9]); 
  Ghat[10] = -0.008333333333333333*((94.8683298050514*fR[45]-94.8683298050514*fL[45]-73.48469228349535*(fR[31]+fL[31])+42.42640687119286*fR[17]-42.42640687119286*fL[17])*amax-18.97366596101028*alphaR[1]*fAvg[17]-21.21320343559643*(alphaR[0]*fAvg[10]+alphaR[1]*fAvg[6])); 
  Ghat[11] = 0.008333333333333333*((73.48469228349536*(fR[32]+fL[32])-42.42640687119286*fR[21]+42.42640687119286*fL[21])*amax+21.21320343559643*alphaR[0]*fAvg[11]+18.97366596101028*alphaR[1]*fAvg[4]); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[34]+fL[34])-42.42640687119286*fR[23]+42.42640687119286*fL[23])*amax+21.21320343559643*alphaR[0]*fAvg[12]+21.21320343559643*alphaR[1]*fAvg[8]); 
  Ghat[13] = 0.008333333333333333*((73.48469228349536*(fR[35]+fL[35])-42.42640687119286*fR[25]+42.42640687119286*fL[25])*amax+21.21320343559643*alphaR[0]*fAvg[13]+18.97366596101028*alphaR[1]*fAvg[5]); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[40]+fL[40])-42.42640687119286*fR[27]+42.42640687119286*fL[27])*amax+21.21320343559643*alphaR[1]*fAvg[18]+21.21320343559643*alphaR[0]*fAvg[14]); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[28]+42.42640687119286*fL[28])*amax+21.21320343559643*alphaR[0]*fAvg[15]+21.21320343559643*alphaR[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[43]+fL[43])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+21.21320343559643*alphaR[1]*fAvg[19]+21.21320343559643*alphaR[0]*fAvg[16]); 
  Ghat[17] = 0.025*((24.49489742783179*(fR[44]+fL[44])-14.14213562373095*fR[37]+14.14213562373095*fL[37])*amax+7.071067811865476*alphaR[0]*fAvg[17]+6.324555320336761*alphaR[1]*fAvg[10]); 
  Ghat[18] = 0.008333333333333333*((73.48469228349536*(fR[46]+fL[46])-42.42640687119286*fR[39]+42.42640687119286*fL[39])*amax+21.21320343559643*alphaR[0]*fAvg[18]+21.21320343559643*alphaR[1]*fAvg[14]); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[42]+42.42640687119286*fL[42])*amax+21.21320343559643*alphaR[0]*fAvg[19]+21.21320343559643*alphaR[1]*fAvg[16]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = -1.224744871391589*Ghat[0]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = 0.7071067811865475*Ghat[4]; 
  incr[7] = -1.224744871391589*Ghat[2]; 
  incr[8] = 0.7071067811865475*Ghat[5]; 
  incr[9] = -1.224744871391589*Ghat[3]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 0.7071067811865475*Ghat[7]; 
  incr[12] = 1.58113883008419*Ghat[0]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = 0.7071067811865475*Ghat[10]; 
  incr[18] = -1.224744871391589*Ghat[6]; 
  incr[19] = -1.224744871391589*Ghat[7]; 
  incr[20] = 1.58113883008419*Ghat[1]; 
  incr[21] = 0.7071067811865475*Ghat[11]; 
  incr[22] = 1.58113883008419*Ghat[2]; 
  incr[23] = 0.7071067811865475*Ghat[12]; 
  incr[24] = -1.224744871391589*Ghat[8]; 
  incr[25] = 0.7071067811865475*Ghat[13]; 
  incr[26] = 1.58113883008419*Ghat[3]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = 0.7071067811865475*Ghat[15]; 
  incr[29] = -1.224744871391589*Ghat[9]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = -1.224744871391589*Ghat[11]; 
  incr[33] = 1.58113883008419*Ghat[4]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = -1.224744871391589*Ghat[13]; 
  incr[36] = 1.58113883008419*Ghat[5]; 
  incr[37] = 0.7071067811865475*Ghat[17]; 
  incr[38] = 1.58113883008419*Ghat[6]; 
  incr[39] = 0.7071067811865475*Ghat[18]; 
  incr[40] = -1.224744871391589*Ghat[14]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = 0.7071067811865475*Ghat[19]; 
  incr[43] = -1.224744871391589*Ghat[16]; 
  incr[44] = -1.224744871391589*Ghat[17]; 
  incr[45] = 1.58113883008419*Ghat[10]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

  outR[0] += incr[0]*rdy2R; 
  outR[1] += incr[1]*rdy2R; 
  outR[2] += incr[2]*rdy2R; 
  outR[3] += incr[3]*rdy2R; 
  outR[4] += incr[4]*rdy2R; 
  outR[5] += incr[5]*rdy2R; 
  outR[6] += incr[6]*rdy2R; 
  outR[7] += incr[7]*rdy2R; 
  outR[8] += incr[8]*rdy2R; 
  outR[9] += incr[9]*rdy2R; 
  outR[10] += incr[10]*rdy2R; 
  outR[11] += incr[11]*rdy2R; 
  outR[12] += incr[12]*rdy2R; 
  outR[13] += incr[13]*rdy2R; 
  outR[14] += incr[14]*rdy2R; 
  outR[15] += incr[15]*rdy2R; 
  outR[16] += incr[16]*rdy2R; 
  outR[17] += incr[17]*rdy2R; 
  outR[18] += incr[18]*rdy2R; 
  outR[19] += incr[19]*rdy2R; 
  outR[20] += incr[20]*rdy2R; 
  outR[21] += incr[21]*rdy2R; 
  outR[22] += incr[22]*rdy2R; 
  outR[23] += incr[23]*rdy2R; 
  outR[24] += incr[24]*rdy2R; 
  outR[25] += incr[25]*rdy2R; 
  outR[26] += incr[26]*rdy2R; 
  outR[27] += incr[27]*rdy2R; 
  outR[28] += incr[28]*rdy2R; 
  outR[29] += incr[29]*rdy2R; 
  outR[30] += incr[30]*rdy2R; 
  outR[31] += incr[31]*rdy2R; 
  outR[32] += incr[32]*rdy2R; 
  outR[33] += incr[33]*rdy2R; 
  outR[34] += incr[34]*rdy2R; 
  outR[35] += incr[35]*rdy2R; 
  outR[36] += incr[36]*rdy2R; 
  outR[37] += incr[37]*rdy2R; 
  outR[38] += incr[38]*rdy2R; 
  outR[39] += incr[39]*rdy2R; 
  outR[40] += incr[40]*rdy2R; 
  outR[41] += incr[41]*rdy2R; 
  outR[42] += incr[42]*rdy2R; 
  outR[43] += incr[43]*rdy2R; 
  outR[44] += incr[44]*rdy2R; 
  outR[45] += incr[45]*rdy2R; 
  outR[46] += incr[46]*rdy2R; 
  outR[47] += incr[47]*rdy2R; 

  outL[0] += -1.0*incr[0]*rdy2L; 
  outL[1] += -1.0*incr[1]*rdy2L; 
  outL[2] += incr[2]*rdy2L; 
  outL[3] += -1.0*incr[3]*rdy2L; 
  outL[4] += -1.0*incr[4]*rdy2L; 
  outL[5] += incr[5]*rdy2L; 
  outL[6] += -1.0*incr[6]*rdy2L; 
  outL[7] += incr[7]*rdy2L; 
  outL[8] += -1.0*incr[8]*rdy2L; 
  outL[9] += incr[9]*rdy2L; 
  outL[10] += -1.0*incr[10]*rdy2L; 
  outL[11] += -1.0*incr[11]*rdy2L; 
  outL[12] += -1.0*incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += -1.0*incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 
  outL[16] += incr[16]*rdy2L; 
  outL[17] += -1.0*incr[17]*rdy2L; 
  outL[18] += incr[18]*rdy2L; 
  outL[19] += incr[19]*rdy2L; 
  outL[20] += -1.0*incr[20]*rdy2L; 
  outL[21] += -1.0*incr[21]*rdy2L; 
  outL[22] += -1.0*incr[22]*rdy2L; 
  outL[23] += -1.0*incr[23]*rdy2L; 
  outL[24] += incr[24]*rdy2L; 
  outL[25] += -1.0*incr[25]*rdy2L; 
  outL[26] += -1.0*incr[26]*rdy2L; 
  outL[27] += -1.0*incr[27]*rdy2L; 
  outL[28] += -1.0*incr[28]*rdy2L; 
  outL[29] += incr[29]*rdy2L; 
  outL[30] += -1.0*incr[30]*rdy2L; 
  outL[31] += incr[31]*rdy2L; 
  outL[32] += incr[32]*rdy2L; 
  outL[33] += -1.0*incr[33]*rdy2L; 
  outL[34] += incr[34]*rdy2L; 
  outL[35] += incr[35]*rdy2L; 
  outL[36] += -1.0*incr[36]*rdy2L; 
  outL[37] += -1.0*incr[37]*rdy2L; 
  outL[38] += -1.0*incr[38]*rdy2L; 
  outL[39] += -1.0*incr[39]*rdy2L; 
  outL[40] += incr[40]*rdy2L; 
  outL[41] += incr[41]*rdy2L; 
  outL[42] += -1.0*incr[42]*rdy2L; 
  outL[43] += incr[43]*rdy2L; 
  outL[44] += incr[44]*rdy2L; 
  outL[45] += -1.0*incr[45]*rdy2L; 
  outL[46] += incr[46]*rdy2L; 
  outL[47] += incr[47]*rdy2L; 
return std::abs(alphaSurfAvgR); 
} 
double GyrokineticGenGeoSurf2x2vSer_vpar_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[11] = 2.0*phi[4]*q_; 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 

  double BstarXdBmagR[48]; 

  double BstarYdBmagR[48]; 

  double alphaR[20]; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = 0.0; 

  double incr[48]; 
  // alpha == 0, so nothing to do.
  return std::abs(alphaSurfAvgR);
}
double GyrokineticGenGeoSurf2x2vSer_x_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[8] = (1.154700538379252*bmag[1])/rdmu2R; 
  hamilR[11] = 2.0*(bmag[4]*wmuR+phi[4]*q_); 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 
  hamilR[25] = (1.154700538379251*bmag[4])/rdmu2R; 

  double BstarXdBmagR[48]; 

  double alphaR[20]; 
  alphaR[0] = -(0.125*(((27.38612787525831*b_z[4]-21.21320343559643*b_z[1]+12.24744871391589*b_z[0])*jacobTotInv[4]+(12.24744871391589*jacobTotInv[0]-21.21320343559643*jacobTotInv[1])*b_z[4]+(16.43167672515499*b_z[1]-9.48683298050514*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(5.477225575051662*b_z[0]-9.48683298050514*b_z[1]))*hamilR[19]+(((-21.21320343559643*b_z[4])+16.43167672515498*b_z[1]-9.48683298050514*b_z[0])*jacobTotInv[4]+(16.43167672515498*jacobTotInv[1]-9.48683298050514*jacobTotInv[0])*b_z[4]+(7.348469228349534*b_z[0]-12.72792206135786*b_z[1])*jacobTotInv[1]+jacobTotInv[0]*(7.348469228349534*b_z[1]-4.242640687119286*b_z[0]))*hamilR[5]+hamilR[2]*((12.24744871391589*b_z[4]-9.48683298050514*b_z[1]+5.477225575051662*b_z[0])*jacobTotInv[4]+(5.477225575051662*jacobTotInv[0]-9.48683298050514*jacobTotInv[1])*b_z[4]+(7.348469228349534*b_z[1]-4.242640687119286*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(2.449489742783178*b_z[0]-4.242640687119286*b_z[1])))*rdy2R)/q_; 
  alphaR[1] = (0.125*(((47.43416490252569*b_z[4]-36.74234614174768*b_z[1]+21.21320343559643*b_z[0])*jacobTotInv[4]+(21.21320343559643*jacobTotInv[0]-36.74234614174768*jacobTotInv[1])*b_z[4]+(28.46049894151541*b_z[1]-16.43167672515499*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(9.48683298050514*b_z[0]-16.43167672515499*b_z[1]))*hamilR[20]+(((-27.38612787525831*b_z[4])+21.21320343559643*b_z[1]-12.24744871391589*b_z[0])*jacobTotInv[4]+(21.21320343559643*jacobTotInv[1]-12.24744871391589*jacobTotInv[0])*b_z[4]+(9.48683298050514*b_z[0]-16.43167672515498*b_z[1])*jacobTotInv[1]+jacobTotInv[0]*(9.48683298050514*b_z[1]-5.477225575051662*b_z[0]))*hamilR[12])*rdy2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(((19.36491673103708*b_z[4]-15.0*b_z[1]+8.660254037844387*b_z[0])*jacobTotInv[4]+(8.660254037844387*jacobTotInv[0]-15.0*jacobTotInv[1])*b_z[4]+(11.61895003862225*b_z[1]-6.708203932499369*b_z[0])*jacobTotInv[1]-6.708203932499369*jacobTotInv[0]*b_z[1]+3.872983346207417*b_z[0]*jacobTotInv[0])*hamilR[19]+(((-15.0*b_z[4])+11.61895003862225*b_z[1]-6.708203932499369*b_z[0])*jacobTotInv[4]+(11.61895003862225*jacobTotInv[1]-6.708203932499369*jacobTotInv[0])*b_z[4]+(5.196152422706631*b_z[0]-9.0*b_z[1])*jacobTotInv[1]+5.196152422706631*jacobTotInv[0]*b_z[1]-3.0*b_z[0]*jacobTotInv[0])*hamilR[5]+(8.660254037844386*hamilR[2]*b_z[4]+(3.872983346207417*b_z[0]-6.708203932499369*b_z[1])*hamilR[2])*jacobTotInv[4]+(3.872983346207417*jacobTotInv[0]-6.708203932499369*jacobTotInv[1])*hamilR[2]*b_z[4]+((5.196152422706631*b_z[1]-3.0*b_z[0])*jacobTotInv[1]-3.0*jacobTotInv[0]*b_z[1]+1.732050807568877*b_z[0]*jacobTotInv[0])*hamilR[2])*rdy2R)/q_; 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[11]+fL[11])+1.732050807568877*(fL[1]-1.0*fR[1])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[19]+fL[19])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[2]+fL[2])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[21]+fL[21])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[25]+fL[25])+3.0*(fL[8]-1.0*fR[8]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[32]+fL[32])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[7]+fL[7]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[35]+fL[35])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[9]+fL[9]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[37]+fL[37])+1.732050807568877*(fL[17]-1.0*fR[17])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[20]-1.0*(8.660254037844387*fL[20]+5.0*(fR[12]+fL[12]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[23]-1.0*(8.660254037844387*fL[23]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[28]-1.0*(8.660254037844387*fL[28]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[44]+fL[44])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[18]+fL[18])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[33]-1.0*(8.660254037844387*fL[33]+5.0*(fR[22]+fL[22]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[24]+fL[24]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[36]-1.0*(8.660254037844387*fL[36]+5.0*(fR[26]+fL[26]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[39]-1.0*(8.660254037844387*fL[39]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[29]+fL[29]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[42]-1.0*(8.660254037844387*fL[42]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[45]-1.0*(8.660254037844387*fL[45]+5.0*(fR[38]+fL[38]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[40]+fL[40]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[43]+fL[43]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[11]-6.324555320336761*fL[11]-4.898979485566357*(fR[1]+fL[1])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[19]-94.8683298050514*fL[19]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[2]-42.42640687119286*fL[2])*amax-18.97366596101028*alphaR[1]*fAvg[7]-21.21320343559643*(alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.04166666666666666*((18.97366596101028*fR[21]-18.97366596101028*fL[21]-14.69693845669907*(fR[6]+fL[6])+8.485281374238571*fR[3]-8.485281374238571*fL[3])*amax-4.242640687119286*(alphaR[1]*fAvg[4]+alphaR[0]*fAvg[2])); 
  Ghat[3] = -0.04166666666666666*((18.97366596101028*fR[25]-18.97366596101028*fL[25]-14.69693845669907*(fR[8]+fL[8])+8.485281374238571*fR[4]-8.485281374238571*fL[4])*amax-4.242640687119286*(alphaR[1]*fAvg[5]+alphaR[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[32]-94.86832980505142*fL[32]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[7]-42.42640687119286*fL[7])*amax-18.97366596101028*alphaR[1]*fAvg[11]-21.21320343559643*(alphaR[0]*fAvg[4]+alphaR[1]*fAvg[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[35]-94.86832980505142*fL[35]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[9]-42.42640687119286*fL[9])*amax-18.97366596101028*alphaR[1]*fAvg[13]-21.21320343559643*(alphaR[0]*fAvg[5]+alphaR[1]*fAvg[3])); 
  Ghat[6] = -0.125*((6.324555320336761*fR[37]-6.324555320336761*fL[37]-4.898979485566357*(fR[17]+fL[17])+2.828427124746191*fR[10]-2.828427124746191*fL[10])*amax-1.414213562373095*(alphaR[1]*fAvg[10]+alphaR[0]*fAvg[6])); 
  Ghat[7] = 0.025*((24.49489742783179*(fR[20]+fL[20])-14.14213562373095*fR[12]+14.14213562373095*fL[12])*amax+7.071067811865476*alphaR[0]*fAvg[7]+6.324555320336761*alphaR[1]*fAvg[1]); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[23]+fL[23])-42.42640687119286*fR[13]+42.42640687119286*fL[13])*amax+21.21320343559643*alphaR[1]*fAvg[12]+21.21320343559643*alphaR[0]*fAvg[8]); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[28]+fL[28])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[1]*fAvg[15]+21.21320343559643*alphaR[0]*fAvg[9]); 
  Ghat[10] = -0.008333333333333333*((94.8683298050514*fR[44]-94.8683298050514*fL[44]-73.48469228349535*(fR[31]+fL[31])+42.42640687119286*fR[18]-42.42640687119286*fL[18])*amax-18.97366596101028*alphaR[1]*fAvg[17]-21.21320343559643*(alphaR[0]*fAvg[10]+alphaR[1]*fAvg[6])); 
  Ghat[11] = 0.008333333333333333*((73.48469228349536*(fR[33]+fL[33])-42.42640687119286*fR[22]+42.42640687119286*fL[22])*amax+21.21320343559643*alphaR[0]*fAvg[11]+18.97366596101028*alphaR[1]*fAvg[4]); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[34]+fL[34])-42.42640687119286*fR[24]+42.42640687119286*fL[24])*amax+21.21320343559643*alphaR[0]*fAvg[12]+21.21320343559643*alphaR[1]*fAvg[8]); 
  Ghat[13] = 0.008333333333333333*((73.48469228349536*(fR[36]+fL[36])-42.42640687119286*fR[26]+42.42640687119286*fL[26])*amax+21.21320343559643*alphaR[0]*fAvg[13]+18.97366596101028*alphaR[1]*fAvg[5]); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[39]+fL[39])-42.42640687119286*fR[27]+42.42640687119286*fL[27])*amax+21.21320343559643*alphaR[1]*fAvg[18]+21.21320343559643*alphaR[0]*fAvg[14]); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[29]+42.42640687119286*fL[29])*amax+21.21320343559643*alphaR[0]*fAvg[15]+21.21320343559643*alphaR[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[42]+fL[42])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+21.21320343559643*alphaR[1]*fAvg[19]+21.21320343559643*alphaR[0]*fAvg[16]); 
  Ghat[17] = 0.025*((24.49489742783179*(fR[45]+fL[45])-14.14213562373095*fR[38]+14.14213562373095*fL[38])*amax+7.071067811865476*alphaR[0]*fAvg[17]+6.324555320336761*alphaR[1]*fAvg[10]); 
  Ghat[18] = 0.008333333333333333*((73.48469228349536*(fR[46]+fL[46])-42.42640687119286*fR[40]+42.42640687119286*fL[40])*amax+21.21320343559643*alphaR[0]*fAvg[18]+21.21320343559643*alphaR[1]*fAvg[14]); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[43]+42.42640687119286*fL[43])*amax+21.21320343559643*alphaR[0]*fAvg[19]+21.21320343559643*alphaR[1]*fAvg[16]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = -1.224744871391589*Ghat[0]; 
  incr[2] = 0.7071067811865475*Ghat[1]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = -1.224744871391589*Ghat[2]; 
  incr[7] = 0.7071067811865475*Ghat[4]; 
  incr[8] = -1.224744871391589*Ghat[3]; 
  incr[9] = 0.7071067811865475*Ghat[5]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 1.58113883008419*Ghat[0]; 
  incr[12] = 0.7071067811865475*Ghat[7]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = -1.224744871391589*Ghat[6]; 
  incr[18] = 0.7071067811865475*Ghat[10]; 
  incr[19] = 1.58113883008419*Ghat[1]; 
  incr[20] = -1.224744871391589*Ghat[7]; 
  incr[21] = 1.58113883008419*Ghat[2]; 
  incr[22] = 0.7071067811865475*Ghat[11]; 
  incr[23] = -1.224744871391589*Ghat[8]; 
  incr[24] = 0.7071067811865475*Ghat[12]; 
  incr[25] = 1.58113883008419*Ghat[3]; 
  incr[26] = 0.7071067811865475*Ghat[13]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = -1.224744871391589*Ghat[9]; 
  incr[29] = 0.7071067811865475*Ghat[15]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = 1.58113883008419*Ghat[4]; 
  incr[33] = -1.224744871391589*Ghat[11]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = 1.58113883008419*Ghat[5]; 
  incr[36] = -1.224744871391589*Ghat[13]; 
  incr[37] = 1.58113883008419*Ghat[6]; 
  incr[38] = 0.7071067811865475*Ghat[17]; 
  incr[39] = -1.224744871391589*Ghat[14]; 
  incr[40] = 0.7071067811865475*Ghat[18]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = -1.224744871391589*Ghat[16]; 
  incr[43] = 0.7071067811865475*Ghat[19]; 
  incr[44] = 1.58113883008419*Ghat[10]; 
  incr[45] = -1.224744871391589*Ghat[17]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 
  outR[8] += incr[8]*rdx2R; 
  outR[9] += incr[9]*rdx2R; 
  outR[10] += incr[10]*rdx2R; 
  outR[11] += incr[11]*rdx2R; 
  outR[12] += incr[12]*rdx2R; 
  outR[13] += incr[13]*rdx2R; 
  outR[14] += incr[14]*rdx2R; 
  outR[15] += incr[15]*rdx2R; 
  outR[16] += incr[16]*rdx2R; 
  outR[17] += incr[17]*rdx2R; 
  outR[18] += incr[18]*rdx2R; 
  outR[19] += incr[19]*rdx2R; 
  outR[20] += incr[20]*rdx2R; 
  outR[21] += incr[21]*rdx2R; 
  outR[22] += incr[22]*rdx2R; 
  outR[23] += incr[23]*rdx2R; 
  outR[24] += incr[24]*rdx2R; 
  outR[25] += incr[25]*rdx2R; 
  outR[26] += incr[26]*rdx2R; 
  outR[27] += incr[27]*rdx2R; 
  outR[28] += incr[28]*rdx2R; 
  outR[29] += incr[29]*rdx2R; 
  outR[30] += incr[30]*rdx2R; 
  outR[31] += incr[31]*rdx2R; 
  outR[32] += incr[32]*rdx2R; 
  outR[33] += incr[33]*rdx2R; 
  outR[34] += incr[34]*rdx2R; 
  outR[35] += incr[35]*rdx2R; 
  outR[36] += incr[36]*rdx2R; 
  outR[37] += incr[37]*rdx2R; 
  outR[38] += incr[38]*rdx2R; 
  outR[39] += incr[39]*rdx2R; 
  outR[40] += incr[40]*rdx2R; 
  outR[41] += incr[41]*rdx2R; 
  outR[42] += incr[42]*rdx2R; 
  outR[43] += incr[43]*rdx2R; 
  outR[44] += incr[44]*rdx2R; 
  outR[45] += incr[45]*rdx2R; 
  outR[46] += incr[46]*rdx2R; 
  outR[47] += incr[47]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += -1.0*incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += incr[6]*rdx2L; 
  outL[7] += -1.0*incr[7]*rdx2L; 
  outL[8] += incr[8]*rdx2L; 
  outL[9] += -1.0*incr[9]*rdx2L; 
  outL[10] += -1.0*incr[10]*rdx2L; 
  outL[11] += -1.0*incr[11]*rdx2L; 
  outL[12] += -1.0*incr[12]*rdx2L; 
  outL[13] += -1.0*incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 
  outL[16] += incr[16]*rdx2L; 
  outL[17] += incr[17]*rdx2L; 
  outL[18] += -1.0*incr[18]*rdx2L; 
  outL[19] += -1.0*incr[19]*rdx2L; 
  outL[20] += incr[20]*rdx2L; 
  outL[21] += -1.0*incr[21]*rdx2L; 
  outL[22] += -1.0*incr[22]*rdx2L; 
  outL[23] += incr[23]*rdx2L; 
  outL[24] += -1.0*incr[24]*rdx2L; 
  outL[25] += -1.0*incr[25]*rdx2L; 
  outL[26] += -1.0*incr[26]*rdx2L; 
  outL[27] += -1.0*incr[27]*rdx2L; 
  outL[28] += incr[28]*rdx2L; 
  outL[29] += -1.0*incr[29]*rdx2L; 
  outL[30] += -1.0*incr[30]*rdx2L; 
  outL[31] += incr[31]*rdx2L; 
  outL[32] += -1.0*incr[32]*rdx2L; 
  outL[33] += incr[33]*rdx2L; 
  outL[34] += incr[34]*rdx2L; 
  outL[35] += -1.0*incr[35]*rdx2L; 
  outL[36] += incr[36]*rdx2L; 
  outL[37] += -1.0*incr[37]*rdx2L; 
  outL[38] += -1.0*incr[38]*rdx2L; 
  outL[39] += incr[39]*rdx2L; 
  outL[40] += -1.0*incr[40]*rdx2L; 
  outL[41] += incr[41]*rdx2L; 
  outL[42] += incr[42]*rdx2L; 
  outL[43] += -1.0*incr[43]*rdx2L; 
  outL[44] += -1.0*incr[44]*rdx2L; 
  outL[45] += incr[45]*rdx2L; 
  outL[46] += incr[46]*rdx2L; 
  outL[47] += incr[47]*rdx2L; 
return std::abs(alphaSurfAvgR); 
} 
double GyrokineticGenGeoSurf2x2vSer_y_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[8] = (1.154700538379252*bmag[1])/rdmu2R; 
  hamilR[11] = 2.0*(bmag[4]*wmuR+phi[4]*q_); 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 
  hamilR[25] = (1.154700538379251*bmag[4])/rdmu2R; 

  double BstarYdBmagR[48]; 
  BstarYdBmagR[0] = -(1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[1] = -(1.732050807568877*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[3] = -(1.0*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[6] = -(1.0*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[11] = -(1.732050807568877*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[21] = -(1.0*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[20]; 
  alphaR[0] = (0.03535533905932736*((19.36491673103708*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamilR[20]+((-30.0*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-33.54101966249684*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[19]+(17.32050807568877*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+19.36491673103709*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[11]+(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*(8.660254037844386*hamilR[1]-15.0*hamilR[5]))*m_*rdx2R+(19.36491673103709*BstarYdBmagR[3]*hamilR[13]+8.660254037844386*BstarYdBmagR[0]*hamilR[3])*q_*rdvpar2R))/(m_*q_); 
  alphaR[1] = (0.005050762722761052*(((121.2435565298214*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+135.5544171172596*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[20]+(((-368.9512162874653*b_z[4])-210.0*b_z[0])*jacobTotInv[4]-210.0*jacobTotInv[0]*b_z[4]-422.6168477474601*b_z[1]*jacobTotInv[1]-234.7871376374779*b_z[0]*jacobTotInv[0])*hamilR[19]+((213.014084041408*b_z[4]+121.2435565298214*b_z[0])*jacobTotInv[4]+121.2435565298214*jacobTotInv[0]*b_z[4]+243.9979508110672*b_z[1]*jacobTotInv[1]+135.5544171172596*b_z[0]*jacobTotInv[0])*hamilR[11]+((-93.91485505499116*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-105.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[5]+hamilR[1]*(54.22176684690384*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+60.6217782649107*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])))*m_*rdx2R+(135.5544171172596*BstarYdBmagR[6]*hamilR[13]+60.6217782649107*BstarYdBmagR[1]*hamilR[3])*q_*rdvpar2R))/(m_*q_); 
  alphaR[2] = (0.1767766952966368*(3.872983346207417*BstarYdBmagR[0]*hamilR[13]+1.732050807568877*BstarYdBmagR[3]*hamilR[3])*rdvpar2R)/m_; 
  alphaR[3] = (0.03535533905932736*((17.32050807568877*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+19.36491673103708*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[25]+8.660254037844386*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamilR[8])*rdx2R)/q_; 
  alphaR[4] = (0.1767766952966368*(3.872983346207417*BstarYdBmagR[1]*hamilR[13]+1.732050807568877*hamilR[3]*BstarYdBmagR[6])*rdvpar2R)/m_; 
  alphaR[5] = (0.005050762722761052*(((213.0140840414079*b_z[4]+121.2435565298214*b_z[0])*jacobTotInv[4]+121.2435565298214*jacobTotInv[0]*b_z[4]+243.9979508110673*b_z[1]*jacobTotInv[1]+135.5544171172596*b_z[0]*jacobTotInv[0])*hamilR[25]+(54.22176684690384*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+60.6217782649107*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[8])*rdx2R)/q_; 
  alphaR[7] = (0.005050762722761052*((((86.60254037844389*b_z[4]+135.5544171172596*b_z[0])*jacobTotInv[4]+135.5544171172596*jacobTotInv[0]*b_z[4]+121.2435565298214*b_z[1]*jacobTotInv[1])*hamilR[20]+((-368.9512162874653*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-210.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[19]+(213.014084041408*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+121.2435565298214*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[11]+(((-67.0820393249937*b_z[4])-105.0*b_z[0])*jacobTotInv[4]-105.0*jacobTotInv[0]*b_z[4]-93.91485505499116*b_z[1]*jacobTotInv[1])*hamilR[5]+hamilR[1]*((38.72983346207418*b_z[4]+60.6217782649107*b_z[0])*jacobTotInv[4]+60.6217782649107*jacobTotInv[0]*b_z[4]+54.22176684690384*b_z[1]*jacobTotInv[1]))*m_*rdx2R+(135.5544171172596*hamilR[13]*BstarYdBmagR[21]+60.6217782649107*hamilR[3]*BstarYdBmagR[11])*q_*rdvpar2R))/(m_*q_); 
  alphaR[8] = (0.6123724356957944*BstarYdBmagR[3]*hamilR[13]*rdvpar2R)/m_; 
  alphaR[11] = (0.04564354645876382*(6.708203932499369*hamilR[3]*BstarYdBmagR[21]+15.0*BstarYdBmagR[11]*hamilR[13])*rdvpar2R)/m_; 
  alphaR[12] = (0.6123724356957944*BstarYdBmagR[6]*hamilR[13]*rdvpar2R)/m_; 
  alphaR[13] = (0.005050762722761052*((213.014084041408*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+121.2435565298214*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[25]+((38.72983346207417*b_z[4]+60.62177826491071*b_z[0])*jacobTotInv[4]+60.62177826491071*jacobTotInv[0]*b_z[4]+54.22176684690384*b_z[1]*jacobTotInv[1])*hamilR[8])*rdx2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.00625*(((19.36491673103708*b_z[4]*jacobTotInv[4]+19.36491673103708*b_z[1]*jacobTotInv[1]+19.36491673103708*b_z[0]*jacobTotInv[0])*hamilR[20]+((-30.0*b_z[1]*jacobTotInv[4])-30.0*jacobTotInv[1]*b_z[4]-33.54101966249684*b_z[0]*jacobTotInv[1]-33.54101966249684*jacobTotInv[0]*b_z[1])*hamilR[19]+(17.32050807568877*b_z[1]*jacobTotInv[4]+17.32050807568877*jacobTotInv[1]*b_z[4]+19.36491673103709*b_z[0]*jacobTotInv[1]+19.36491673103709*jacobTotInv[0]*b_z[1])*hamilR[11]+((-15.0*b_z[4]*jacobTotInv[4])-15.0*b_z[1]*jacobTotInv[1]-15.0*b_z[0]*jacobTotInv[0])*hamilR[5]+8.660254037844386*hamilR[1]*b_z[4]*jacobTotInv[4]+8.660254037844386*b_z[1]*hamilR[1]*jacobTotInv[1]+8.660254037844386*b_z[0]*jacobTotInv[0]*hamilR[1])*m_*rdx2R+(19.36491673103709*BstarYdBmagR[3]*hamilR[13]+8.660254037844386*BstarYdBmagR[0]*hamilR[3])*q_*rdvpar2R))/(m_*q_); 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[12]+fL[12])+1.732050807568877*(fL[2]-1.0*fR[2])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[20]+fL[20])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[22]+fL[22])+3.0*(fL[7]-1.0*fR[7]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[26]+fL[26])+3.0*(fL[9]-1.0*fR[9]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[33]+fL[33])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[6]+fL[6]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[36]+fL[36])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[8]+fL[8]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[38]+fL[38])+1.732050807568877*(fL[18]-1.0*fR[18])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[19]-1.0*(8.660254037844387*fL[19]+5.0*(fR[11]+fL[11]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[24]-1.0*(8.660254037844387*fL[24]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[29]-1.0*(8.660254037844387*fL[29]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[45]+fL[45])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[17]+fL[17])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[32]-1.0*(8.660254037844387*fL[32]+5.0*(fR[21]+fL[21]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[23]+fL[23]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[35]-1.0*(8.660254037844387*fL[35]+5.0*(fR[25]+fL[25]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[40]-1.0*(8.660254037844387*fL[40]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[28]+fL[28]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[43]-1.0*(8.660254037844387*fL[43]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[44]-1.0*(8.660254037844387*fL[44]+5.0*(fR[37]+fL[37]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[39]+fL[39]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[42]+fL[42]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[12]-6.324555320336761*fL[12]-4.898979485566357*(fR[2]+fL[2])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[13]*fAvg[13]+alphaR[12]*fAvg[12]+alphaR[11]*fAvg[11]+alphaR[8]*fAvg[8]+alphaR[7]*fAvg[7]+alphaR[5]*fAvg[5]+alphaR[4]*fAvg[4]+alphaR[3]*fAvg[3]+alphaR[2]*fAvg[2]+alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[20]-94.8683298050514*fL[20]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[1]-42.42640687119286*fL[1])*amax-18.97366596101028*(alphaR[5]*fAvg[13]+fAvg[5]*alphaR[13])-21.21320343559643*(alphaR[8]*fAvg[12]+fAvg[8]*alphaR[12])-18.97366596101028*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11]+alphaR[1]*fAvg[7]+fAvg[1]*alphaR[7])-21.21320343559643*(alphaR[3]*fAvg[5]+fAvg[3]*alphaR[5]+alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4]+alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[22]-94.8683298050514*fL[22]-73.48469228349535*(fR[7]+fL[7])+42.42640687119286*fR[3]-42.42640687119286*fL[3])*amax-21.21320343559643*alphaR[13]*fAvg[17]-18.97366596101028*(alphaR[4]*fAvg[12]+fAvg[4]*alphaR[12])-21.21320343559643*(alphaR[7]*fAvg[11]+fAvg[7]*alphaR[11])-21.21320343559643*alphaR[5]*fAvg[10]-18.97366596101028*(alphaR[2]*fAvg[8]+fAvg[2]*alphaR[8])-21.21320343559643*(alphaR[3]*fAvg[6]+alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4]+alphaR[0]*fAvg[2]+fAvg[0]*alphaR[2])); 
  Ghat[3] = -0.008333333333333333*((94.8683298050514*fR[26]-94.8683298050514*fL[26]-73.48469228349535*(fR[9]+fL[9])+42.42640687119286*fR[4]-42.42640687119286*fL[4])*amax-21.21320343559643*(alphaR[12]*fAvg[18]+alphaR[11]*fAvg[17])-18.97366596101028*alphaR[5]*fAvg[15]-21.21320343559643*(alphaR[8]*fAvg[14]+alphaR[7]*fAvg[13]+fAvg[7]*alphaR[13])-21.21320343559643*alphaR[4]*fAvg[10]-18.97366596101028*alphaR[3]*fAvg[9]-21.21320343559643*(alphaR[2]*fAvg[6]+alphaR[1]*fAvg[5]+fAvg[1]*alphaR[5]+alphaR[0]*fAvg[3]+fAvg[0]*alphaR[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[33]-94.86832980505142*fL[33]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[6]-42.42640687119286*fL[6])*amax-18.97366596101028*(alphaR[5]*fAvg[17]+fAvg[10]*alphaR[13])+((-16.97056274847715*alphaR[11])-18.97366596101028*alphaR[2])*fAvg[12]+((-16.97056274847715*fAvg[11])-18.97366596101028*fAvg[2])*alphaR[12]-18.97366596101028*(alphaR[1]*fAvg[11]+fAvg[1]*alphaR[11])-21.21320343559643*alphaR[3]*fAvg[10]-18.97366596101028*(alphaR[4]*fAvg[8]+fAvg[4]*alphaR[8]+alphaR[4]*fAvg[7]+fAvg[4]*alphaR[7])-21.21320343559643*(alphaR[5]*fAvg[6]+alphaR[0]*fAvg[4]+fAvg[0]*alphaR[4]+alphaR[1]*fAvg[2]+fAvg[1]*alphaR[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[36]-94.86832980505142*fL[36]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[8]-42.42640687119286*fL[8])*amax-21.21320343559643*alphaR[8]*fAvg[18]-18.97366596101028*alphaR[4]*fAvg[17]+((-16.97056274847715*alphaR[13])-18.97366596101028*alphaR[3])*fAvg[15]-21.21320343559643*alphaR[12]*fAvg[14]-18.97366596101028*(alphaR[1]*fAvg[13]+fAvg[1]*alphaR[13])+fAvg[10]*((-18.97366596101028*alphaR[11])-21.21320343559643*alphaR[2])-18.97366596101028*(alphaR[5]*(fAvg[9]+fAvg[7])+fAvg[5]*alphaR[7])-21.21320343559643*(alphaR[4]*fAvg[6]+alphaR[0]*fAvg[5]+fAvg[0]*alphaR[5]+alphaR[1]*fAvg[3]+fAvg[1]*alphaR[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[38]-94.86832980505142*fL[38]-73.48469228349535*(fR[18]+fL[18])+42.42640687119286*fR[10]-42.42640687119286*fL[10])*amax-18.97366596101028*(alphaR[5]*fAvg[19]+alphaR[4]*fAvg[18])-21.21320343559643*alphaR[7]*fAvg[17]-18.97366596101028*(alphaR[3]*fAvg[16]+alphaR[2]*fAvg[14])-21.21320343559643*(alphaR[11]*fAvg[13]+fAvg[11]*alphaR[13])+fAvg[10]*((-18.97366596101028*alphaR[12])-21.21320343559643*alphaR[1])-18.97366596101028*fAvg[6]*alphaR[8]-21.21320343559643*(alphaR[0]*fAvg[6]+alphaR[4]*fAvg[5]+fAvg[4]*alphaR[5]+alphaR[2]*fAvg[3]+fAvg[2]*alphaR[3])); 
  Ghat[7] = 0.001190476190476191*((514.3928459844675*(fR[19]+fL[19])-296.98484809835*fR[11]+296.98484809835*fL[11])*amax+(94.86832980505142*alphaR[13]+148.492424049175*alphaR[3])*fAvg[13]+148.492424049175*fAvg[3]*alphaR[13]+132.815661727072*alphaR[12]*fAvg[12]+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[11]+148.492424049175*fAvg[2]*alphaR[11]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[7]+148.492424049175*fAvg[0]*alphaR[7]+132.815661727072*(alphaR[5]*fAvg[5]+alphaR[4]*fAvg[4]+alphaR[1]*fAvg[1])); 
  Ghat[8] = 0.001190476190476191*((514.3928459844675*(fR[24]+fL[24])-296.98484809835*fR[13]+296.98484809835*fL[13])*amax+148.492424049175*(alphaR[5]*fAvg[18]+alphaR[3]*fAvg[14])+(94.86832980505142*alphaR[12]+148.492424049175*alphaR[1])*fAvg[12]+148.492424049175*fAvg[1]*alphaR[12]+132.815661727072*alphaR[11]*fAvg[11]+(94.86832980505142*alphaR[8]+148.492424049175*alphaR[0])*fAvg[8]+148.492424049175*fAvg[0]*alphaR[8]+132.815661727072*(alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[29]+fL[29])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[4]*fAvg[19]+21.21320343559643*(alphaR[2]*fAvg[16]+alphaR[1]*fAvg[15])+18.97366596101028*alphaR[13]*fAvg[13]+21.21320343559643*alphaR[0]*fAvg[9]+18.97366596101028*(alphaR[5]*fAvg[5]+alphaR[3]*fAvg[3])); 
  Ghat[10] = -0.001666666666666667*((474.341649025257*fR[45]-474.341649025257*fL[45]-367.4234614174767*(fR[31]+fL[31])+212.1320343559643*fR[17]-212.1320343559643*fL[17])*amax+((-84.85281374238573*alphaR[13])-94.86832980505142*alphaR[3])*fAvg[19]+((-84.85281374238573*alphaR[11])-94.86832980505142*alphaR[2])*fAvg[18]+((-84.85281374238573*alphaR[12])-94.86832980505142*alphaR[1])*fAvg[17]-94.8683298050514*(alphaR[5]*fAvg[16]+alphaR[4]*(fAvg[14]+fAvg[13])+fAvg[4]*alphaR[13]+fAvg[6]*alphaR[12]+alphaR[5]*fAvg[11]+fAvg[5]*alphaR[11])+((-94.86832980505142*(alphaR[8]+alphaR[7]))-106.0660171779821*alphaR[0])*fAvg[10]-106.0660171779821*(alphaR[1]*fAvg[6]+alphaR[2]*fAvg[5]+fAvg[2]*alphaR[5]+alphaR[3]*fAvg[4]+fAvg[3]*alphaR[4])); 
  Ghat[11] = 0.001190476190476191*((514.3928459844675*(fR[32]+fL[32])-296.98484809835*fR[21]+296.98484809835*fL[21])*amax+(94.86832980505142*alphaR[13]+148.492424049175*alphaR[3])*fAvg[17]+148.492424049175*fAvg[6]*alphaR[13]+118.79393923934*(alphaR[4]*fAvg[12]+fAvg[4]*alphaR[12])+(132.815661727072*alphaR[8]+94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[11]+(132.815661727072*fAvg[8]+94.86832980505142*fAvg[7]+148.492424049175*fAvg[0])*alphaR[11]+132.815661727072*alphaR[5]*fAvg[10]+148.492424049175*(alphaR[2]*fAvg[7]+fAvg[2]*alphaR[7])+132.815661727072*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4])); 
  Ghat[12] = 0.001190476190476191*((514.3928459844675*(fR[34]+fL[34])-296.98484809835*fR[23]+296.98484809835*fL[23])*amax+(132.815661727072*alphaR[13]+148.492424049175*alphaR[3])*fAvg[18]+148.492424049175*alphaR[5]*fAvg[14]+(94.86832980505142*alphaR[8]+132.815661727072*alphaR[7]+148.492424049175*alphaR[0])*fAvg[12]+(94.86832980505142*fAvg[8]+132.815661727072*fAvg[7]+148.492424049175*fAvg[0])*alphaR[12]+118.79393923934*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11])+148.492424049175*(alphaR[1]*fAvg[8]+fAvg[1]*alphaR[8])+132.815661727072*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4])); 
  Ghat[13] = 0.001190476190476191*((514.3928459844675*(fR[35]+fL[35])-296.98484809835*fR[25]+296.98484809835*fL[25])*amax+132.815661727072*alphaR[12]*fAvg[18]+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[17]+118.79393923934*alphaR[5]*fAvg[15]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[13]+(132.815661727072*fAvg[9]+94.86832980505142*fAvg[7]+148.492424049175*fAvg[0])*alphaR[13]+148.492424049175*fAvg[6]*alphaR[11]+132.815661727072*alphaR[4]*fAvg[10]+148.492424049175*(alphaR[3]*fAvg[7]+fAvg[3]*alphaR[7])+132.815661727072*(alphaR[1]*fAvg[5]+fAvg[1]*alphaR[5])); 
  Ghat[14] = 0.001190476190476191*((514.3928459844675*(fR[40]+fL[40])-296.98484809835*fR[27]+296.98484809835*fL[27])*amax+(94.86832980505142*alphaR[12]+148.492424049175*alphaR[1])*fAvg[18]+132.815661727072*alphaR[11]*fAvg[17]+(94.86832980505142*alphaR[8]+148.492424049175*alphaR[0])*fAvg[14]+148.492424049175*(alphaR[5]*fAvg[12]+fAvg[5]*alphaR[12])+132.815661727072*alphaR[4]*fAvg[10]+148.492424049175*(alphaR[3]*fAvg[8]+fAvg[3]*alphaR[8])+132.815661727072*alphaR[2]*fAvg[6]); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[28]+42.42640687119286*fL[28])*amax+(18.97366596101028*alphaR[11]+21.21320343559643*alphaR[2])*fAvg[19]+21.21320343559643*alphaR[4]*fAvg[16]+(18.97366596101028*alphaR[7]+21.21320343559643*alphaR[0])*fAvg[15]+16.97056274847715*(alphaR[5]*fAvg[13]+fAvg[5]*alphaR[13])+21.21320343559643*alphaR[1]*fAvg[9]+18.97366596101028*(alphaR[3]*fAvg[5]+fAvg[3]*alphaR[5])); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[43]+fL[43])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+(18.97366596101028*alphaR[12]+21.21320343559643*alphaR[1])*fAvg[19]+18.97366596101028*alphaR[13]*fAvg[17]+(18.97366596101028*alphaR[8]+21.21320343559643*alphaR[0])*fAvg[16]+21.21320343559643*alphaR[4]*fAvg[15]+18.97366596101028*alphaR[5]*fAvg[10]+21.21320343559643*alphaR[2]*fAvg[9]+18.97366596101028*alphaR[3]*fAvg[6]); 
  Ghat[17] = 2.380952380952381e-4*((2571.964229922338*(fR[44]+fL[44])-1484.92424049175*fR[37]+1484.92424049175*fL[37])*amax+593.9696961967002*(alphaR[5]*fAvg[19]+alphaR[4]*fAvg[18])+(664.0783086353599*alphaR[8]+474.3416490252571*alphaR[7]+742.462120245875*alphaR[0])*fAvg[17]+664.0783086353599*(alphaR[13]*fAvg[16]+alphaR[11]*fAvg[14])+(474.3416490252571*alphaR[11]+742.462120245875*alphaR[2])*fAvg[13]+(474.3416490252571*fAvg[11]+742.462120245875*fAvg[2])*alphaR[13]+593.9696961967002*fAvg[10]*alphaR[12]+742.462120245875*(alphaR[3]*fAvg[11]+fAvg[3]*alphaR[11])+664.0783086353599*alphaR[1]*fAvg[10]+742.462120245875*fAvg[6]*alphaR[7]+664.0783086353599*(alphaR[4]*fAvg[5]+fAvg[4]*alphaR[5])); 
  Ghat[18] = 2.380952380952381e-4*((2571.964229922338*(fR[46]+fL[46])-1484.92424049175*fR[39]+1484.92424049175*fL[39])*amax+(474.3416490252571*alphaR[8]+664.0783086353599*alphaR[7]+742.462120245875*alphaR[0])*fAvg[18]+593.9696961967002*alphaR[4]*fAvg[17]+(474.3416490252571*alphaR[12]+742.462120245875*alphaR[1])*fAvg[14]+664.0783086353599*(alphaR[12]*fAvg[13]+fAvg[12]*alphaR[13])+742.462120245875*(alphaR[3]*fAvg[12]+fAvg[3]*alphaR[12])+fAvg[10]*(593.9696961967002*alphaR[11]+664.0783086353599*alphaR[2])+742.462120245875*(alphaR[5]*fAvg[8]+fAvg[5]*alphaR[8])+664.0783086353599*alphaR[4]*fAvg[6]); 
  Ghat[19] = 0.001666666666666667*((367.4234614174769*(fR[47]+fL[47])-212.1320343559643*fR[42]+212.1320343559643*fL[42])*amax+(94.86832980505142*(alphaR[8]+alphaR[7])+106.0660171779821*alphaR[0])*fAvg[19]+84.85281374238573*alphaR[5]*fAvg[17]+(94.86832980505142*alphaR[12]+106.0660171779822*alphaR[1])*fAvg[16]+(94.86832980505142*alphaR[11]+106.0660171779822*alphaR[2])*fAvg[15]+fAvg[10]*(84.85281374238573*alphaR[13]+94.86832980505142*alphaR[3])+106.0660171779821*alphaR[4]*fAvg[9]+94.86832980505142*alphaR[5]*fAvg[6]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = -1.224744871391589*Ghat[0]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = 0.7071067811865475*Ghat[4]; 
  incr[7] = -1.224744871391589*Ghat[2]; 
  incr[8] = 0.7071067811865475*Ghat[5]; 
  incr[9] = -1.224744871391589*Ghat[3]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 0.7071067811865475*Ghat[7]; 
  incr[12] = 1.58113883008419*Ghat[0]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = 0.7071067811865475*Ghat[10]; 
  incr[18] = -1.224744871391589*Ghat[6]; 
  incr[19] = -1.224744871391589*Ghat[7]; 
  incr[20] = 1.58113883008419*Ghat[1]; 
  incr[21] = 0.7071067811865475*Ghat[11]; 
  incr[22] = 1.58113883008419*Ghat[2]; 
  incr[23] = 0.7071067811865475*Ghat[12]; 
  incr[24] = -1.224744871391589*Ghat[8]; 
  incr[25] = 0.7071067811865475*Ghat[13]; 
  incr[26] = 1.58113883008419*Ghat[3]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = 0.7071067811865475*Ghat[15]; 
  incr[29] = -1.224744871391589*Ghat[9]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = -1.224744871391589*Ghat[11]; 
  incr[33] = 1.58113883008419*Ghat[4]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = -1.224744871391589*Ghat[13]; 
  incr[36] = 1.58113883008419*Ghat[5]; 
  incr[37] = 0.7071067811865475*Ghat[17]; 
  incr[38] = 1.58113883008419*Ghat[6]; 
  incr[39] = 0.7071067811865475*Ghat[18]; 
  incr[40] = -1.224744871391589*Ghat[14]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = 0.7071067811865475*Ghat[19]; 
  incr[43] = -1.224744871391589*Ghat[16]; 
  incr[44] = -1.224744871391589*Ghat[17]; 
  incr[45] = 1.58113883008419*Ghat[10]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

  outR[0] += incr[0]*rdy2R; 
  outR[1] += incr[1]*rdy2R; 
  outR[2] += incr[2]*rdy2R; 
  outR[3] += incr[3]*rdy2R; 
  outR[4] += incr[4]*rdy2R; 
  outR[5] += incr[5]*rdy2R; 
  outR[6] += incr[6]*rdy2R; 
  outR[7] += incr[7]*rdy2R; 
  outR[8] += incr[8]*rdy2R; 
  outR[9] += incr[9]*rdy2R; 
  outR[10] += incr[10]*rdy2R; 
  outR[11] += incr[11]*rdy2R; 
  outR[12] += incr[12]*rdy2R; 
  outR[13] += incr[13]*rdy2R; 
  outR[14] += incr[14]*rdy2R; 
  outR[15] += incr[15]*rdy2R; 
  outR[16] += incr[16]*rdy2R; 
  outR[17] += incr[17]*rdy2R; 
  outR[18] += incr[18]*rdy2R; 
  outR[19] += incr[19]*rdy2R; 
  outR[20] += incr[20]*rdy2R; 
  outR[21] += incr[21]*rdy2R; 
  outR[22] += incr[22]*rdy2R; 
  outR[23] += incr[23]*rdy2R; 
  outR[24] += incr[24]*rdy2R; 
  outR[25] += incr[25]*rdy2R; 
  outR[26] += incr[26]*rdy2R; 
  outR[27] += incr[27]*rdy2R; 
  outR[28] += incr[28]*rdy2R; 
  outR[29] += incr[29]*rdy2R; 
  outR[30] += incr[30]*rdy2R; 
  outR[31] += incr[31]*rdy2R; 
  outR[32] += incr[32]*rdy2R; 
  outR[33] += incr[33]*rdy2R; 
  outR[34] += incr[34]*rdy2R; 
  outR[35] += incr[35]*rdy2R; 
  outR[36] += incr[36]*rdy2R; 
  outR[37] += incr[37]*rdy2R; 
  outR[38] += incr[38]*rdy2R; 
  outR[39] += incr[39]*rdy2R; 
  outR[40] += incr[40]*rdy2R; 
  outR[41] += incr[41]*rdy2R; 
  outR[42] += incr[42]*rdy2R; 
  outR[43] += incr[43]*rdy2R; 
  outR[44] += incr[44]*rdy2R; 
  outR[45] += incr[45]*rdy2R; 
  outR[46] += incr[46]*rdy2R; 
  outR[47] += incr[47]*rdy2R; 

  outL[0] += -1.0*incr[0]*rdy2L; 
  outL[1] += -1.0*incr[1]*rdy2L; 
  outL[2] += incr[2]*rdy2L; 
  outL[3] += -1.0*incr[3]*rdy2L; 
  outL[4] += -1.0*incr[4]*rdy2L; 
  outL[5] += incr[5]*rdy2L; 
  outL[6] += -1.0*incr[6]*rdy2L; 
  outL[7] += incr[7]*rdy2L; 
  outL[8] += -1.0*incr[8]*rdy2L; 
  outL[9] += incr[9]*rdy2L; 
  outL[10] += -1.0*incr[10]*rdy2L; 
  outL[11] += -1.0*incr[11]*rdy2L; 
  outL[12] += -1.0*incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += -1.0*incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 
  outL[16] += incr[16]*rdy2L; 
  outL[17] += -1.0*incr[17]*rdy2L; 
  outL[18] += incr[18]*rdy2L; 
  outL[19] += incr[19]*rdy2L; 
  outL[20] += -1.0*incr[20]*rdy2L; 
  outL[21] += -1.0*incr[21]*rdy2L; 
  outL[22] += -1.0*incr[22]*rdy2L; 
  outL[23] += -1.0*incr[23]*rdy2L; 
  outL[24] += incr[24]*rdy2L; 
  outL[25] += -1.0*incr[25]*rdy2L; 
  outL[26] += -1.0*incr[26]*rdy2L; 
  outL[27] += -1.0*incr[27]*rdy2L; 
  outL[28] += -1.0*incr[28]*rdy2L; 
  outL[29] += incr[29]*rdy2L; 
  outL[30] += -1.0*incr[30]*rdy2L; 
  outL[31] += incr[31]*rdy2L; 
  outL[32] += incr[32]*rdy2L; 
  outL[33] += -1.0*incr[33]*rdy2L; 
  outL[34] += incr[34]*rdy2L; 
  outL[35] += incr[35]*rdy2L; 
  outL[36] += -1.0*incr[36]*rdy2L; 
  outL[37] += -1.0*incr[37]*rdy2L; 
  outL[38] += -1.0*incr[38]*rdy2L; 
  outL[39] += -1.0*incr[39]*rdy2L; 
  outL[40] += incr[40]*rdy2L; 
  outL[41] += incr[41]*rdy2L; 
  outL[42] += -1.0*incr[42]*rdy2L; 
  outL[43] += incr[43]*rdy2L; 
  outL[44] += incr[44]*rdy2L; 
  outL[45] += -1.0*incr[45]*rdy2L; 
  outL[46] += incr[46]*rdy2L; 
  outL[47] += incr[47]*rdy2L; 
return std::abs(alphaSurfAvgR); 
} 
double GyrokineticGenGeoSurf2x2vSer_vpar_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[8] = (1.154700538379252*bmag[1])/rdmu2R; 
  hamilR[11] = 2.0*(bmag[4]*wmuR+phi[4]*q_); 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 
  hamilR[25] = (1.154700538379251*bmag[4])/rdmu2R; 

  double BstarXdBmagR[48]; 

  double BstarYdBmagR[48]; 
  BstarYdBmagR[0] = -(1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[1] = -(1.732050807568877*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[3] = -(1.0*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[6] = -(1.0*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[11] = -(1.732050807568877*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[21] = -(1.0*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[20]; 
  alphaR[0] = (0.03535533905932736*(hamilR[19]*(15.0*BstarYdBmagR[21]-8.660254037844387*BstarYdBmagR[11])+hamilR[5]*(15.0*BstarYdBmagR[6]-8.660254037844386*BstarYdBmagR[1])+hamilR[2]*(15.0*BstarYdBmagR[3]-8.660254037844386*BstarYdBmagR[0]))*rdy2R)/m_; 
  alphaR[1] = (0.03535533905932736*(13.41640786499874*hamilR[5]*BstarYdBmagR[21]+(13.41640786499874*BstarYdBmagR[6]-7.745966692414834*BstarYdBmagR[1])*hamilR[19]-7.745966692414834*hamilR[5]*BstarYdBmagR[11]+15.0*hamilR[2]*BstarYdBmagR[6]+(15.0*BstarYdBmagR[3]-8.660254037844386*BstarYdBmagR[0])*hamilR[5]-8.660254037844386*BstarYdBmagR[1]*hamilR[2])*rdy2R)/m_; 
  alphaR[2] = (0.1767766952966368*((6.708203932499369*BstarYdBmagR[6]-3.872983346207417*BstarYdBmagR[1])*hamilR[20]+(6.708203932499369*BstarYdBmagR[3]-3.872983346207417*BstarYdBmagR[0])*hamilR[12])*rdy2R)/m_; 
  alphaR[4] = (0.03535533905932736*(hamilR[20]*(30.0*BstarYdBmagR[21]-17.32050807568877*BstarYdBmagR[11]+33.54101966249684*BstarYdBmagR[3]-19.36491673103708*BstarYdBmagR[0])+(33.54101966249685*BstarYdBmagR[6]-19.36491673103709*BstarYdBmagR[1])*hamilR[12])*rdy2R)/m_; 
  alphaR[7] = (0.005050762722761052*((67.0820393249937*hamilR[19]+105.0*hamilR[2])*BstarYdBmagR[21]+((-38.72983346207417*BstarYdBmagR[11])+105.0*BstarYdBmagR[3]-60.62177826491071*BstarYdBmagR[0])*hamilR[19]-60.6217782649107*hamilR[2]*BstarYdBmagR[11]+hamilR[5]*(93.91485505499116*BstarYdBmagR[6]-54.22176684690384*BstarYdBmagR[1]))*rdy2R)/m_; 
  alphaR[11] = (0.1767766952966368*(6.708203932499369*hamilR[12]*BstarYdBmagR[21]+(6.0*BstarYdBmagR[6]-3.464101615137754*BstarYdBmagR[1])*hamilR[20]-3.872983346207417*BstarYdBmagR[11]*hamilR[12])*rdy2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.00625*(15.0*hamilR[19]*BstarYdBmagR[21]-8.660254037844387*BstarYdBmagR[11]*hamilR[19]+15.0*hamilR[5]*BstarYdBmagR[6]-8.660254037844386*BstarYdBmagR[1]*hamilR[5]+15.0*hamilR[2]*BstarYdBmagR[3]-8.660254037844386*BstarYdBmagR[0]*hamilR[2])*rdy2R)/m_; 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[13]+fL[13])+1.732050807568877*(fL[3]-1.0*fR[3])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[23]+fL[23])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[24]+fL[24])+3.0*(fL[7]-1.0*fR[7]))+3.0*(fR[2]+fL[2])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[27]+fL[27])+3.0*(fL[10]-1.0*fR[10]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[34]+fL[34])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[5]+fL[5]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[39]+fL[39])+1.732050807568877*(fL[17]-1.0*fR[17])+fR[8]+fL[8]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[40]+fL[40])+1.732050807568877*(fL[18]-1.0*fR[18])+fR[9]+fL[9]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[21]-1.0*(8.660254037844387*fL[21]+5.0*(fR[11]+fL[11]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[22]-1.0*(8.660254037844387*fL[22]+5.0*(fR[12]+fL[12]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[30]-1.0*(8.660254037844387*fL[30]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[46]+fL[46])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[16]+fL[16])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[32]-1.0*(8.660254037844387*fL[32]+5.0*(fR[19]+fL[19]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[33]-1.0*(8.660254037844387*fL[33]+5.0*(fR[20]+fL[20]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[37]-1.0*(8.660254037844387*fL[37]+5.0*(fR[25]+fL[25]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[38]-1.0*(8.660254037844387*fL[38]+5.0*(fR[26]+fL[26]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[42]-1.0*(8.660254037844387*fL[42]+5.0*(fR[28]+fL[28]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[43]-1.0*(8.660254037844387*fL[43]+5.0*(fR[29]+fL[29]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[44]-1.0*(8.660254037844387*fL[44]+5.0*(fR[35]+fL[35]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[45]-1.0*(8.660254037844387*fL[45]+5.0*(fR[36]+fL[36]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[41]+fL[41]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[13]-6.324555320336761*fL[13]-4.898979485566357*(fR[3]+fL[3])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[11]*fAvg[11]+alphaR[7]*fAvg[7]+alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2]+alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[23]-94.8683298050514*fL[23]-73.48469228349535*(fR[6]+fL[6])+42.42640687119286*fR[1]-42.42640687119286*fL[1])*amax-18.97366596101028*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11]+alphaR[1]*fAvg[7]+fAvg[1]*alphaR[7])-21.21320343559643*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4]+alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[24]-94.8683298050514*fL[24]-73.48469228349535*(fR[7]+fL[7])+42.42640687119286*fR[2]-42.42640687119286*fL[2])*amax-18.97366596101028*alphaR[4]*fAvg[12]-21.21320343559643*(alphaR[7]*fAvg[11]+fAvg[7]*alphaR[11])-18.97366596101028*alphaR[2]*fAvg[8]-21.21320343559643*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4]+alphaR[0]*fAvg[2]+fAvg[0]*alphaR[2])); 
  Ghat[3] = -0.008333333333333333*((94.8683298050514*fR[27]-94.8683298050514*fL[27]-73.48469228349535*(fR[10]+fL[10])+42.42640687119286*fR[4]-42.42640687119286*fL[4])*amax-21.21320343559643*(alphaR[11]*fAvg[17]+alphaR[7]*fAvg[13])-21.21320343559643*(alphaR[4]*fAvg[10]+alphaR[2]*fAvg[6]+alphaR[1]*fAvg[5]+alphaR[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[34]-94.86832980505142*fL[34]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[5]-42.42640687119286*fL[5])*amax+((-16.97056274847715*alphaR[11])-18.97366596101028*alphaR[2])*fAvg[12]-18.97366596101028*(alphaR[1]*fAvg[11]+fAvg[1]*alphaR[11]+alphaR[4]*(fAvg[8]+fAvg[7])+fAvg[4]*alphaR[7])-21.21320343559643*(alphaR[0]*fAvg[4]+fAvg[0]*alphaR[4]+alphaR[1]*fAvg[2]+fAvg[1]*alphaR[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[39]-94.86832980505142*fL[39]-73.48469228349535*(fR[17]+fL[17])+42.42640687119286*fR[8]-42.42640687119286*fL[8])*amax-18.97366596101028*(alphaR[4]*fAvg[17]+alphaR[1]*fAvg[13])+fAvg[10]*((-18.97366596101028*alphaR[11])-21.21320343559643*alphaR[2])-18.97366596101028*fAvg[5]*alphaR[7]-21.21320343559643*(alphaR[4]*fAvg[6]+alphaR[0]*fAvg[5]+alphaR[1]*fAvg[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[40]-94.86832980505142*fL[40]-73.48469228349535*(fR[18]+fL[18])+42.42640687119286*fR[9]-42.42640687119286*fL[9])*amax-18.97366596101028*alphaR[4]*fAvg[18]-21.21320343559643*alphaR[7]*fAvg[17]-18.97366596101028*alphaR[2]*fAvg[14]-21.21320343559643*(alphaR[11]*fAvg[13]+alphaR[1]*fAvg[10]+alphaR[0]*fAvg[6]+alphaR[4]*fAvg[5]+alphaR[2]*fAvg[3])); 
  Ghat[7] = 0.001190476190476191*((514.3928459844675*(fR[21]+fL[21])-296.98484809835*fR[11]+296.98484809835*fL[11])*amax+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[11]+148.492424049175*fAvg[2]*alphaR[11]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[7]+148.492424049175*fAvg[0]*alphaR[7]+132.815661727072*(alphaR[4]*fAvg[4]+alphaR[1]*fAvg[1])); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[22]+fL[22])-42.42640687119286*fR[12]+42.42640687119286*fL[12])*amax+21.21320343559643*alphaR[1]*fAvg[12]+18.97366596101028*alphaR[11]*fAvg[11]+21.21320343559643*alphaR[0]*fAvg[8]+18.97366596101028*(alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[30]+fL[30])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[4]*fAvg[19]+21.21320343559643*(alphaR[2]*fAvg[16]+alphaR[1]*fAvg[15])+21.21320343559643*alphaR[0]*fAvg[9]); 
  Ghat[10] = -0.001666666666666667*((474.341649025257*fR[46]-474.341649025257*fL[46]-367.4234614174767*(fR[31]+fL[31])+212.1320343559643*fR[16]-212.1320343559643*fL[16])*amax+((-84.85281374238573*alphaR[11])-94.86832980505142*alphaR[2])*fAvg[18]-94.86832980505142*alphaR[1]*fAvg[17]-94.8683298050514*(alphaR[4]*(fAvg[14]+fAvg[13])+fAvg[5]*alphaR[11])+((-94.86832980505142*alphaR[7])-106.0660171779821*alphaR[0])*fAvg[10]-106.0660171779821*(alphaR[1]*fAvg[6]+alphaR[2]*fAvg[5]+fAvg[3]*alphaR[4])); 
  Ghat[11] = 0.001190476190476191*((514.3928459844675*(fR[32]+fL[32])-296.98484809835*fR[19]+296.98484809835*fL[19])*amax+118.79393923934*alphaR[4]*fAvg[12]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[11]+(132.815661727072*fAvg[8]+94.86832980505142*fAvg[7]+148.492424049175*fAvg[0])*alphaR[11]+148.492424049175*(alphaR[2]*fAvg[7]+fAvg[2]*alphaR[7])+132.815661727072*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4])); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[33]+fL[33])-42.42640687119286*fR[20]+42.42640687119286*fL[20])*amax+(18.97366596101028*alphaR[7]+21.21320343559643*alphaR[0])*fAvg[12]+16.97056274847715*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11])+21.21320343559643*alphaR[1]*fAvg[8]+18.97366596101028*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4])); 
  Ghat[13] = 0.001190476190476191*((514.3928459844675*(fR[37]+fL[37])-296.98484809835*fR[25]+296.98484809835*fL[25])*amax+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[17]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[13]+148.492424049175*fAvg[6]*alphaR[11]+132.815661727072*alphaR[4]*fAvg[10]+148.492424049175*fAvg[3]*alphaR[7]+132.815661727072*alphaR[1]*fAvg[5]); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[38]+fL[38])-42.42640687119286*fR[26]+42.42640687119286*fL[26])*amax+21.21320343559643*alphaR[1]*fAvg[18]+18.97366596101028*alphaR[11]*fAvg[17]+21.21320343559643*alphaR[0]*fAvg[14]+18.97366596101028*(alphaR[4]*fAvg[10]+alphaR[2]*fAvg[6])); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[42]+fL[42])-42.42640687119286*fR[28]+42.42640687119286*fL[28])*amax+(18.97366596101028*alphaR[11]+21.21320343559643*alphaR[2])*fAvg[19]+21.21320343559643*alphaR[4]*fAvg[16]+(18.97366596101028*alphaR[7]+21.21320343559643*alphaR[0])*fAvg[15]+21.21320343559643*alphaR[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[43]+fL[43])-42.42640687119286*fR[29]+42.42640687119286*fL[29])*amax+21.21320343559643*alphaR[1]*fAvg[19]+21.21320343559643*(alphaR[0]*fAvg[16]+alphaR[4]*fAvg[15])+21.21320343559643*alphaR[2]*fAvg[9]); 
  Ghat[17] = 0.001190476190476191*((514.3928459844675*(fR[44]+fL[44])-296.98484809835*fR[35]+296.98484809835*fL[35])*amax+118.79393923934*alphaR[4]*fAvg[18]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[17]+132.815661727072*alphaR[11]*fAvg[14]+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[13]+148.492424049175*fAvg[3]*alphaR[11]+132.815661727072*alphaR[1]*fAvg[10]+148.492424049175*fAvg[6]*alphaR[7]+132.815661727072*alphaR[4]*fAvg[5]); 
  Ghat[18] = 0.001666666666666667*((367.4234614174769*(fR[45]+fL[45])-212.1320343559643*fR[36]+212.1320343559643*fL[36])*amax+(94.86832980505142*alphaR[7]+106.0660171779821*alphaR[0])*fAvg[18]+84.85281374238573*alphaR[4]*fAvg[17]+106.0660171779822*alphaR[1]*fAvg[14]+84.85281374238573*fAvg[10]*alphaR[11]+94.86832980505142*(alphaR[2]*fAvg[10]+alphaR[4]*fAvg[6])); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[41]+42.42640687119286*fL[41])*amax+(18.97366596101028*alphaR[7]+21.21320343559643*alphaR[0])*fAvg[19]+21.21320343559643*alphaR[1]*fAvg[16]+(18.97366596101028*alphaR[11]+21.21320343559643*alphaR[2])*fAvg[15]+21.21320343559643*alphaR[4]*fAvg[9]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = 0.7071067811865475*Ghat[2]; 
  incr[3] = -1.224744871391589*Ghat[0]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = 0.7071067811865475*Ghat[4]; 
  incr[6] = -1.224744871391589*Ghat[1]; 
  incr[7] = -1.224744871391589*Ghat[2]; 
  incr[8] = 0.7071067811865475*Ghat[5]; 
  incr[9] = 0.7071067811865475*Ghat[6]; 
  incr[10] = -1.224744871391589*Ghat[3]; 
  incr[11] = 0.7071067811865475*Ghat[7]; 
  incr[12] = 0.7071067811865475*Ghat[8]; 
  incr[13] = 1.58113883008419*Ghat[0]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = 0.7071067811865475*Ghat[10]; 
  incr[17] = -1.224744871391589*Ghat[5]; 
  incr[18] = -1.224744871391589*Ghat[6]; 
  incr[19] = 0.7071067811865475*Ghat[11]; 
  incr[20] = 0.7071067811865475*Ghat[12]; 
  incr[21] = -1.224744871391589*Ghat[7]; 
  incr[22] = -1.224744871391589*Ghat[8]; 
  incr[23] = 1.58113883008419*Ghat[1]; 
  incr[24] = 1.58113883008419*Ghat[2]; 
  incr[25] = 0.7071067811865475*Ghat[13]; 
  incr[26] = 0.7071067811865475*Ghat[14]; 
  incr[27] = 1.58113883008419*Ghat[3]; 
  incr[28] = 0.7071067811865475*Ghat[15]; 
  incr[29] = 0.7071067811865475*Ghat[16]; 
  incr[30] = -1.224744871391589*Ghat[9]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = -1.224744871391589*Ghat[11]; 
  incr[33] = -1.224744871391589*Ghat[12]; 
  incr[34] = 1.58113883008419*Ghat[4]; 
  incr[35] = 0.7071067811865475*Ghat[17]; 
  incr[36] = 0.7071067811865475*Ghat[18]; 
  incr[37] = -1.224744871391589*Ghat[13]; 
  incr[38] = -1.224744871391589*Ghat[14]; 
  incr[39] = 1.58113883008419*Ghat[5]; 
  incr[40] = 1.58113883008419*Ghat[6]; 
  incr[41] = 0.7071067811865475*Ghat[19]; 
  incr[42] = -1.224744871391589*Ghat[15]; 
  incr[43] = -1.224744871391589*Ghat[16]; 
  incr[44] = -1.224744871391589*Ghat[17]; 
  incr[45] = -1.224744871391589*Ghat[18]; 
  incr[46] = 1.58113883008419*Ghat[10]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 
  outR[8] += incr[8]*rdvpar2R; 
  outR[9] += incr[9]*rdvpar2R; 
  outR[10] += incr[10]*rdvpar2R; 
  outR[11] += incr[11]*rdvpar2R; 
  outR[12] += incr[12]*rdvpar2R; 
  outR[13] += incr[13]*rdvpar2R; 
  outR[14] += incr[14]*rdvpar2R; 
  outR[15] += incr[15]*rdvpar2R; 
  outR[16] += incr[16]*rdvpar2R; 
  outR[17] += incr[17]*rdvpar2R; 
  outR[18] += incr[18]*rdvpar2R; 
  outR[19] += incr[19]*rdvpar2R; 
  outR[20] += incr[20]*rdvpar2R; 
  outR[21] += incr[21]*rdvpar2R; 
  outR[22] += incr[22]*rdvpar2R; 
  outR[23] += incr[23]*rdvpar2R; 
  outR[24] += incr[24]*rdvpar2R; 
  outR[25] += incr[25]*rdvpar2R; 
  outR[26] += incr[26]*rdvpar2R; 
  outR[27] += incr[27]*rdvpar2R; 
  outR[28] += incr[28]*rdvpar2R; 
  outR[29] += incr[29]*rdvpar2R; 
  outR[30] += incr[30]*rdvpar2R; 
  outR[31] += incr[31]*rdvpar2R; 
  outR[32] += incr[32]*rdvpar2R; 
  outR[33] += incr[33]*rdvpar2R; 
  outR[34] += incr[34]*rdvpar2R; 
  outR[35] += incr[35]*rdvpar2R; 
  outR[36] += incr[36]*rdvpar2R; 
  outR[37] += incr[37]*rdvpar2R; 
  outR[38] += incr[38]*rdvpar2R; 
  outR[39] += incr[39]*rdvpar2R; 
  outR[40] += incr[40]*rdvpar2R; 
  outR[41] += incr[41]*rdvpar2R; 
  outR[42] += incr[42]*rdvpar2R; 
  outR[43] += incr[43]*rdvpar2R; 
  outR[44] += incr[44]*rdvpar2R; 
  outR[45] += incr[45]*rdvpar2R; 
  outR[46] += incr[46]*rdvpar2R; 
  outR[47] += incr[47]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += -1.0*incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 
  outL[4] += -1.0*incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 
  outL[8] += -1.0*incr[8]*rdvpar2L; 
  outL[9] += -1.0*incr[9]*rdvpar2L; 
  outL[10] += incr[10]*rdvpar2L; 
  outL[11] += -1.0*incr[11]*rdvpar2L; 
  outL[12] += -1.0*incr[12]*rdvpar2L; 
  outL[13] += -1.0*incr[13]*rdvpar2L; 
  outL[14] += -1.0*incr[14]*rdvpar2L; 
  outL[15] += incr[15]*rdvpar2L; 
  outL[16] += -1.0*incr[16]*rdvpar2L; 
  outL[17] += incr[17]*rdvpar2L; 
  outL[18] += incr[18]*rdvpar2L; 
  outL[19] += -1.0*incr[19]*rdvpar2L; 
  outL[20] += -1.0*incr[20]*rdvpar2L; 
  outL[21] += incr[21]*rdvpar2L; 
  outL[22] += incr[22]*rdvpar2L; 
  outL[23] += -1.0*incr[23]*rdvpar2L; 
  outL[24] += -1.0*incr[24]*rdvpar2L; 
  outL[25] += -1.0*incr[25]*rdvpar2L; 
  outL[26] += -1.0*incr[26]*rdvpar2L; 
  outL[27] += -1.0*incr[27]*rdvpar2L; 
  outL[28] += -1.0*incr[28]*rdvpar2L; 
  outL[29] += -1.0*incr[29]*rdvpar2L; 
  outL[30] += incr[30]*rdvpar2L; 
  outL[31] += incr[31]*rdvpar2L; 
  outL[32] += incr[32]*rdvpar2L; 
  outL[33] += incr[33]*rdvpar2L; 
  outL[34] += -1.0*incr[34]*rdvpar2L; 
  outL[35] += -1.0*incr[35]*rdvpar2L; 
  outL[36] += -1.0*incr[36]*rdvpar2L; 
  outL[37] += incr[37]*rdvpar2L; 
  outL[38] += incr[38]*rdvpar2L; 
  outL[39] += -1.0*incr[39]*rdvpar2L; 
  outL[40] += -1.0*incr[40]*rdvpar2L; 
  outL[41] += -1.0*incr[41]*rdvpar2L; 
  outL[42] += incr[42]*rdvpar2L; 
  outL[43] += incr[43]*rdvpar2L; 
  outL[44] += incr[44]*rdvpar2L; 
  outL[45] += incr[45]*rdvpar2L; 
  outL[46] += -1.0*incr[46]*rdvpar2L; 
  outL[47] += incr[47]*rdvpar2L; 
return std::abs(alphaSurfAvgR); 
} 
double GyrokineticGenGeoSurf2x2vSer_x_P2_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[11] = 2.0*phi[4]*q_; 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 

  double BstarXdBmagR[48]; 

  double alphaR[20]; 
  alphaR[0] = -(0.1767766952966368*b_z[0]*jacobTotInv[0]*(3.872983346207417*hamilR[19]-3.0*hamilR[5]+1.732050807568877*hamilR[2])*rdy2R)/q_; 
  alphaR[1] = (0.1767766952966368*b_z[0]*jacobTotInv[0]*(6.708203932499369*hamilR[20]-3.872983346207417*hamilR[12])*rdy2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(3.872983346207417*b_z[0]*jacobTotInv[0]*hamilR[19]-3.0*b_z[0]*jacobTotInv[0]*hamilR[5]+1.732050807568877*b_z[0]*jacobTotInv[0]*hamilR[2])*rdy2R)/q_; 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[11]+fL[11])+1.732050807568877*(fL[1]-1.0*fR[1])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[19]+fL[19])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[2]+fL[2])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[21]+fL[21])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[25]+fL[25])+3.0*(fL[8]-1.0*fR[8]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[32]+fL[32])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[7]+fL[7]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[35]+fL[35])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[9]+fL[9]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[37]+fL[37])+1.732050807568877*(fL[17]-1.0*fR[17])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[20]-1.0*(8.660254037844387*fL[20]+5.0*(fR[12]+fL[12]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[23]-1.0*(8.660254037844387*fL[23]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[28]-1.0*(8.660254037844387*fL[28]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[44]+fL[44])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[18]+fL[18])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[33]-1.0*(8.660254037844387*fL[33]+5.0*(fR[22]+fL[22]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[24]+fL[24]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[36]-1.0*(8.660254037844387*fL[36]+5.0*(fR[26]+fL[26]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[39]-1.0*(8.660254037844387*fL[39]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[29]+fL[29]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[42]-1.0*(8.660254037844387*fL[42]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[45]-1.0*(8.660254037844387*fL[45]+5.0*(fR[38]+fL[38]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[40]+fL[40]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[43]+fL[43]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[11]-6.324555320336761*fL[11]-4.898979485566357*(fR[1]+fL[1])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[19]-94.8683298050514*fL[19]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[2]-42.42640687119286*fL[2])*amax-18.97366596101028*alphaR[1]*fAvg[7]-21.21320343559643*(alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.04166666666666666*((18.97366596101028*fR[21]-18.97366596101028*fL[21]-14.69693845669907*(fR[6]+fL[6])+8.485281374238571*fR[3]-8.485281374238571*fL[3])*amax-4.242640687119286*(alphaR[1]*fAvg[4]+alphaR[0]*fAvg[2])); 
  Ghat[3] = -0.04166666666666666*((18.97366596101028*fR[25]-18.97366596101028*fL[25]-14.69693845669907*(fR[8]+fL[8])+8.485281374238571*fR[4]-8.485281374238571*fL[4])*amax-4.242640687119286*(alphaR[1]*fAvg[5]+alphaR[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[32]-94.86832980505142*fL[32]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[7]-42.42640687119286*fL[7])*amax-18.97366596101028*alphaR[1]*fAvg[11]-21.21320343559643*(alphaR[0]*fAvg[4]+alphaR[1]*fAvg[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[35]-94.86832980505142*fL[35]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[9]-42.42640687119286*fL[9])*amax-18.97366596101028*alphaR[1]*fAvg[13]-21.21320343559643*(alphaR[0]*fAvg[5]+alphaR[1]*fAvg[3])); 
  Ghat[6] = -0.125*((6.324555320336761*fR[37]-6.324555320336761*fL[37]-4.898979485566357*(fR[17]+fL[17])+2.828427124746191*fR[10]-2.828427124746191*fL[10])*amax-1.414213562373095*(alphaR[1]*fAvg[10]+alphaR[0]*fAvg[6])); 
  Ghat[7] = 0.025*((24.49489742783179*(fR[20]+fL[20])-14.14213562373095*fR[12]+14.14213562373095*fL[12])*amax+7.071067811865476*alphaR[0]*fAvg[7]+6.324555320336761*alphaR[1]*fAvg[1]); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[23]+fL[23])-42.42640687119286*fR[13]+42.42640687119286*fL[13])*amax+21.21320343559643*alphaR[1]*fAvg[12]+21.21320343559643*alphaR[0]*fAvg[8]); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[28]+fL[28])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[1]*fAvg[15]+21.21320343559643*alphaR[0]*fAvg[9]); 
  Ghat[10] = -0.008333333333333333*((94.8683298050514*fR[44]-94.8683298050514*fL[44]-73.48469228349535*(fR[31]+fL[31])+42.42640687119286*fR[18]-42.42640687119286*fL[18])*amax-18.97366596101028*alphaR[1]*fAvg[17]-21.21320343559643*(alphaR[0]*fAvg[10]+alphaR[1]*fAvg[6])); 
  Ghat[11] = 0.008333333333333333*((73.48469228349536*(fR[33]+fL[33])-42.42640687119286*fR[22]+42.42640687119286*fL[22])*amax+21.21320343559643*alphaR[0]*fAvg[11]+18.97366596101028*alphaR[1]*fAvg[4]); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[34]+fL[34])-42.42640687119286*fR[24]+42.42640687119286*fL[24])*amax+21.21320343559643*alphaR[0]*fAvg[12]+21.21320343559643*alphaR[1]*fAvg[8]); 
  Ghat[13] = 0.008333333333333333*((73.48469228349536*(fR[36]+fL[36])-42.42640687119286*fR[26]+42.42640687119286*fL[26])*amax+21.21320343559643*alphaR[0]*fAvg[13]+18.97366596101028*alphaR[1]*fAvg[5]); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[39]+fL[39])-42.42640687119286*fR[27]+42.42640687119286*fL[27])*amax+21.21320343559643*alphaR[1]*fAvg[18]+21.21320343559643*alphaR[0]*fAvg[14]); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[29]+42.42640687119286*fL[29])*amax+21.21320343559643*alphaR[0]*fAvg[15]+21.21320343559643*alphaR[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[42]+fL[42])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+21.21320343559643*alphaR[1]*fAvg[19]+21.21320343559643*alphaR[0]*fAvg[16]); 
  Ghat[17] = 0.025*((24.49489742783179*(fR[45]+fL[45])-14.14213562373095*fR[38]+14.14213562373095*fL[38])*amax+7.071067811865476*alphaR[0]*fAvg[17]+6.324555320336761*alphaR[1]*fAvg[10]); 
  Ghat[18] = 0.008333333333333333*((73.48469228349536*(fR[46]+fL[46])-42.42640687119286*fR[40]+42.42640687119286*fL[40])*amax+21.21320343559643*alphaR[0]*fAvg[18]+21.21320343559643*alphaR[1]*fAvg[14]); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[43]+42.42640687119286*fL[43])*amax+21.21320343559643*alphaR[0]*fAvg[19]+21.21320343559643*alphaR[1]*fAvg[16]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = -1.224744871391589*Ghat[0]; 
  incr[2] = 0.7071067811865475*Ghat[1]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = -1.224744871391589*Ghat[2]; 
  incr[7] = 0.7071067811865475*Ghat[4]; 
  incr[8] = -1.224744871391589*Ghat[3]; 
  incr[9] = 0.7071067811865475*Ghat[5]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 1.58113883008419*Ghat[0]; 
  incr[12] = 0.7071067811865475*Ghat[7]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = -1.224744871391589*Ghat[6]; 
  incr[18] = 0.7071067811865475*Ghat[10]; 
  incr[19] = 1.58113883008419*Ghat[1]; 
  incr[20] = -1.224744871391589*Ghat[7]; 
  incr[21] = 1.58113883008419*Ghat[2]; 
  incr[22] = 0.7071067811865475*Ghat[11]; 
  incr[23] = -1.224744871391589*Ghat[8]; 
  incr[24] = 0.7071067811865475*Ghat[12]; 
  incr[25] = 1.58113883008419*Ghat[3]; 
  incr[26] = 0.7071067811865475*Ghat[13]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = -1.224744871391589*Ghat[9]; 
  incr[29] = 0.7071067811865475*Ghat[15]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = 1.58113883008419*Ghat[4]; 
  incr[33] = -1.224744871391589*Ghat[11]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = 1.58113883008419*Ghat[5]; 
  incr[36] = -1.224744871391589*Ghat[13]; 
  incr[37] = 1.58113883008419*Ghat[6]; 
  incr[38] = 0.7071067811865475*Ghat[17]; 
  incr[39] = -1.224744871391589*Ghat[14]; 
  incr[40] = 0.7071067811865475*Ghat[18]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = -1.224744871391589*Ghat[16]; 
  incr[43] = 0.7071067811865475*Ghat[19]; 
  incr[44] = 1.58113883008419*Ghat[10]; 
  incr[45] = -1.224744871391589*Ghat[17]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 
  outR[8] += incr[8]*rdx2R; 
  outR[9] += incr[9]*rdx2R; 
  outR[10] += incr[10]*rdx2R; 
  outR[11] += incr[11]*rdx2R; 
  outR[12] += incr[12]*rdx2R; 
  outR[13] += incr[13]*rdx2R; 
  outR[14] += incr[14]*rdx2R; 
  outR[15] += incr[15]*rdx2R; 
  outR[16] += incr[16]*rdx2R; 
  outR[17] += incr[17]*rdx2R; 
  outR[18] += incr[18]*rdx2R; 
  outR[19] += incr[19]*rdx2R; 
  outR[20] += incr[20]*rdx2R; 
  outR[21] += incr[21]*rdx2R; 
  outR[22] += incr[22]*rdx2R; 
  outR[23] += incr[23]*rdx2R; 
  outR[24] += incr[24]*rdx2R; 
  outR[25] += incr[25]*rdx2R; 
  outR[26] += incr[26]*rdx2R; 
  outR[27] += incr[27]*rdx2R; 
  outR[28] += incr[28]*rdx2R; 
  outR[29] += incr[29]*rdx2R; 
  outR[30] += incr[30]*rdx2R; 
  outR[31] += incr[31]*rdx2R; 
  outR[32] += incr[32]*rdx2R; 
  outR[33] += incr[33]*rdx2R; 
  outR[34] += incr[34]*rdx2R; 
  outR[35] += incr[35]*rdx2R; 
  outR[36] += incr[36]*rdx2R; 
  outR[37] += incr[37]*rdx2R; 
  outR[38] += incr[38]*rdx2R; 
  outR[39] += incr[39]*rdx2R; 
  outR[40] += incr[40]*rdx2R; 
  outR[41] += incr[41]*rdx2R; 
  outR[42] += incr[42]*rdx2R; 
  outR[43] += incr[43]*rdx2R; 
  outR[44] += incr[44]*rdx2R; 
  outR[45] += incr[45]*rdx2R; 
  outR[46] += incr[46]*rdx2R; 
  outR[47] += incr[47]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += -1.0*incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += incr[6]*rdx2L; 
  outL[7] += -1.0*incr[7]*rdx2L; 
  outL[8] += incr[8]*rdx2L; 
  outL[9] += -1.0*incr[9]*rdx2L; 
  outL[10] += -1.0*incr[10]*rdx2L; 
  outL[11] += -1.0*incr[11]*rdx2L; 
  outL[12] += -1.0*incr[12]*rdx2L; 
  outL[13] += -1.0*incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 
  outL[16] += incr[16]*rdx2L; 
  outL[17] += incr[17]*rdx2L; 
  outL[18] += -1.0*incr[18]*rdx2L; 
  outL[19] += -1.0*incr[19]*rdx2L; 
  outL[20] += incr[20]*rdx2L; 
  outL[21] += -1.0*incr[21]*rdx2L; 
  outL[22] += -1.0*incr[22]*rdx2L; 
  outL[23] += incr[23]*rdx2L; 
  outL[24] += -1.0*incr[24]*rdx2L; 
  outL[25] += -1.0*incr[25]*rdx2L; 
  outL[26] += -1.0*incr[26]*rdx2L; 
  outL[27] += -1.0*incr[27]*rdx2L; 
  outL[28] += incr[28]*rdx2L; 
  outL[29] += -1.0*incr[29]*rdx2L; 
  outL[30] += -1.0*incr[30]*rdx2L; 
  outL[31] += incr[31]*rdx2L; 
  outL[32] += -1.0*incr[32]*rdx2L; 
  outL[33] += incr[33]*rdx2L; 
  outL[34] += incr[34]*rdx2L; 
  outL[35] += -1.0*incr[35]*rdx2L; 
  outL[36] += incr[36]*rdx2L; 
  outL[37] += -1.0*incr[37]*rdx2L; 
  outL[38] += -1.0*incr[38]*rdx2L; 
  outL[39] += incr[39]*rdx2L; 
  outL[40] += -1.0*incr[40]*rdx2L; 
  outL[41] += incr[41]*rdx2L; 
  outL[42] += incr[42]*rdx2L; 
  outL[43] += -1.0*incr[43]*rdx2L; 
  outL[44] += -1.0*incr[44]*rdx2L; 
  outL[45] += incr[45]*rdx2L; 
  outL[46] += incr[46]*rdx2L; 
  outL[47] += incr[47]*rdx2L; 
return std::abs(alphaSurfAvgR); 
} 
double GyrokineticGenGeoSurf2x2vSer_y_P2_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[11] = 2.0*phi[4]*q_; 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 

  double BstarYdBmagR[48]; 

  double alphaR[20]; 
  alphaR[0] = (0.1767766952966368*b_z[0]*jacobTotInv[0]*(3.872983346207417*hamilR[20]-3.0*hamilR[5]+1.732050807568877*hamilR[1])*rdx2R)/q_; 
  alphaR[1] = -(0.1767766952966368*b_z[0]*jacobTotInv[0]*(6.708203932499369*hamilR[19]-3.872983346207417*hamilR[11])*rdx2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*(3.872983346207417*b_z[0]*jacobTotInv[0]*hamilR[20]-3.0*b_z[0]*jacobTotInv[0]*hamilR[5]+1.732050807568877*b_z[0]*jacobTotInv[0]*hamilR[1])*rdx2R)/q_; 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[12]+fL[12])+1.732050807568877*(fL[2]-1.0*fR[2])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[20]+fL[20])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[22]+fL[22])+3.0*(fL[7]-1.0*fR[7]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[26]+fL[26])+3.0*(fL[9]-1.0*fR[9]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[33]+fL[33])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[6]+fL[6]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[36]+fL[36])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[8]+fL[8]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[38]+fL[38])+1.732050807568877*(fL[18]-1.0*fR[18])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[19]-1.0*(8.660254037844387*fL[19]+5.0*(fR[11]+fL[11]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[24]-1.0*(8.660254037844387*fL[24]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[29]-1.0*(8.660254037844387*fL[29]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[45]+fL[45])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[17]+fL[17])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[32]-1.0*(8.660254037844387*fL[32]+5.0*(fR[21]+fL[21]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[23]+fL[23]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[35]-1.0*(8.660254037844387*fL[35]+5.0*(fR[25]+fL[25]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[40]-1.0*(8.660254037844387*fL[40]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[28]+fL[28]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[43]-1.0*(8.660254037844387*fL[43]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[44]-1.0*(8.660254037844387*fL[44]+5.0*(fR[37]+fL[37]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[39]+fL[39]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[42]+fL[42]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[12]-6.324555320336761*fL[12]-4.898979485566357*(fR[2]+fL[2])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[20]-94.8683298050514*fL[20]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[1]-42.42640687119286*fL[1])*amax-18.97366596101028*alphaR[1]*fAvg[7]-21.21320343559643*(alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.04166666666666666*((18.97366596101028*fR[22]-18.97366596101028*fL[22]-14.69693845669907*(fR[7]+fL[7])+8.485281374238571*fR[3]-8.485281374238571*fL[3])*amax-4.242640687119286*(alphaR[1]*fAvg[4]+alphaR[0]*fAvg[2])); 
  Ghat[3] = -0.04166666666666666*((18.97366596101028*fR[26]-18.97366596101028*fL[26]-14.69693845669907*(fR[9]+fL[9])+8.485281374238571*fR[4]-8.485281374238571*fL[4])*amax-4.242640687119286*(alphaR[1]*fAvg[5]+alphaR[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[33]-94.86832980505142*fL[33]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[6]-42.42640687119286*fL[6])*amax-18.97366596101028*alphaR[1]*fAvg[11]-21.21320343559643*(alphaR[0]*fAvg[4]+alphaR[1]*fAvg[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[36]-94.86832980505142*fL[36]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[8]-42.42640687119286*fL[8])*amax-18.97366596101028*alphaR[1]*fAvg[13]-21.21320343559643*(alphaR[0]*fAvg[5]+alphaR[1]*fAvg[3])); 
  Ghat[6] = -0.125*((6.324555320336761*fR[38]-6.324555320336761*fL[38]-4.898979485566357*(fR[18]+fL[18])+2.828427124746191*fR[10]-2.828427124746191*fL[10])*amax-1.414213562373095*(alphaR[1]*fAvg[10]+alphaR[0]*fAvg[6])); 
  Ghat[7] = 0.025*((24.49489742783179*(fR[19]+fL[19])-14.14213562373095*fR[11]+14.14213562373095*fL[11])*amax+7.071067811865476*alphaR[0]*fAvg[7]+6.324555320336761*alphaR[1]*fAvg[1]); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[24]+fL[24])-42.42640687119286*fR[13]+42.42640687119286*fL[13])*amax+21.21320343559643*alphaR[1]*fAvg[12]+21.21320343559643*alphaR[0]*fAvg[8]); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[29]+fL[29])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[1]*fAvg[15]+21.21320343559643*alphaR[0]*fAvg[9]); 
  Ghat[10] = -0.008333333333333333*((94.8683298050514*fR[45]-94.8683298050514*fL[45]-73.48469228349535*(fR[31]+fL[31])+42.42640687119286*fR[17]-42.42640687119286*fL[17])*amax-18.97366596101028*alphaR[1]*fAvg[17]-21.21320343559643*(alphaR[0]*fAvg[10]+alphaR[1]*fAvg[6])); 
  Ghat[11] = 0.008333333333333333*((73.48469228349536*(fR[32]+fL[32])-42.42640687119286*fR[21]+42.42640687119286*fL[21])*amax+21.21320343559643*alphaR[0]*fAvg[11]+18.97366596101028*alphaR[1]*fAvg[4]); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[34]+fL[34])-42.42640687119286*fR[23]+42.42640687119286*fL[23])*amax+21.21320343559643*alphaR[0]*fAvg[12]+21.21320343559643*alphaR[1]*fAvg[8]); 
  Ghat[13] = 0.008333333333333333*((73.48469228349536*(fR[35]+fL[35])-42.42640687119286*fR[25]+42.42640687119286*fL[25])*amax+21.21320343559643*alphaR[0]*fAvg[13]+18.97366596101028*alphaR[1]*fAvg[5]); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[40]+fL[40])-42.42640687119286*fR[27]+42.42640687119286*fL[27])*amax+21.21320343559643*alphaR[1]*fAvg[18]+21.21320343559643*alphaR[0]*fAvg[14]); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[28]+42.42640687119286*fL[28])*amax+21.21320343559643*alphaR[0]*fAvg[15]+21.21320343559643*alphaR[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[43]+fL[43])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+21.21320343559643*alphaR[1]*fAvg[19]+21.21320343559643*alphaR[0]*fAvg[16]); 
  Ghat[17] = 0.025*((24.49489742783179*(fR[44]+fL[44])-14.14213562373095*fR[37]+14.14213562373095*fL[37])*amax+7.071067811865476*alphaR[0]*fAvg[17]+6.324555320336761*alphaR[1]*fAvg[10]); 
  Ghat[18] = 0.008333333333333333*((73.48469228349536*(fR[46]+fL[46])-42.42640687119286*fR[39]+42.42640687119286*fL[39])*amax+21.21320343559643*alphaR[0]*fAvg[18]+21.21320343559643*alphaR[1]*fAvg[14]); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[42]+42.42640687119286*fL[42])*amax+21.21320343559643*alphaR[0]*fAvg[19]+21.21320343559643*alphaR[1]*fAvg[16]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = -1.224744871391589*Ghat[0]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = 0.7071067811865475*Ghat[4]; 
  incr[7] = -1.224744871391589*Ghat[2]; 
  incr[8] = 0.7071067811865475*Ghat[5]; 
  incr[9] = -1.224744871391589*Ghat[3]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 0.7071067811865475*Ghat[7]; 
  incr[12] = 1.58113883008419*Ghat[0]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = 0.7071067811865475*Ghat[10]; 
  incr[18] = -1.224744871391589*Ghat[6]; 
  incr[19] = -1.224744871391589*Ghat[7]; 
  incr[20] = 1.58113883008419*Ghat[1]; 
  incr[21] = 0.7071067811865475*Ghat[11]; 
  incr[22] = 1.58113883008419*Ghat[2]; 
  incr[23] = 0.7071067811865475*Ghat[12]; 
  incr[24] = -1.224744871391589*Ghat[8]; 
  incr[25] = 0.7071067811865475*Ghat[13]; 
  incr[26] = 1.58113883008419*Ghat[3]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = 0.7071067811865475*Ghat[15]; 
  incr[29] = -1.224744871391589*Ghat[9]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = -1.224744871391589*Ghat[11]; 
  incr[33] = 1.58113883008419*Ghat[4]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = -1.224744871391589*Ghat[13]; 
  incr[36] = 1.58113883008419*Ghat[5]; 
  incr[37] = 0.7071067811865475*Ghat[17]; 
  incr[38] = 1.58113883008419*Ghat[6]; 
  incr[39] = 0.7071067811865475*Ghat[18]; 
  incr[40] = -1.224744871391589*Ghat[14]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = 0.7071067811865475*Ghat[19]; 
  incr[43] = -1.224744871391589*Ghat[16]; 
  incr[44] = -1.224744871391589*Ghat[17]; 
  incr[45] = 1.58113883008419*Ghat[10]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

  outR[0] += incr[0]*rdy2R; 
  outR[1] += incr[1]*rdy2R; 
  outR[2] += incr[2]*rdy2R; 
  outR[3] += incr[3]*rdy2R; 
  outR[4] += incr[4]*rdy2R; 
  outR[5] += incr[5]*rdy2R; 
  outR[6] += incr[6]*rdy2R; 
  outR[7] += incr[7]*rdy2R; 
  outR[8] += incr[8]*rdy2R; 
  outR[9] += incr[9]*rdy2R; 
  outR[10] += incr[10]*rdy2R; 
  outR[11] += incr[11]*rdy2R; 
  outR[12] += incr[12]*rdy2R; 
  outR[13] += incr[13]*rdy2R; 
  outR[14] += incr[14]*rdy2R; 
  outR[15] += incr[15]*rdy2R; 
  outR[16] += incr[16]*rdy2R; 
  outR[17] += incr[17]*rdy2R; 
  outR[18] += incr[18]*rdy2R; 
  outR[19] += incr[19]*rdy2R; 
  outR[20] += incr[20]*rdy2R; 
  outR[21] += incr[21]*rdy2R; 
  outR[22] += incr[22]*rdy2R; 
  outR[23] += incr[23]*rdy2R; 
  outR[24] += incr[24]*rdy2R; 
  outR[25] += incr[25]*rdy2R; 
  outR[26] += incr[26]*rdy2R; 
  outR[27] += incr[27]*rdy2R; 
  outR[28] += incr[28]*rdy2R; 
  outR[29] += incr[29]*rdy2R; 
  outR[30] += incr[30]*rdy2R; 
  outR[31] += incr[31]*rdy2R; 
  outR[32] += incr[32]*rdy2R; 
  outR[33] += incr[33]*rdy2R; 
  outR[34] += incr[34]*rdy2R; 
  outR[35] += incr[35]*rdy2R; 
  outR[36] += incr[36]*rdy2R; 
  outR[37] += incr[37]*rdy2R; 
  outR[38] += incr[38]*rdy2R; 
  outR[39] += incr[39]*rdy2R; 
  outR[40] += incr[40]*rdy2R; 
  outR[41] += incr[41]*rdy2R; 
  outR[42] += incr[42]*rdy2R; 
  outR[43] += incr[43]*rdy2R; 
  outR[44] += incr[44]*rdy2R; 
  outR[45] += incr[45]*rdy2R; 
  outR[46] += incr[46]*rdy2R; 
  outR[47] += incr[47]*rdy2R; 

  outL[0] += -1.0*incr[0]*rdy2L; 
  outL[1] += -1.0*incr[1]*rdy2L; 
  outL[2] += incr[2]*rdy2L; 
  outL[3] += -1.0*incr[3]*rdy2L; 
  outL[4] += -1.0*incr[4]*rdy2L; 
  outL[5] += incr[5]*rdy2L; 
  outL[6] += -1.0*incr[6]*rdy2L; 
  outL[7] += incr[7]*rdy2L; 
  outL[8] += -1.0*incr[8]*rdy2L; 
  outL[9] += incr[9]*rdy2L; 
  outL[10] += -1.0*incr[10]*rdy2L; 
  outL[11] += -1.0*incr[11]*rdy2L; 
  outL[12] += -1.0*incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += -1.0*incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 
  outL[16] += incr[16]*rdy2L; 
  outL[17] += -1.0*incr[17]*rdy2L; 
  outL[18] += incr[18]*rdy2L; 
  outL[19] += incr[19]*rdy2L; 
  outL[20] += -1.0*incr[20]*rdy2L; 
  outL[21] += -1.0*incr[21]*rdy2L; 
  outL[22] += -1.0*incr[22]*rdy2L; 
  outL[23] += -1.0*incr[23]*rdy2L; 
  outL[24] += incr[24]*rdy2L; 
  outL[25] += -1.0*incr[25]*rdy2L; 
  outL[26] += -1.0*incr[26]*rdy2L; 
  outL[27] += -1.0*incr[27]*rdy2L; 
  outL[28] += -1.0*incr[28]*rdy2L; 
  outL[29] += incr[29]*rdy2L; 
  outL[30] += -1.0*incr[30]*rdy2L; 
  outL[31] += incr[31]*rdy2L; 
  outL[32] += incr[32]*rdy2L; 
  outL[33] += -1.0*incr[33]*rdy2L; 
  outL[34] += incr[34]*rdy2L; 
  outL[35] += incr[35]*rdy2L; 
  outL[36] += -1.0*incr[36]*rdy2L; 
  outL[37] += -1.0*incr[37]*rdy2L; 
  outL[38] += -1.0*incr[38]*rdy2L; 
  outL[39] += -1.0*incr[39]*rdy2L; 
  outL[40] += incr[40]*rdy2L; 
  outL[41] += incr[41]*rdy2L; 
  outL[42] += -1.0*incr[42]*rdy2L; 
  outL[43] += incr[43]*rdy2L; 
  outL[44] += incr[44]*rdy2L; 
  outL[45] += -1.0*incr[45]*rdy2L; 
  outL[46] += incr[46]*rdy2L; 
  outL[47] += incr[47]*rdy2L; 
return std::abs(alphaSurfAvgR); 
} 
double GyrokineticGenGeoSurf2x2vSer_vpar_P2_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[11] = 2.0*phi[4]*q_; 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 

  double BstarXdBmagR[48]; 

  double BstarYdBmagR[48]; 

  double alphaR[20]; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = 0.0; 

  double incr[48]; 
  // alpha == 0, so nothing to do.
  return std::abs(alphaSurfAvgR);
}
double GyrokineticGenGeoSurf2x2vSer_x_P2_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[8] = (1.154700538379252*bmag[1])/rdmu2R; 
  hamilR[11] = 2.0*(bmag[4]*wmuR+phi[4]*q_); 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 
  hamilR[25] = (1.154700538379251*bmag[4])/rdmu2R; 

  double BstarXdBmagR[48]; 

  double alphaR[20]; 
  alphaR[0] = -(0.125*(((27.38612787525831*b_z[4]-21.21320343559643*b_z[1]+12.24744871391589*b_z[0])*jacobTotInv[4]+(12.24744871391589*jacobTotInv[0]-21.21320343559643*jacobTotInv[1])*b_z[4]+(16.43167672515499*b_z[1]-9.48683298050514*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(5.477225575051662*b_z[0]-9.48683298050514*b_z[1]))*hamilR[19]+(((-21.21320343559643*b_z[4])+16.43167672515498*b_z[1]-9.48683298050514*b_z[0])*jacobTotInv[4]+(16.43167672515498*jacobTotInv[1]-9.48683298050514*jacobTotInv[0])*b_z[4]+(7.348469228349534*b_z[0]-12.72792206135786*b_z[1])*jacobTotInv[1]+jacobTotInv[0]*(7.348469228349534*b_z[1]-4.242640687119286*b_z[0]))*hamilR[5]+hamilR[2]*((12.24744871391589*b_z[4]-9.48683298050514*b_z[1]+5.477225575051662*b_z[0])*jacobTotInv[4]+(5.477225575051662*jacobTotInv[0]-9.48683298050514*jacobTotInv[1])*b_z[4]+(7.348469228349534*b_z[1]-4.242640687119286*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(2.449489742783178*b_z[0]-4.242640687119286*b_z[1])))*rdy2R)/q_; 
  alphaR[1] = (0.125*(((47.43416490252569*b_z[4]-36.74234614174768*b_z[1]+21.21320343559643*b_z[0])*jacobTotInv[4]+(21.21320343559643*jacobTotInv[0]-36.74234614174768*jacobTotInv[1])*b_z[4]+(28.46049894151541*b_z[1]-16.43167672515499*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(9.48683298050514*b_z[0]-16.43167672515499*b_z[1]))*hamilR[20]+(((-27.38612787525831*b_z[4])+21.21320343559643*b_z[1]-12.24744871391589*b_z[0])*jacobTotInv[4]+(21.21320343559643*jacobTotInv[1]-12.24744871391589*jacobTotInv[0])*b_z[4]+(9.48683298050514*b_z[0]-16.43167672515498*b_z[1])*jacobTotInv[1]+jacobTotInv[0]*(9.48683298050514*b_z[1]-5.477225575051662*b_z[0]))*hamilR[12])*rdy2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(((19.36491673103708*b_z[4]-15.0*b_z[1]+8.660254037844387*b_z[0])*jacobTotInv[4]+(8.660254037844387*jacobTotInv[0]-15.0*jacobTotInv[1])*b_z[4]+(11.61895003862225*b_z[1]-6.708203932499369*b_z[0])*jacobTotInv[1]-6.708203932499369*jacobTotInv[0]*b_z[1]+3.872983346207417*b_z[0]*jacobTotInv[0])*hamilR[19]+(((-15.0*b_z[4])+11.61895003862225*b_z[1]-6.708203932499369*b_z[0])*jacobTotInv[4]+(11.61895003862225*jacobTotInv[1]-6.708203932499369*jacobTotInv[0])*b_z[4]+(5.196152422706631*b_z[0]-9.0*b_z[1])*jacobTotInv[1]+5.196152422706631*jacobTotInv[0]*b_z[1]-3.0*b_z[0]*jacobTotInv[0])*hamilR[5]+(8.660254037844386*hamilR[2]*b_z[4]+(3.872983346207417*b_z[0]-6.708203932499369*b_z[1])*hamilR[2])*jacobTotInv[4]+(3.872983346207417*jacobTotInv[0]-6.708203932499369*jacobTotInv[1])*hamilR[2]*b_z[4]+((5.196152422706631*b_z[1]-3.0*b_z[0])*jacobTotInv[1]-3.0*jacobTotInv[0]*b_z[1]+1.732050807568877*b_z[0]*jacobTotInv[0])*hamilR[2])*rdy2R)/q_; 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[11]+fL[11])+1.732050807568877*(fL[1]-1.0*fR[1])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[19]+fL[19])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[2]+fL[2])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[21]+fL[21])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[25]+fL[25])+3.0*(fL[8]-1.0*fR[8]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[32]+fL[32])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[7]+fL[7]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[35]+fL[35])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[9]+fL[9]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[37]+fL[37])+1.732050807568877*(fL[17]-1.0*fR[17])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[20]-1.0*(8.660254037844387*fL[20]+5.0*(fR[12]+fL[12]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[23]-1.0*(8.660254037844387*fL[23]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[28]-1.0*(8.660254037844387*fL[28]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[44]+fL[44])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[18]+fL[18])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[33]-1.0*(8.660254037844387*fL[33]+5.0*(fR[22]+fL[22]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[24]+fL[24]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[36]-1.0*(8.660254037844387*fL[36]+5.0*(fR[26]+fL[26]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[39]-1.0*(8.660254037844387*fL[39]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[29]+fL[29]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[42]-1.0*(8.660254037844387*fL[42]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[45]-1.0*(8.660254037844387*fL[45]+5.0*(fR[38]+fL[38]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[40]+fL[40]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[43]+fL[43]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[11]-6.324555320336761*fL[11]-4.898979485566357*(fR[1]+fL[1])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[19]-94.8683298050514*fL[19]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[2]-42.42640687119286*fL[2])*amax-18.97366596101028*alphaR[1]*fAvg[7]-21.21320343559643*(alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.04166666666666666*((18.97366596101028*fR[21]-18.97366596101028*fL[21]-14.69693845669907*(fR[6]+fL[6])+8.485281374238571*fR[3]-8.485281374238571*fL[3])*amax-4.242640687119286*(alphaR[1]*fAvg[4]+alphaR[0]*fAvg[2])); 
  Ghat[3] = -0.04166666666666666*((18.97366596101028*fR[25]-18.97366596101028*fL[25]-14.69693845669907*(fR[8]+fL[8])+8.485281374238571*fR[4]-8.485281374238571*fL[4])*amax-4.242640687119286*(alphaR[1]*fAvg[5]+alphaR[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[32]-94.86832980505142*fL[32]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[7]-42.42640687119286*fL[7])*amax-18.97366596101028*alphaR[1]*fAvg[11]-21.21320343559643*(alphaR[0]*fAvg[4]+alphaR[1]*fAvg[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[35]-94.86832980505142*fL[35]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[9]-42.42640687119286*fL[9])*amax-18.97366596101028*alphaR[1]*fAvg[13]-21.21320343559643*(alphaR[0]*fAvg[5]+alphaR[1]*fAvg[3])); 
  Ghat[6] = -0.125*((6.324555320336761*fR[37]-6.324555320336761*fL[37]-4.898979485566357*(fR[17]+fL[17])+2.828427124746191*fR[10]-2.828427124746191*fL[10])*amax-1.414213562373095*(alphaR[1]*fAvg[10]+alphaR[0]*fAvg[6])); 
  Ghat[7] = 0.025*((24.49489742783179*(fR[20]+fL[20])-14.14213562373095*fR[12]+14.14213562373095*fL[12])*amax+7.071067811865476*alphaR[0]*fAvg[7]+6.324555320336761*alphaR[1]*fAvg[1]); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[23]+fL[23])-42.42640687119286*fR[13]+42.42640687119286*fL[13])*amax+21.21320343559643*alphaR[1]*fAvg[12]+21.21320343559643*alphaR[0]*fAvg[8]); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[28]+fL[28])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[1]*fAvg[15]+21.21320343559643*alphaR[0]*fAvg[9]); 
  Ghat[10] = -0.008333333333333333*((94.8683298050514*fR[44]-94.8683298050514*fL[44]-73.48469228349535*(fR[31]+fL[31])+42.42640687119286*fR[18]-42.42640687119286*fL[18])*amax-18.97366596101028*alphaR[1]*fAvg[17]-21.21320343559643*(alphaR[0]*fAvg[10]+alphaR[1]*fAvg[6])); 
  Ghat[11] = 0.008333333333333333*((73.48469228349536*(fR[33]+fL[33])-42.42640687119286*fR[22]+42.42640687119286*fL[22])*amax+21.21320343559643*alphaR[0]*fAvg[11]+18.97366596101028*alphaR[1]*fAvg[4]); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[34]+fL[34])-42.42640687119286*fR[24]+42.42640687119286*fL[24])*amax+21.21320343559643*alphaR[0]*fAvg[12]+21.21320343559643*alphaR[1]*fAvg[8]); 
  Ghat[13] = 0.008333333333333333*((73.48469228349536*(fR[36]+fL[36])-42.42640687119286*fR[26]+42.42640687119286*fL[26])*amax+21.21320343559643*alphaR[0]*fAvg[13]+18.97366596101028*alphaR[1]*fAvg[5]); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[39]+fL[39])-42.42640687119286*fR[27]+42.42640687119286*fL[27])*amax+21.21320343559643*alphaR[1]*fAvg[18]+21.21320343559643*alphaR[0]*fAvg[14]); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[29]+42.42640687119286*fL[29])*amax+21.21320343559643*alphaR[0]*fAvg[15]+21.21320343559643*alphaR[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[42]+fL[42])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+21.21320343559643*alphaR[1]*fAvg[19]+21.21320343559643*alphaR[0]*fAvg[16]); 
  Ghat[17] = 0.025*((24.49489742783179*(fR[45]+fL[45])-14.14213562373095*fR[38]+14.14213562373095*fL[38])*amax+7.071067811865476*alphaR[0]*fAvg[17]+6.324555320336761*alphaR[1]*fAvg[10]); 
  Ghat[18] = 0.008333333333333333*((73.48469228349536*(fR[46]+fL[46])-42.42640687119286*fR[40]+42.42640687119286*fL[40])*amax+21.21320343559643*alphaR[0]*fAvg[18]+21.21320343559643*alphaR[1]*fAvg[14]); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[43]+42.42640687119286*fL[43])*amax+21.21320343559643*alphaR[0]*fAvg[19]+21.21320343559643*alphaR[1]*fAvg[16]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = -1.224744871391589*Ghat[0]; 
  incr[2] = 0.7071067811865475*Ghat[1]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = -1.224744871391589*Ghat[2]; 
  incr[7] = 0.7071067811865475*Ghat[4]; 
  incr[8] = -1.224744871391589*Ghat[3]; 
  incr[9] = 0.7071067811865475*Ghat[5]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 1.58113883008419*Ghat[0]; 
  incr[12] = 0.7071067811865475*Ghat[7]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = -1.224744871391589*Ghat[6]; 
  incr[18] = 0.7071067811865475*Ghat[10]; 
  incr[19] = 1.58113883008419*Ghat[1]; 
  incr[20] = -1.224744871391589*Ghat[7]; 
  incr[21] = 1.58113883008419*Ghat[2]; 
  incr[22] = 0.7071067811865475*Ghat[11]; 
  incr[23] = -1.224744871391589*Ghat[8]; 
  incr[24] = 0.7071067811865475*Ghat[12]; 
  incr[25] = 1.58113883008419*Ghat[3]; 
  incr[26] = 0.7071067811865475*Ghat[13]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = -1.224744871391589*Ghat[9]; 
  incr[29] = 0.7071067811865475*Ghat[15]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = 1.58113883008419*Ghat[4]; 
  incr[33] = -1.224744871391589*Ghat[11]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = 1.58113883008419*Ghat[5]; 
  incr[36] = -1.224744871391589*Ghat[13]; 
  incr[37] = 1.58113883008419*Ghat[6]; 
  incr[38] = 0.7071067811865475*Ghat[17]; 
  incr[39] = -1.224744871391589*Ghat[14]; 
  incr[40] = 0.7071067811865475*Ghat[18]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = -1.224744871391589*Ghat[16]; 
  incr[43] = 0.7071067811865475*Ghat[19]; 
  incr[44] = 1.58113883008419*Ghat[10]; 
  incr[45] = -1.224744871391589*Ghat[17]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 
  outR[8] += incr[8]*rdx2R; 
  outR[9] += incr[9]*rdx2R; 
  outR[10] += incr[10]*rdx2R; 
  outR[11] += incr[11]*rdx2R; 
  outR[12] += incr[12]*rdx2R; 
  outR[13] += incr[13]*rdx2R; 
  outR[14] += incr[14]*rdx2R; 
  outR[15] += incr[15]*rdx2R; 
  outR[16] += incr[16]*rdx2R; 
  outR[17] += incr[17]*rdx2R; 
  outR[18] += incr[18]*rdx2R; 
  outR[19] += incr[19]*rdx2R; 
  outR[20] += incr[20]*rdx2R; 
  outR[21] += incr[21]*rdx2R; 
  outR[22] += incr[22]*rdx2R; 
  outR[23] += incr[23]*rdx2R; 
  outR[24] += incr[24]*rdx2R; 
  outR[25] += incr[25]*rdx2R; 
  outR[26] += incr[26]*rdx2R; 
  outR[27] += incr[27]*rdx2R; 
  outR[28] += incr[28]*rdx2R; 
  outR[29] += incr[29]*rdx2R; 
  outR[30] += incr[30]*rdx2R; 
  outR[31] += incr[31]*rdx2R; 
  outR[32] += incr[32]*rdx2R; 
  outR[33] += incr[33]*rdx2R; 
  outR[34] += incr[34]*rdx2R; 
  outR[35] += incr[35]*rdx2R; 
  outR[36] += incr[36]*rdx2R; 
  outR[37] += incr[37]*rdx2R; 
  outR[38] += incr[38]*rdx2R; 
  outR[39] += incr[39]*rdx2R; 
  outR[40] += incr[40]*rdx2R; 
  outR[41] += incr[41]*rdx2R; 
  outR[42] += incr[42]*rdx2R; 
  outR[43] += incr[43]*rdx2R; 
  outR[44] += incr[44]*rdx2R; 
  outR[45] += incr[45]*rdx2R; 
  outR[46] += incr[46]*rdx2R; 
  outR[47] += incr[47]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += -1.0*incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += incr[6]*rdx2L; 
  outL[7] += -1.0*incr[7]*rdx2L; 
  outL[8] += incr[8]*rdx2L; 
  outL[9] += -1.0*incr[9]*rdx2L; 
  outL[10] += -1.0*incr[10]*rdx2L; 
  outL[11] += -1.0*incr[11]*rdx2L; 
  outL[12] += -1.0*incr[12]*rdx2L; 
  outL[13] += -1.0*incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 
  outL[16] += incr[16]*rdx2L; 
  outL[17] += incr[17]*rdx2L; 
  outL[18] += -1.0*incr[18]*rdx2L; 
  outL[19] += -1.0*incr[19]*rdx2L; 
  outL[20] += incr[20]*rdx2L; 
  outL[21] += -1.0*incr[21]*rdx2L; 
  outL[22] += -1.0*incr[22]*rdx2L; 
  outL[23] += incr[23]*rdx2L; 
  outL[24] += -1.0*incr[24]*rdx2L; 
  outL[25] += -1.0*incr[25]*rdx2L; 
  outL[26] += -1.0*incr[26]*rdx2L; 
  outL[27] += -1.0*incr[27]*rdx2L; 
  outL[28] += incr[28]*rdx2L; 
  outL[29] += -1.0*incr[29]*rdx2L; 
  outL[30] += -1.0*incr[30]*rdx2L; 
  outL[31] += incr[31]*rdx2L; 
  outL[32] += -1.0*incr[32]*rdx2L; 
  outL[33] += incr[33]*rdx2L; 
  outL[34] += incr[34]*rdx2L; 
  outL[35] += -1.0*incr[35]*rdx2L; 
  outL[36] += incr[36]*rdx2L; 
  outL[37] += -1.0*incr[37]*rdx2L; 
  outL[38] += -1.0*incr[38]*rdx2L; 
  outL[39] += incr[39]*rdx2L; 
  outL[40] += -1.0*incr[40]*rdx2L; 
  outL[41] += incr[41]*rdx2L; 
  outL[42] += incr[42]*rdx2L; 
  outL[43] += -1.0*incr[43]*rdx2L; 
  outL[44] += -1.0*incr[44]*rdx2L; 
  outL[45] += incr[45]*rdx2L; 
  outL[46] += incr[46]*rdx2L; 
  outL[47] += incr[47]*rdx2L; 
return std::abs(alphaSurfAvgR); 
} 
double GyrokineticGenGeoSurf2x2vSer_y_P2_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[8] = (1.154700538379252*bmag[1])/rdmu2R; 
  hamilR[11] = 2.0*(bmag[4]*wmuR+phi[4]*q_); 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 
  hamilR[25] = (1.154700538379251*bmag[4])/rdmu2R; 

  double BstarYdBmagR[48]; 
  BstarYdBmagR[0] = -(1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[1] = -(1.732050807568877*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[3] = -(1.0*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[6] = -(1.0*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[11] = -(1.732050807568877*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[21] = -(1.0*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[20]; 
  alphaR[0] = (0.03535533905932736*((19.36491673103708*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamilR[20]+((-30.0*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-33.54101966249684*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[19]+(17.32050807568877*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+19.36491673103709*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[11]+(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*(8.660254037844386*hamilR[1]-15.0*hamilR[5]))*m_*rdx2R+(19.36491673103709*BstarYdBmagR[3]*hamilR[13]+8.660254037844386*BstarYdBmagR[0]*hamilR[3])*q_*rdvpar2R))/(m_*q_); 
  alphaR[1] = (0.005050762722761052*(((121.2435565298214*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+135.5544171172596*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[20]+(((-368.9512162874653*b_z[4])-210.0*b_z[0])*jacobTotInv[4]-210.0*jacobTotInv[0]*b_z[4]-422.6168477474601*b_z[1]*jacobTotInv[1]-234.7871376374779*b_z[0]*jacobTotInv[0])*hamilR[19]+((213.014084041408*b_z[4]+121.2435565298214*b_z[0])*jacobTotInv[4]+121.2435565298214*jacobTotInv[0]*b_z[4]+243.9979508110672*b_z[1]*jacobTotInv[1]+135.5544171172596*b_z[0]*jacobTotInv[0])*hamilR[11]+((-93.91485505499116*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-105.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[5]+hamilR[1]*(54.22176684690384*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+60.6217782649107*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])))*m_*rdx2R+(135.5544171172596*BstarYdBmagR[6]*hamilR[13]+60.6217782649107*BstarYdBmagR[1]*hamilR[3])*q_*rdvpar2R))/(m_*q_); 
  alphaR[2] = (0.1767766952966368*(3.872983346207417*BstarYdBmagR[0]*hamilR[13]+1.732050807568877*BstarYdBmagR[3]*hamilR[3])*rdvpar2R)/m_; 
  alphaR[3] = (0.03535533905932736*((17.32050807568877*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+19.36491673103708*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[25]+8.660254037844386*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamilR[8])*rdx2R)/q_; 
  alphaR[4] = (0.1767766952966368*(3.872983346207417*BstarYdBmagR[1]*hamilR[13]+1.732050807568877*hamilR[3]*BstarYdBmagR[6])*rdvpar2R)/m_; 
  alphaR[5] = (0.005050762722761052*(((213.0140840414079*b_z[4]+121.2435565298214*b_z[0])*jacobTotInv[4]+121.2435565298214*jacobTotInv[0]*b_z[4]+243.9979508110673*b_z[1]*jacobTotInv[1]+135.5544171172596*b_z[0]*jacobTotInv[0])*hamilR[25]+(54.22176684690384*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+60.6217782649107*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[8])*rdx2R)/q_; 
  alphaR[7] = (0.005050762722761052*((((86.60254037844389*b_z[4]+135.5544171172596*b_z[0])*jacobTotInv[4]+135.5544171172596*jacobTotInv[0]*b_z[4]+121.2435565298214*b_z[1]*jacobTotInv[1])*hamilR[20]+((-368.9512162874653*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-210.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[19]+(213.014084041408*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+121.2435565298214*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[11]+(((-67.0820393249937*b_z[4])-105.0*b_z[0])*jacobTotInv[4]-105.0*jacobTotInv[0]*b_z[4]-93.91485505499116*b_z[1]*jacobTotInv[1])*hamilR[5]+hamilR[1]*((38.72983346207418*b_z[4]+60.6217782649107*b_z[0])*jacobTotInv[4]+60.6217782649107*jacobTotInv[0]*b_z[4]+54.22176684690384*b_z[1]*jacobTotInv[1]))*m_*rdx2R+(135.5544171172596*hamilR[13]*BstarYdBmagR[21]+60.6217782649107*hamilR[3]*BstarYdBmagR[11])*q_*rdvpar2R))/(m_*q_); 
  alphaR[8] = (0.6123724356957944*BstarYdBmagR[3]*hamilR[13]*rdvpar2R)/m_; 
  alphaR[11] = (0.04564354645876382*(6.708203932499369*hamilR[3]*BstarYdBmagR[21]+15.0*BstarYdBmagR[11]*hamilR[13])*rdvpar2R)/m_; 
  alphaR[12] = (0.6123724356957944*BstarYdBmagR[6]*hamilR[13]*rdvpar2R)/m_; 
  alphaR[13] = (0.005050762722761052*((213.014084041408*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+121.2435565298214*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[25]+((38.72983346207417*b_z[4]+60.62177826491071*b_z[0])*jacobTotInv[4]+60.62177826491071*jacobTotInv[0]*b_z[4]+54.22176684690384*b_z[1]*jacobTotInv[1])*hamilR[8])*rdx2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.00625*(((19.36491673103708*b_z[4]*jacobTotInv[4]+19.36491673103708*b_z[1]*jacobTotInv[1]+19.36491673103708*b_z[0]*jacobTotInv[0])*hamilR[20]+((-30.0*b_z[1]*jacobTotInv[4])-30.0*jacobTotInv[1]*b_z[4]-33.54101966249684*b_z[0]*jacobTotInv[1]-33.54101966249684*jacobTotInv[0]*b_z[1])*hamilR[19]+(17.32050807568877*b_z[1]*jacobTotInv[4]+17.32050807568877*jacobTotInv[1]*b_z[4]+19.36491673103709*b_z[0]*jacobTotInv[1]+19.36491673103709*jacobTotInv[0]*b_z[1])*hamilR[11]+((-15.0*b_z[4]*jacobTotInv[4])-15.0*b_z[1]*jacobTotInv[1]-15.0*b_z[0]*jacobTotInv[0])*hamilR[5]+8.660254037844386*hamilR[1]*b_z[4]*jacobTotInv[4]+8.660254037844386*b_z[1]*hamilR[1]*jacobTotInv[1]+8.660254037844386*b_z[0]*jacobTotInv[0]*hamilR[1])*m_*rdx2R+(19.36491673103709*BstarYdBmagR[3]*hamilR[13]+8.660254037844386*BstarYdBmagR[0]*hamilR[3])*q_*rdvpar2R))/(m_*q_); 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[12]+fL[12])+1.732050807568877*(fL[2]-1.0*fR[2])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[20]+fL[20])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[22]+fL[22])+3.0*(fL[7]-1.0*fR[7]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[26]+fL[26])+3.0*(fL[9]-1.0*fR[9]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[33]+fL[33])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[6]+fL[6]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[36]+fL[36])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[8]+fL[8]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[38]+fL[38])+1.732050807568877*(fL[18]-1.0*fR[18])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[19]-1.0*(8.660254037844387*fL[19]+5.0*(fR[11]+fL[11]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[24]-1.0*(8.660254037844387*fL[24]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[29]-1.0*(8.660254037844387*fL[29]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[45]+fL[45])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[17]+fL[17])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[32]-1.0*(8.660254037844387*fL[32]+5.0*(fR[21]+fL[21]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[23]+fL[23]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[35]-1.0*(8.660254037844387*fL[35]+5.0*(fR[25]+fL[25]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[40]-1.0*(8.660254037844387*fL[40]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[28]+fL[28]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[43]-1.0*(8.660254037844387*fL[43]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[44]-1.0*(8.660254037844387*fL[44]+5.0*(fR[37]+fL[37]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[39]+fL[39]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[42]+fL[42]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[12]-6.324555320336761*fL[12]-4.898979485566357*(fR[2]+fL[2])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[13]*fAvg[13]+alphaR[12]*fAvg[12]+alphaR[11]*fAvg[11]+alphaR[8]*fAvg[8]+alphaR[7]*fAvg[7]+alphaR[5]*fAvg[5]+alphaR[4]*fAvg[4]+alphaR[3]*fAvg[3]+alphaR[2]*fAvg[2]+alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[20]-94.8683298050514*fL[20]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[1]-42.42640687119286*fL[1])*amax-18.97366596101028*(alphaR[5]*fAvg[13]+fAvg[5]*alphaR[13])-21.21320343559643*(alphaR[8]*fAvg[12]+fAvg[8]*alphaR[12])-18.97366596101028*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11]+alphaR[1]*fAvg[7]+fAvg[1]*alphaR[7])-21.21320343559643*(alphaR[3]*fAvg[5]+fAvg[3]*alphaR[5]+alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4]+alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[22]-94.8683298050514*fL[22]-73.48469228349535*(fR[7]+fL[7])+42.42640687119286*fR[3]-42.42640687119286*fL[3])*amax-21.21320343559643*alphaR[13]*fAvg[17]-18.97366596101028*(alphaR[4]*fAvg[12]+fAvg[4]*alphaR[12])-21.21320343559643*(alphaR[7]*fAvg[11]+fAvg[7]*alphaR[11])-21.21320343559643*alphaR[5]*fAvg[10]-18.97366596101028*(alphaR[2]*fAvg[8]+fAvg[2]*alphaR[8])-21.21320343559643*(alphaR[3]*fAvg[6]+alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4]+alphaR[0]*fAvg[2]+fAvg[0]*alphaR[2])); 
  Ghat[3] = -0.008333333333333333*((94.8683298050514*fR[26]-94.8683298050514*fL[26]-73.48469228349535*(fR[9]+fL[9])+42.42640687119286*fR[4]-42.42640687119286*fL[4])*amax-21.21320343559643*(alphaR[12]*fAvg[18]+alphaR[11]*fAvg[17])-18.97366596101028*alphaR[5]*fAvg[15]-21.21320343559643*(alphaR[8]*fAvg[14]+alphaR[7]*fAvg[13]+fAvg[7]*alphaR[13])-21.21320343559643*alphaR[4]*fAvg[10]-18.97366596101028*alphaR[3]*fAvg[9]-21.21320343559643*(alphaR[2]*fAvg[6]+alphaR[1]*fAvg[5]+fAvg[1]*alphaR[5]+alphaR[0]*fAvg[3]+fAvg[0]*alphaR[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[33]-94.86832980505142*fL[33]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[6]-42.42640687119286*fL[6])*amax-18.97366596101028*(alphaR[5]*fAvg[17]+fAvg[10]*alphaR[13])+((-16.97056274847715*alphaR[11])-18.97366596101028*alphaR[2])*fAvg[12]+((-16.97056274847715*fAvg[11])-18.97366596101028*fAvg[2])*alphaR[12]-18.97366596101028*(alphaR[1]*fAvg[11]+fAvg[1]*alphaR[11])-21.21320343559643*alphaR[3]*fAvg[10]-18.97366596101028*(alphaR[4]*fAvg[8]+fAvg[4]*alphaR[8]+alphaR[4]*fAvg[7]+fAvg[4]*alphaR[7])-21.21320343559643*(alphaR[5]*fAvg[6]+alphaR[0]*fAvg[4]+fAvg[0]*alphaR[4]+alphaR[1]*fAvg[2]+fAvg[1]*alphaR[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[36]-94.86832980505142*fL[36]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[8]-42.42640687119286*fL[8])*amax-21.21320343559643*alphaR[8]*fAvg[18]-18.97366596101028*alphaR[4]*fAvg[17]+((-16.97056274847715*alphaR[13])-18.97366596101028*alphaR[3])*fAvg[15]-21.21320343559643*alphaR[12]*fAvg[14]-18.97366596101028*(alphaR[1]*fAvg[13]+fAvg[1]*alphaR[13])+fAvg[10]*((-18.97366596101028*alphaR[11])-21.21320343559643*alphaR[2])-18.97366596101028*(alphaR[5]*(fAvg[9]+fAvg[7])+fAvg[5]*alphaR[7])-21.21320343559643*(alphaR[4]*fAvg[6]+alphaR[0]*fAvg[5]+fAvg[0]*alphaR[5]+alphaR[1]*fAvg[3]+fAvg[1]*alphaR[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[38]-94.86832980505142*fL[38]-73.48469228349535*(fR[18]+fL[18])+42.42640687119286*fR[10]-42.42640687119286*fL[10])*amax-18.97366596101028*(alphaR[5]*fAvg[19]+alphaR[4]*fAvg[18])-21.21320343559643*alphaR[7]*fAvg[17]-18.97366596101028*(alphaR[3]*fAvg[16]+alphaR[2]*fAvg[14])-21.21320343559643*(alphaR[11]*fAvg[13]+fAvg[11]*alphaR[13])+fAvg[10]*((-18.97366596101028*alphaR[12])-21.21320343559643*alphaR[1])-18.97366596101028*fAvg[6]*alphaR[8]-21.21320343559643*(alphaR[0]*fAvg[6]+alphaR[4]*fAvg[5]+fAvg[4]*alphaR[5]+alphaR[2]*fAvg[3]+fAvg[2]*alphaR[3])); 
  Ghat[7] = 0.001190476190476191*((514.3928459844675*(fR[19]+fL[19])-296.98484809835*fR[11]+296.98484809835*fL[11])*amax+(94.86832980505142*alphaR[13]+148.492424049175*alphaR[3])*fAvg[13]+148.492424049175*fAvg[3]*alphaR[13]+132.815661727072*alphaR[12]*fAvg[12]+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[11]+148.492424049175*fAvg[2]*alphaR[11]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[7]+148.492424049175*fAvg[0]*alphaR[7]+132.815661727072*(alphaR[5]*fAvg[5]+alphaR[4]*fAvg[4]+alphaR[1]*fAvg[1])); 
  Ghat[8] = 0.001190476190476191*((514.3928459844675*(fR[24]+fL[24])-296.98484809835*fR[13]+296.98484809835*fL[13])*amax+148.492424049175*(alphaR[5]*fAvg[18]+alphaR[3]*fAvg[14])+(94.86832980505142*alphaR[12]+148.492424049175*alphaR[1])*fAvg[12]+148.492424049175*fAvg[1]*alphaR[12]+132.815661727072*alphaR[11]*fAvg[11]+(94.86832980505142*alphaR[8]+148.492424049175*alphaR[0])*fAvg[8]+148.492424049175*fAvg[0]*alphaR[8]+132.815661727072*(alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[29]+fL[29])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[4]*fAvg[19]+21.21320343559643*(alphaR[2]*fAvg[16]+alphaR[1]*fAvg[15])+18.97366596101028*alphaR[13]*fAvg[13]+21.21320343559643*alphaR[0]*fAvg[9]+18.97366596101028*(alphaR[5]*fAvg[5]+alphaR[3]*fAvg[3])); 
  Ghat[10] = -0.001666666666666667*((474.341649025257*fR[45]-474.341649025257*fL[45]-367.4234614174767*(fR[31]+fL[31])+212.1320343559643*fR[17]-212.1320343559643*fL[17])*amax+((-84.85281374238573*alphaR[13])-94.86832980505142*alphaR[3])*fAvg[19]+((-84.85281374238573*alphaR[11])-94.86832980505142*alphaR[2])*fAvg[18]+((-84.85281374238573*alphaR[12])-94.86832980505142*alphaR[1])*fAvg[17]-94.8683298050514*(alphaR[5]*fAvg[16]+alphaR[4]*(fAvg[14]+fAvg[13])+fAvg[4]*alphaR[13]+fAvg[6]*alphaR[12]+alphaR[5]*fAvg[11]+fAvg[5]*alphaR[11])+((-94.86832980505142*(alphaR[8]+alphaR[7]))-106.0660171779821*alphaR[0])*fAvg[10]-106.0660171779821*(alphaR[1]*fAvg[6]+alphaR[2]*fAvg[5]+fAvg[2]*alphaR[5]+alphaR[3]*fAvg[4]+fAvg[3]*alphaR[4])); 
  Ghat[11] = 0.001190476190476191*((514.3928459844675*(fR[32]+fL[32])-296.98484809835*fR[21]+296.98484809835*fL[21])*amax+(94.86832980505142*alphaR[13]+148.492424049175*alphaR[3])*fAvg[17]+148.492424049175*fAvg[6]*alphaR[13]+118.79393923934*(alphaR[4]*fAvg[12]+fAvg[4]*alphaR[12])+(132.815661727072*alphaR[8]+94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[11]+(132.815661727072*fAvg[8]+94.86832980505142*fAvg[7]+148.492424049175*fAvg[0])*alphaR[11]+132.815661727072*alphaR[5]*fAvg[10]+148.492424049175*(alphaR[2]*fAvg[7]+fAvg[2]*alphaR[7])+132.815661727072*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4])); 
  Ghat[12] = 0.001190476190476191*((514.3928459844675*(fR[34]+fL[34])-296.98484809835*fR[23]+296.98484809835*fL[23])*amax+(132.815661727072*alphaR[13]+148.492424049175*alphaR[3])*fAvg[18]+148.492424049175*alphaR[5]*fAvg[14]+(94.86832980505142*alphaR[8]+132.815661727072*alphaR[7]+148.492424049175*alphaR[0])*fAvg[12]+(94.86832980505142*fAvg[8]+132.815661727072*fAvg[7]+148.492424049175*fAvg[0])*alphaR[12]+118.79393923934*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11])+148.492424049175*(alphaR[1]*fAvg[8]+fAvg[1]*alphaR[8])+132.815661727072*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4])); 
  Ghat[13] = 0.001190476190476191*((514.3928459844675*(fR[35]+fL[35])-296.98484809835*fR[25]+296.98484809835*fL[25])*amax+132.815661727072*alphaR[12]*fAvg[18]+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[17]+118.79393923934*alphaR[5]*fAvg[15]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[13]+(132.815661727072*fAvg[9]+94.86832980505142*fAvg[7]+148.492424049175*fAvg[0])*alphaR[13]+148.492424049175*fAvg[6]*alphaR[11]+132.815661727072*alphaR[4]*fAvg[10]+148.492424049175*(alphaR[3]*fAvg[7]+fAvg[3]*alphaR[7])+132.815661727072*(alphaR[1]*fAvg[5]+fAvg[1]*alphaR[5])); 
  Ghat[14] = 0.001190476190476191*((514.3928459844675*(fR[40]+fL[40])-296.98484809835*fR[27]+296.98484809835*fL[27])*amax+(94.86832980505142*alphaR[12]+148.492424049175*alphaR[1])*fAvg[18]+132.815661727072*alphaR[11]*fAvg[17]+(94.86832980505142*alphaR[8]+148.492424049175*alphaR[0])*fAvg[14]+148.492424049175*(alphaR[5]*fAvg[12]+fAvg[5]*alphaR[12])+132.815661727072*alphaR[4]*fAvg[10]+148.492424049175*(alphaR[3]*fAvg[8]+fAvg[3]*alphaR[8])+132.815661727072*alphaR[2]*fAvg[6]); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[28]+42.42640687119286*fL[28])*amax+(18.97366596101028*alphaR[11]+21.21320343559643*alphaR[2])*fAvg[19]+21.21320343559643*alphaR[4]*fAvg[16]+(18.97366596101028*alphaR[7]+21.21320343559643*alphaR[0])*fAvg[15]+16.97056274847715*(alphaR[5]*fAvg[13]+fAvg[5]*alphaR[13])+21.21320343559643*alphaR[1]*fAvg[9]+18.97366596101028*(alphaR[3]*fAvg[5]+fAvg[3]*alphaR[5])); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[43]+fL[43])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+(18.97366596101028*alphaR[12]+21.21320343559643*alphaR[1])*fAvg[19]+18.97366596101028*alphaR[13]*fAvg[17]+(18.97366596101028*alphaR[8]+21.21320343559643*alphaR[0])*fAvg[16]+21.21320343559643*alphaR[4]*fAvg[15]+18.97366596101028*alphaR[5]*fAvg[10]+21.21320343559643*alphaR[2]*fAvg[9]+18.97366596101028*alphaR[3]*fAvg[6]); 
  Ghat[17] = 2.380952380952381e-4*((2571.964229922338*(fR[44]+fL[44])-1484.92424049175*fR[37]+1484.92424049175*fL[37])*amax+593.9696961967002*(alphaR[5]*fAvg[19]+alphaR[4]*fAvg[18])+(664.0783086353599*alphaR[8]+474.3416490252571*alphaR[7]+742.462120245875*alphaR[0])*fAvg[17]+664.0783086353599*(alphaR[13]*fAvg[16]+alphaR[11]*fAvg[14])+(474.3416490252571*alphaR[11]+742.462120245875*alphaR[2])*fAvg[13]+(474.3416490252571*fAvg[11]+742.462120245875*fAvg[2])*alphaR[13]+593.9696961967002*fAvg[10]*alphaR[12]+742.462120245875*(alphaR[3]*fAvg[11]+fAvg[3]*alphaR[11])+664.0783086353599*alphaR[1]*fAvg[10]+742.462120245875*fAvg[6]*alphaR[7]+664.0783086353599*(alphaR[4]*fAvg[5]+fAvg[4]*alphaR[5])); 
  Ghat[18] = 2.380952380952381e-4*((2571.964229922338*(fR[46]+fL[46])-1484.92424049175*fR[39]+1484.92424049175*fL[39])*amax+(474.3416490252571*alphaR[8]+664.0783086353599*alphaR[7]+742.462120245875*alphaR[0])*fAvg[18]+593.9696961967002*alphaR[4]*fAvg[17]+(474.3416490252571*alphaR[12]+742.462120245875*alphaR[1])*fAvg[14]+664.0783086353599*(alphaR[12]*fAvg[13]+fAvg[12]*alphaR[13])+742.462120245875*(alphaR[3]*fAvg[12]+fAvg[3]*alphaR[12])+fAvg[10]*(593.9696961967002*alphaR[11]+664.0783086353599*alphaR[2])+742.462120245875*(alphaR[5]*fAvg[8]+fAvg[5]*alphaR[8])+664.0783086353599*alphaR[4]*fAvg[6]); 
  Ghat[19] = 0.001666666666666667*((367.4234614174769*(fR[47]+fL[47])-212.1320343559643*fR[42]+212.1320343559643*fL[42])*amax+(94.86832980505142*(alphaR[8]+alphaR[7])+106.0660171779821*alphaR[0])*fAvg[19]+84.85281374238573*alphaR[5]*fAvg[17]+(94.86832980505142*alphaR[12]+106.0660171779822*alphaR[1])*fAvg[16]+(94.86832980505142*alphaR[11]+106.0660171779822*alphaR[2])*fAvg[15]+fAvg[10]*(84.85281374238573*alphaR[13]+94.86832980505142*alphaR[3])+106.0660171779821*alphaR[4]*fAvg[9]+94.86832980505142*alphaR[5]*fAvg[6]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = -1.224744871391589*Ghat[0]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = 0.7071067811865475*Ghat[4]; 
  incr[7] = -1.224744871391589*Ghat[2]; 
  incr[8] = 0.7071067811865475*Ghat[5]; 
  incr[9] = -1.224744871391589*Ghat[3]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 0.7071067811865475*Ghat[7]; 
  incr[12] = 1.58113883008419*Ghat[0]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = 0.7071067811865475*Ghat[10]; 
  incr[18] = -1.224744871391589*Ghat[6]; 
  incr[19] = -1.224744871391589*Ghat[7]; 
  incr[20] = 1.58113883008419*Ghat[1]; 
  incr[21] = 0.7071067811865475*Ghat[11]; 
  incr[22] = 1.58113883008419*Ghat[2]; 
  incr[23] = 0.7071067811865475*Ghat[12]; 
  incr[24] = -1.224744871391589*Ghat[8]; 
  incr[25] = 0.7071067811865475*Ghat[13]; 
  incr[26] = 1.58113883008419*Ghat[3]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = 0.7071067811865475*Ghat[15]; 
  incr[29] = -1.224744871391589*Ghat[9]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = -1.224744871391589*Ghat[11]; 
  incr[33] = 1.58113883008419*Ghat[4]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = -1.224744871391589*Ghat[13]; 
  incr[36] = 1.58113883008419*Ghat[5]; 
  incr[37] = 0.7071067811865475*Ghat[17]; 
  incr[38] = 1.58113883008419*Ghat[6]; 
  incr[39] = 0.7071067811865475*Ghat[18]; 
  incr[40] = -1.224744871391589*Ghat[14]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = 0.7071067811865475*Ghat[19]; 
  incr[43] = -1.224744871391589*Ghat[16]; 
  incr[44] = -1.224744871391589*Ghat[17]; 
  incr[45] = 1.58113883008419*Ghat[10]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

  outR[0] += incr[0]*rdy2R; 
  outR[1] += incr[1]*rdy2R; 
  outR[2] += incr[2]*rdy2R; 
  outR[3] += incr[3]*rdy2R; 
  outR[4] += incr[4]*rdy2R; 
  outR[5] += incr[5]*rdy2R; 
  outR[6] += incr[6]*rdy2R; 
  outR[7] += incr[7]*rdy2R; 
  outR[8] += incr[8]*rdy2R; 
  outR[9] += incr[9]*rdy2R; 
  outR[10] += incr[10]*rdy2R; 
  outR[11] += incr[11]*rdy2R; 
  outR[12] += incr[12]*rdy2R; 
  outR[13] += incr[13]*rdy2R; 
  outR[14] += incr[14]*rdy2R; 
  outR[15] += incr[15]*rdy2R; 
  outR[16] += incr[16]*rdy2R; 
  outR[17] += incr[17]*rdy2R; 
  outR[18] += incr[18]*rdy2R; 
  outR[19] += incr[19]*rdy2R; 
  outR[20] += incr[20]*rdy2R; 
  outR[21] += incr[21]*rdy2R; 
  outR[22] += incr[22]*rdy2R; 
  outR[23] += incr[23]*rdy2R; 
  outR[24] += incr[24]*rdy2R; 
  outR[25] += incr[25]*rdy2R; 
  outR[26] += incr[26]*rdy2R; 
  outR[27] += incr[27]*rdy2R; 
  outR[28] += incr[28]*rdy2R; 
  outR[29] += incr[29]*rdy2R; 
  outR[30] += incr[30]*rdy2R; 
  outR[31] += incr[31]*rdy2R; 
  outR[32] += incr[32]*rdy2R; 
  outR[33] += incr[33]*rdy2R; 
  outR[34] += incr[34]*rdy2R; 
  outR[35] += incr[35]*rdy2R; 
  outR[36] += incr[36]*rdy2R; 
  outR[37] += incr[37]*rdy2R; 
  outR[38] += incr[38]*rdy2R; 
  outR[39] += incr[39]*rdy2R; 
  outR[40] += incr[40]*rdy2R; 
  outR[41] += incr[41]*rdy2R; 
  outR[42] += incr[42]*rdy2R; 
  outR[43] += incr[43]*rdy2R; 
  outR[44] += incr[44]*rdy2R; 
  outR[45] += incr[45]*rdy2R; 
  outR[46] += incr[46]*rdy2R; 
  outR[47] += incr[47]*rdy2R; 

  outL[0] += -1.0*incr[0]*rdy2L; 
  outL[1] += -1.0*incr[1]*rdy2L; 
  outL[2] += incr[2]*rdy2L; 
  outL[3] += -1.0*incr[3]*rdy2L; 
  outL[4] += -1.0*incr[4]*rdy2L; 
  outL[5] += incr[5]*rdy2L; 
  outL[6] += -1.0*incr[6]*rdy2L; 
  outL[7] += incr[7]*rdy2L; 
  outL[8] += -1.0*incr[8]*rdy2L; 
  outL[9] += incr[9]*rdy2L; 
  outL[10] += -1.0*incr[10]*rdy2L; 
  outL[11] += -1.0*incr[11]*rdy2L; 
  outL[12] += -1.0*incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += -1.0*incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 
  outL[16] += incr[16]*rdy2L; 
  outL[17] += -1.0*incr[17]*rdy2L; 
  outL[18] += incr[18]*rdy2L; 
  outL[19] += incr[19]*rdy2L; 
  outL[20] += -1.0*incr[20]*rdy2L; 
  outL[21] += -1.0*incr[21]*rdy2L; 
  outL[22] += -1.0*incr[22]*rdy2L; 
  outL[23] += -1.0*incr[23]*rdy2L; 
  outL[24] += incr[24]*rdy2L; 
  outL[25] += -1.0*incr[25]*rdy2L; 
  outL[26] += -1.0*incr[26]*rdy2L; 
  outL[27] += -1.0*incr[27]*rdy2L; 
  outL[28] += -1.0*incr[28]*rdy2L; 
  outL[29] += incr[29]*rdy2L; 
  outL[30] += -1.0*incr[30]*rdy2L; 
  outL[31] += incr[31]*rdy2L; 
  outL[32] += incr[32]*rdy2L; 
  outL[33] += -1.0*incr[33]*rdy2L; 
  outL[34] += incr[34]*rdy2L; 
  outL[35] += incr[35]*rdy2L; 
  outL[36] += -1.0*incr[36]*rdy2L; 
  outL[37] += -1.0*incr[37]*rdy2L; 
  outL[38] += -1.0*incr[38]*rdy2L; 
  outL[39] += -1.0*incr[39]*rdy2L; 
  outL[40] += incr[40]*rdy2L; 
  outL[41] += incr[41]*rdy2L; 
  outL[42] += -1.0*incr[42]*rdy2L; 
  outL[43] += incr[43]*rdy2L; 
  outL[44] += incr[44]*rdy2L; 
  outL[45] += -1.0*incr[45]*rdy2L; 
  outL[46] += incr[46]*rdy2L; 
  outL[47] += incr[47]*rdy2L; 
return std::abs(alphaSurfAvgR); 
} 
double GyrokineticGenGeoSurf2x2vSer_vpar_P2_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[8] = (1.154700538379252*bmag[1])/rdmu2R; 
  hamilR[11] = 2.0*(bmag[4]*wmuR+phi[4]*q_); 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 
  hamilR[25] = (1.154700538379251*bmag[4])/rdmu2R; 

  double BstarXdBmagR[48]; 

  double BstarYdBmagR[48]; 
  BstarYdBmagR[0] = -(1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[1] = -(1.732050807568877*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[3] = -(1.0*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[6] = -(1.0*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[11] = -(1.732050807568877*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[21] = -(1.0*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[20]; 
  alphaR[0] = (0.03535533905932736*(hamilR[19]*(15.0*BstarYdBmagR[21]-8.660254037844387*BstarYdBmagR[11])+hamilR[5]*(15.0*BstarYdBmagR[6]-8.660254037844386*BstarYdBmagR[1])+hamilR[2]*(15.0*BstarYdBmagR[3]-8.660254037844386*BstarYdBmagR[0]))*rdy2R)/m_; 
  alphaR[1] = (0.03535533905932736*(13.41640786499874*hamilR[5]*BstarYdBmagR[21]+(13.41640786499874*BstarYdBmagR[6]-7.745966692414834*BstarYdBmagR[1])*hamilR[19]-7.745966692414834*hamilR[5]*BstarYdBmagR[11]+15.0*hamilR[2]*BstarYdBmagR[6]+(15.0*BstarYdBmagR[3]-8.660254037844386*BstarYdBmagR[0])*hamilR[5]-8.660254037844386*BstarYdBmagR[1]*hamilR[2])*rdy2R)/m_; 
  alphaR[2] = (0.1767766952966368*((6.708203932499369*BstarYdBmagR[6]-3.872983346207417*BstarYdBmagR[1])*hamilR[20]+(6.708203932499369*BstarYdBmagR[3]-3.872983346207417*BstarYdBmagR[0])*hamilR[12])*rdy2R)/m_; 
  alphaR[4] = (0.03535533905932736*(hamilR[20]*(30.0*BstarYdBmagR[21]-17.32050807568877*BstarYdBmagR[11]+33.54101966249684*BstarYdBmagR[3]-19.36491673103708*BstarYdBmagR[0])+(33.54101966249685*BstarYdBmagR[6]-19.36491673103709*BstarYdBmagR[1])*hamilR[12])*rdy2R)/m_; 
  alphaR[7] = (0.005050762722761052*((67.0820393249937*hamilR[19]+105.0*hamilR[2])*BstarYdBmagR[21]+((-38.72983346207417*BstarYdBmagR[11])+105.0*BstarYdBmagR[3]-60.62177826491071*BstarYdBmagR[0])*hamilR[19]-60.6217782649107*hamilR[2]*BstarYdBmagR[11]+hamilR[5]*(93.91485505499116*BstarYdBmagR[6]-54.22176684690384*BstarYdBmagR[1]))*rdy2R)/m_; 
  alphaR[11] = (0.1767766952966368*(6.708203932499369*hamilR[12]*BstarYdBmagR[21]+(6.0*BstarYdBmagR[6]-3.464101615137754*BstarYdBmagR[1])*hamilR[20]-3.872983346207417*BstarYdBmagR[11]*hamilR[12])*rdy2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.00625*(15.0*hamilR[19]*BstarYdBmagR[21]-8.660254037844387*BstarYdBmagR[11]*hamilR[19]+15.0*hamilR[5]*BstarYdBmagR[6]-8.660254037844386*BstarYdBmagR[1]*hamilR[5]+15.0*hamilR[2]*BstarYdBmagR[3]-8.660254037844386*BstarYdBmagR[0]*hamilR[2])*rdy2R)/m_; 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[13]+fL[13])+1.732050807568877*(fL[3]-1.0*fR[3])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[23]+fL[23])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[24]+fL[24])+3.0*(fL[7]-1.0*fR[7]))+3.0*(fR[2]+fL[2])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[27]+fL[27])+3.0*(fL[10]-1.0*fR[10]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[34]+fL[34])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[5]+fL[5]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[39]+fL[39])+1.732050807568877*(fL[17]-1.0*fR[17])+fR[8]+fL[8]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[40]+fL[40])+1.732050807568877*(fL[18]-1.0*fR[18])+fR[9]+fL[9]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[21]-1.0*(8.660254037844387*fL[21]+5.0*(fR[11]+fL[11]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[22]-1.0*(8.660254037844387*fL[22]+5.0*(fR[12]+fL[12]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[30]-1.0*(8.660254037844387*fL[30]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[46]+fL[46])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[16]+fL[16])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[32]-1.0*(8.660254037844387*fL[32]+5.0*(fR[19]+fL[19]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[33]-1.0*(8.660254037844387*fL[33]+5.0*(fR[20]+fL[20]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[37]-1.0*(8.660254037844387*fL[37]+5.0*(fR[25]+fL[25]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[38]-1.0*(8.660254037844387*fL[38]+5.0*(fR[26]+fL[26]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[42]-1.0*(8.660254037844387*fL[42]+5.0*(fR[28]+fL[28]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[43]-1.0*(8.660254037844387*fL[43]+5.0*(fR[29]+fL[29]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[44]-1.0*(8.660254037844387*fL[44]+5.0*(fR[35]+fL[35]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[45]-1.0*(8.660254037844387*fL[45]+5.0*(fR[36]+fL[36]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[41]+fL[41]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[13]-6.324555320336761*fL[13]-4.898979485566357*(fR[3]+fL[3])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[11]*fAvg[11]+alphaR[7]*fAvg[7]+alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2]+alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[23]-94.8683298050514*fL[23]-73.48469228349535*(fR[6]+fL[6])+42.42640687119286*fR[1]-42.42640687119286*fL[1])*amax-18.97366596101028*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11]+alphaR[1]*fAvg[7]+fAvg[1]*alphaR[7])-21.21320343559643*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4]+alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[24]-94.8683298050514*fL[24]-73.48469228349535*(fR[7]+fL[7])+42.42640687119286*fR[2]-42.42640687119286*fL[2])*amax-18.97366596101028*alphaR[4]*fAvg[12]-21.21320343559643*(alphaR[7]*fAvg[11]+fAvg[7]*alphaR[11])-18.97366596101028*alphaR[2]*fAvg[8]-21.21320343559643*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4]+alphaR[0]*fAvg[2]+fAvg[0]*alphaR[2])); 
  Ghat[3] = -0.008333333333333333*((94.8683298050514*fR[27]-94.8683298050514*fL[27]-73.48469228349535*(fR[10]+fL[10])+42.42640687119286*fR[4]-42.42640687119286*fL[4])*amax-21.21320343559643*(alphaR[11]*fAvg[17]+alphaR[7]*fAvg[13])-21.21320343559643*(alphaR[4]*fAvg[10]+alphaR[2]*fAvg[6]+alphaR[1]*fAvg[5]+alphaR[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[34]-94.86832980505142*fL[34]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[5]-42.42640687119286*fL[5])*amax+((-16.97056274847715*alphaR[11])-18.97366596101028*alphaR[2])*fAvg[12]-18.97366596101028*(alphaR[1]*fAvg[11]+fAvg[1]*alphaR[11]+alphaR[4]*(fAvg[8]+fAvg[7])+fAvg[4]*alphaR[7])-21.21320343559643*(alphaR[0]*fAvg[4]+fAvg[0]*alphaR[4]+alphaR[1]*fAvg[2]+fAvg[1]*alphaR[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[39]-94.86832980505142*fL[39]-73.48469228349535*(fR[17]+fL[17])+42.42640687119286*fR[8]-42.42640687119286*fL[8])*amax-18.97366596101028*(alphaR[4]*fAvg[17]+alphaR[1]*fAvg[13])+fAvg[10]*((-18.97366596101028*alphaR[11])-21.21320343559643*alphaR[2])-18.97366596101028*fAvg[5]*alphaR[7]-21.21320343559643*(alphaR[4]*fAvg[6]+alphaR[0]*fAvg[5]+alphaR[1]*fAvg[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[40]-94.86832980505142*fL[40]-73.48469228349535*(fR[18]+fL[18])+42.42640687119286*fR[9]-42.42640687119286*fL[9])*amax-18.97366596101028*alphaR[4]*fAvg[18]-21.21320343559643*alphaR[7]*fAvg[17]-18.97366596101028*alphaR[2]*fAvg[14]-21.21320343559643*(alphaR[11]*fAvg[13]+alphaR[1]*fAvg[10]+alphaR[0]*fAvg[6]+alphaR[4]*fAvg[5]+alphaR[2]*fAvg[3])); 
  Ghat[7] = 0.001190476190476191*((514.3928459844675*(fR[21]+fL[21])-296.98484809835*fR[11]+296.98484809835*fL[11])*amax+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[11]+148.492424049175*fAvg[2]*alphaR[11]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[7]+148.492424049175*fAvg[0]*alphaR[7]+132.815661727072*(alphaR[4]*fAvg[4]+alphaR[1]*fAvg[1])); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[22]+fL[22])-42.42640687119286*fR[12]+42.42640687119286*fL[12])*amax+21.21320343559643*alphaR[1]*fAvg[12]+18.97366596101028*alphaR[11]*fAvg[11]+21.21320343559643*alphaR[0]*fAvg[8]+18.97366596101028*(alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[30]+fL[30])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[4]*fAvg[19]+21.21320343559643*(alphaR[2]*fAvg[16]+alphaR[1]*fAvg[15])+21.21320343559643*alphaR[0]*fAvg[9]); 
  Ghat[10] = -0.001666666666666667*((474.341649025257*fR[46]-474.341649025257*fL[46]-367.4234614174767*(fR[31]+fL[31])+212.1320343559643*fR[16]-212.1320343559643*fL[16])*amax+((-84.85281374238573*alphaR[11])-94.86832980505142*alphaR[2])*fAvg[18]-94.86832980505142*alphaR[1]*fAvg[17]-94.8683298050514*(alphaR[4]*(fAvg[14]+fAvg[13])+fAvg[5]*alphaR[11])+((-94.86832980505142*alphaR[7])-106.0660171779821*alphaR[0])*fAvg[10]-106.0660171779821*(alphaR[1]*fAvg[6]+alphaR[2]*fAvg[5]+fAvg[3]*alphaR[4])); 
  Ghat[11] = 0.001190476190476191*((514.3928459844675*(fR[32]+fL[32])-296.98484809835*fR[19]+296.98484809835*fL[19])*amax+118.79393923934*alphaR[4]*fAvg[12]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[11]+(132.815661727072*fAvg[8]+94.86832980505142*fAvg[7]+148.492424049175*fAvg[0])*alphaR[11]+148.492424049175*(alphaR[2]*fAvg[7]+fAvg[2]*alphaR[7])+132.815661727072*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4])); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[33]+fL[33])-42.42640687119286*fR[20]+42.42640687119286*fL[20])*amax+(18.97366596101028*alphaR[7]+21.21320343559643*alphaR[0])*fAvg[12]+16.97056274847715*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11])+21.21320343559643*alphaR[1]*fAvg[8]+18.97366596101028*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4])); 
  Ghat[13] = 0.001190476190476191*((514.3928459844675*(fR[37]+fL[37])-296.98484809835*fR[25]+296.98484809835*fL[25])*amax+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[17]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[13]+148.492424049175*fAvg[6]*alphaR[11]+132.815661727072*alphaR[4]*fAvg[10]+148.492424049175*fAvg[3]*alphaR[7]+132.815661727072*alphaR[1]*fAvg[5]); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[38]+fL[38])-42.42640687119286*fR[26]+42.42640687119286*fL[26])*amax+21.21320343559643*alphaR[1]*fAvg[18]+18.97366596101028*alphaR[11]*fAvg[17]+21.21320343559643*alphaR[0]*fAvg[14]+18.97366596101028*(alphaR[4]*fAvg[10]+alphaR[2]*fAvg[6])); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[42]+fL[42])-42.42640687119286*fR[28]+42.42640687119286*fL[28])*amax+(18.97366596101028*alphaR[11]+21.21320343559643*alphaR[2])*fAvg[19]+21.21320343559643*alphaR[4]*fAvg[16]+(18.97366596101028*alphaR[7]+21.21320343559643*alphaR[0])*fAvg[15]+21.21320343559643*alphaR[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[43]+fL[43])-42.42640687119286*fR[29]+42.42640687119286*fL[29])*amax+21.21320343559643*alphaR[1]*fAvg[19]+21.21320343559643*(alphaR[0]*fAvg[16]+alphaR[4]*fAvg[15])+21.21320343559643*alphaR[2]*fAvg[9]); 
  Ghat[17] = 0.001190476190476191*((514.3928459844675*(fR[44]+fL[44])-296.98484809835*fR[35]+296.98484809835*fL[35])*amax+118.79393923934*alphaR[4]*fAvg[18]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[17]+132.815661727072*alphaR[11]*fAvg[14]+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[13]+148.492424049175*fAvg[3]*alphaR[11]+132.815661727072*alphaR[1]*fAvg[10]+148.492424049175*fAvg[6]*alphaR[7]+132.815661727072*alphaR[4]*fAvg[5]); 
  Ghat[18] = 0.001666666666666667*((367.4234614174769*(fR[45]+fL[45])-212.1320343559643*fR[36]+212.1320343559643*fL[36])*amax+(94.86832980505142*alphaR[7]+106.0660171779821*alphaR[0])*fAvg[18]+84.85281374238573*alphaR[4]*fAvg[17]+106.0660171779822*alphaR[1]*fAvg[14]+84.85281374238573*fAvg[10]*alphaR[11]+94.86832980505142*(alphaR[2]*fAvg[10]+alphaR[4]*fAvg[6])); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[41]+42.42640687119286*fL[41])*amax+(18.97366596101028*alphaR[7]+21.21320343559643*alphaR[0])*fAvg[19]+21.21320343559643*alphaR[1]*fAvg[16]+(18.97366596101028*alphaR[11]+21.21320343559643*alphaR[2])*fAvg[15]+21.21320343559643*alphaR[4]*fAvg[9]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = 0.7071067811865475*Ghat[2]; 
  incr[3] = -1.224744871391589*Ghat[0]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = 0.7071067811865475*Ghat[4]; 
  incr[6] = -1.224744871391589*Ghat[1]; 
  incr[7] = -1.224744871391589*Ghat[2]; 
  incr[8] = 0.7071067811865475*Ghat[5]; 
  incr[9] = 0.7071067811865475*Ghat[6]; 
  incr[10] = -1.224744871391589*Ghat[3]; 
  incr[11] = 0.7071067811865475*Ghat[7]; 
  incr[12] = 0.7071067811865475*Ghat[8]; 
  incr[13] = 1.58113883008419*Ghat[0]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = 0.7071067811865475*Ghat[10]; 
  incr[17] = -1.224744871391589*Ghat[5]; 
  incr[18] = -1.224744871391589*Ghat[6]; 
  incr[19] = 0.7071067811865475*Ghat[11]; 
  incr[20] = 0.7071067811865475*Ghat[12]; 
  incr[21] = -1.224744871391589*Ghat[7]; 
  incr[22] = -1.224744871391589*Ghat[8]; 
  incr[23] = 1.58113883008419*Ghat[1]; 
  incr[24] = 1.58113883008419*Ghat[2]; 
  incr[25] = 0.7071067811865475*Ghat[13]; 
  incr[26] = 0.7071067811865475*Ghat[14]; 
  incr[27] = 1.58113883008419*Ghat[3]; 
  incr[28] = 0.7071067811865475*Ghat[15]; 
  incr[29] = 0.7071067811865475*Ghat[16]; 
  incr[30] = -1.224744871391589*Ghat[9]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = -1.224744871391589*Ghat[11]; 
  incr[33] = -1.224744871391589*Ghat[12]; 
  incr[34] = 1.58113883008419*Ghat[4]; 
  incr[35] = 0.7071067811865475*Ghat[17]; 
  incr[36] = 0.7071067811865475*Ghat[18]; 
  incr[37] = -1.224744871391589*Ghat[13]; 
  incr[38] = -1.224744871391589*Ghat[14]; 
  incr[39] = 1.58113883008419*Ghat[5]; 
  incr[40] = 1.58113883008419*Ghat[6]; 
  incr[41] = 0.7071067811865475*Ghat[19]; 
  incr[42] = -1.224744871391589*Ghat[15]; 
  incr[43] = -1.224744871391589*Ghat[16]; 
  incr[44] = -1.224744871391589*Ghat[17]; 
  incr[45] = -1.224744871391589*Ghat[18]; 
  incr[46] = 1.58113883008419*Ghat[10]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 
  outR[8] += incr[8]*rdvpar2R; 
  outR[9] += incr[9]*rdvpar2R; 
  outR[10] += incr[10]*rdvpar2R; 
  outR[11] += incr[11]*rdvpar2R; 
  outR[12] += incr[12]*rdvpar2R; 
  outR[13] += incr[13]*rdvpar2R; 
  outR[14] += incr[14]*rdvpar2R; 
  outR[15] += incr[15]*rdvpar2R; 
  outR[16] += incr[16]*rdvpar2R; 
  outR[17] += incr[17]*rdvpar2R; 
  outR[18] += incr[18]*rdvpar2R; 
  outR[19] += incr[19]*rdvpar2R; 
  outR[20] += incr[20]*rdvpar2R; 
  outR[21] += incr[21]*rdvpar2R; 
  outR[22] += incr[22]*rdvpar2R; 
  outR[23] += incr[23]*rdvpar2R; 
  outR[24] += incr[24]*rdvpar2R; 
  outR[25] += incr[25]*rdvpar2R; 
  outR[26] += incr[26]*rdvpar2R; 
  outR[27] += incr[27]*rdvpar2R; 
  outR[28] += incr[28]*rdvpar2R; 
  outR[29] += incr[29]*rdvpar2R; 
  outR[30] += incr[30]*rdvpar2R; 
  outR[31] += incr[31]*rdvpar2R; 
  outR[32] += incr[32]*rdvpar2R; 
  outR[33] += incr[33]*rdvpar2R; 
  outR[34] += incr[34]*rdvpar2R; 
  outR[35] += incr[35]*rdvpar2R; 
  outR[36] += incr[36]*rdvpar2R; 
  outR[37] += incr[37]*rdvpar2R; 
  outR[38] += incr[38]*rdvpar2R; 
  outR[39] += incr[39]*rdvpar2R; 
  outR[40] += incr[40]*rdvpar2R; 
  outR[41] += incr[41]*rdvpar2R; 
  outR[42] += incr[42]*rdvpar2R; 
  outR[43] += incr[43]*rdvpar2R; 
  outR[44] += incr[44]*rdvpar2R; 
  outR[45] += incr[45]*rdvpar2R; 
  outR[46] += incr[46]*rdvpar2R; 
  outR[47] += incr[47]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += -1.0*incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 
  outL[4] += -1.0*incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 
  outL[8] += -1.0*incr[8]*rdvpar2L; 
  outL[9] += -1.0*incr[9]*rdvpar2L; 
  outL[10] += incr[10]*rdvpar2L; 
  outL[11] += -1.0*incr[11]*rdvpar2L; 
  outL[12] += -1.0*incr[12]*rdvpar2L; 
  outL[13] += -1.0*incr[13]*rdvpar2L; 
  outL[14] += -1.0*incr[14]*rdvpar2L; 
  outL[15] += incr[15]*rdvpar2L; 
  outL[16] += -1.0*incr[16]*rdvpar2L; 
  outL[17] += incr[17]*rdvpar2L; 
  outL[18] += incr[18]*rdvpar2L; 
  outL[19] += -1.0*incr[19]*rdvpar2L; 
  outL[20] += -1.0*incr[20]*rdvpar2L; 
  outL[21] += incr[21]*rdvpar2L; 
  outL[22] += incr[22]*rdvpar2L; 
  outL[23] += -1.0*incr[23]*rdvpar2L; 
  outL[24] += -1.0*incr[24]*rdvpar2L; 
  outL[25] += -1.0*incr[25]*rdvpar2L; 
  outL[26] += -1.0*incr[26]*rdvpar2L; 
  outL[27] += -1.0*incr[27]*rdvpar2L; 
  outL[28] += -1.0*incr[28]*rdvpar2L; 
  outL[29] += -1.0*incr[29]*rdvpar2L; 
  outL[30] += incr[30]*rdvpar2L; 
  outL[31] += incr[31]*rdvpar2L; 
  outL[32] += incr[32]*rdvpar2L; 
  outL[33] += incr[33]*rdvpar2L; 
  outL[34] += -1.0*incr[34]*rdvpar2L; 
  outL[35] += -1.0*incr[35]*rdvpar2L; 
  outL[36] += -1.0*incr[36]*rdvpar2L; 
  outL[37] += incr[37]*rdvpar2L; 
  outL[38] += incr[38]*rdvpar2L; 
  outL[39] += -1.0*incr[39]*rdvpar2L; 
  outL[40] += -1.0*incr[40]*rdvpar2L; 
  outL[41] += -1.0*incr[41]*rdvpar2L; 
  outL[42] += incr[42]*rdvpar2L; 
  outL[43] += incr[43]*rdvpar2L; 
  outL[44] += incr[44]*rdvpar2L; 
  outL[45] += incr[45]*rdvpar2L; 
  outL[46] += -1.0*incr[46]*rdvpar2L; 
  outL[47] += incr[47]*rdvpar2L; 
return std::abs(alphaSurfAvgR); 
} 