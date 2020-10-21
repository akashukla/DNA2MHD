#include <VlasovModDecl.h> 
__host__ __device__ double VlasovSurfElcMag2x2vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[2]; 
  double dv10r = 2/dxvr[2]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double *B2 = &EM[50]; 

  double Ghat[20]; 
  double favg[20]; 
  double alpha[20]; 

  favg[0] = (-1.870828693386971*fr[33])+1.870828693386971*fl[33]+1.58113883008419*fr[13]+1.58113883008419*fl[13]-1.224744871391589*fr[3]+1.224744871391589*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = 1.58113883008419*fr[23]+1.58113883008419*fl[23]-1.224744871391589*fr[6]+1.224744871391589*fl[6]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = 1.58113883008419*fr[24]+1.58113883008419*fl[24]-1.224744871391589*fr[7]+1.224744871391589*fl[7]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = 1.58113883008419*fr[27]+1.58113883008419*fl[27]-1.224744871391589*fr[10]+1.224744871391589*fl[10]+0.7071067811865475*fr[4]+0.7071067811865475*fl[4]; 
  favg[4] = (-1.224744871391589*fr[15])+1.224744871391589*fl[15]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = (-1.224744871391589*fr[17])+1.224744871391589*fl[17]+0.7071067811865475*fr[8]+0.7071067811865475*fl[8]; 
  favg[6] = (-1.224744871391589*fr[18])+1.224744871391589*fl[18]+0.7071067811865475*fr[9]+0.7071067811865475*fl[9]; 
  favg[7] = (-1.224744871391589*fr[21])+1.224744871391589*fl[21]+0.7071067811865475*fr[11]+0.7071067811865475*fl[11]; 
  favg[8] = (-1.224744871391589*fr[22])+1.224744871391589*fl[22]+0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 
  favg[9] = (-1.224744871391589*fr[30])+1.224744871391589*fl[30]+0.7071067811865475*fr[14]+0.7071067811865475*fl[14]; 
  favg[10] = 0.7071067811865475*fr[16]+0.7071067811865475*fl[16]; 
  favg[11] = 0.7071067811865475*fr[19]+0.7071067811865475*fl[19]; 
  favg[12] = 0.7071067811865475*fr[20]+0.7071067811865475*fl[20]; 
  favg[13] = 0.7071067811865475*fr[25]+0.7071067811865475*fl[25]; 
  favg[14] = 0.7071067811865475*fr[26]+0.7071067811865475*fl[26]; 
  favg[15] = 0.7071067811865475*fr[28]+0.7071067811865475*fl[28]; 
  favg[16] = 0.7071067811865475*fr[29]+0.7071067811865475*fl[29]; 
  favg[17] = 0.7071067811865475*fr[31]+0.7071067811865475*fl[31]; 
  favg[18] = 0.7071067811865475*fr[32]+0.7071067811865475*fl[32]; 
  favg[19] = 0.7071067811865475*fr[34]+0.7071067811865475*fl[34]; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 1.414213562373095*(B2[2]*wv2+E0[2]); 
  alpha[3] = 0.408248290463863*B2[0]*dv2; 
  alpha[4] = 1.414213562373095*(B2[3]*wv2+E0[3]); 
  alpha[5] = 0.408248290463863*B2[1]*dv2; 
  alpha[6] = 0.408248290463863*B2[2]*dv2; 
  alpha[7] = 1.414213562373095*(B2[4]*wv2+E0[4]); 
  alpha[8] = 1.414213562373095*(B2[5]*wv2+E0[5]); 
  alpha[10] = 0.408248290463863*B2[3]*dv2; 
  alpha[11] = 1.414213562373095*(B2[6]*wv2+E0[6]); 
  alpha[12] = 1.414213562373095*(B2[7]*wv2+E0[7]); 
  alpha[13] = 0.408248290463863*B2[4]*dv2; 
  alpha[14] = 0.408248290463863*B2[5]*dv2; 
  alpha[17] = 1.414213562373095*(B2[8]*wv2+E0[8]); 
  alpha[18] = 1.414213562373095*(B2[9]*wv2+E0[9]); 

  const double amid = (-0.3952847075210473*alpha[8])-0.3952847075210473*alpha[7]+0.3535533905932737*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(2.645751311064591*(fr[33]+fl[33])+2.23606797749979*(fl[13]-1.0*fr[13])+1.732050807568877*(fr[3]+fl[3])-1.0*fr[0]+fl[0])*amax+0.1767766952966368*(alpha[18]*favg[18]+alpha[17]*favg[17]+alpha[14]*favg[14]+alpha[13]*favg[13]+alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[10]*favg[10]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.001683587574253684*(92.22255689363637*(alpha[7]*favg[17]+favg[7]*alpha[17])+93.91485505499116*(alpha[5]*favg[13]+favg[5]*alpha[13])+105.0*(alpha[8]*favg[12]+favg[8]*alpha[12])+93.91485505499116*(alpha[4]*favg[11]+favg[4]*alpha[11])+105.0*(alpha[6]*favg[10]+favg[6]*alpha[10])+93.91485505499116*(alpha[1]*favg[7]+favg[1]*alpha[7])+105.0*(alpha[3]*favg[5]+favg[3]*alpha[5]+alpha[2]*favg[4]+favg[2]*alpha[4]+alpha[0]*favg[1]+favg[0]*alpha[1]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[23]-1.0*(3.872983346207417*fl[23]+3.0*(fr[6]+fl[6])))+3.0*(fr[1]-1.0*fl[1]))*amax; 
  Ghat[2] = 0.001683587574253684*(92.22255689363637*(alpha[8]*favg[18]+favg[8]*alpha[18])+93.91485505499116*(alpha[6]*favg[14]+favg[6]*alpha[14]+alpha[4]*favg[12]+favg[4]*alpha[12])+105.0*(alpha[7]*favg[11]+favg[7]*alpha[11])+105.0*(alpha[5]*favg[10]+favg[5]*alpha[10])+93.91485505499116*(alpha[2]*favg[8]+favg[2]*alpha[8])+105.0*(alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[1]*favg[4]+favg[1]*alpha[4]+alpha[0]*favg[2]+favg[0]*alpha[2]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[24]-1.0*(3.872983346207417*fl[24]+3.0*(fr[7]+fl[7])))+3.0*(fr[2]-1.0*fl[2]))*amax; 
  Ghat[3] = 0.01178511301977579*(13.41640786499874*(alpha[6]*favg[16]+alpha[5]*favg[15])+15.0*(alpha[8]*favg[14]+favg[8]*alpha[14]+alpha[7]*favg[13]+favg[7]*alpha[13])+15.0*(alpha[4]*favg[10]+favg[4]*alpha[10])+13.41640786499874*alpha[3]*favg[9]+15.0*(alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[3]+favg[0]*alpha[3]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[27]-1.0*(3.872983346207417*fl[27]+3.0*(fr[10]+fl[10])))+3.0*(fr[4]-1.0*fl[4]))*amax; 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[15]+fl[15])-1.0*fr[5]+fl[5])*amax+0.001683587574253684*(92.2225568936364*(alpha[12]*favg[18]+favg[12]*alpha[18]+alpha[11]*favg[17]+favg[11]*alpha[17])+93.91485505499116*(alpha[10]*favg[14]+favg[10]*alpha[14]+alpha[10]*favg[13]+favg[10]*alpha[13])+(84.0*alpha[11]+93.91485505499116*alpha[2])*favg[12]+(84.0*favg[11]+93.91485505499116*favg[2])*alpha[12]+93.91485505499116*(alpha[1]*favg[11]+favg[1]*alpha[11])+105.0*(alpha[3]*favg[10]+favg[3]*alpha[10])+93.91485505499116*(alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[4]*favg[7]+favg[4]*alpha[7])+105.0*(alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[0]*favg[4]+favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2])); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[17]+fl[17])-1.0*fr[8]+fl[8])*amax+0.001683587574253684*(92.2225568936364*(alpha[13]*favg[17]+favg[13]*alpha[17])+93.91485505499116*alpha[10]*favg[16]+(84.0*alpha[13]+93.91485505499116*alpha[3])*favg[15]+105.0*(alpha[12]*favg[14]+favg[12]*alpha[14])+93.91485505499116*(alpha[1]*favg[13]+favg[1]*alpha[13]+alpha[10]*favg[11]+favg[10]*alpha[11])+105.0*(alpha[2]*favg[10]+favg[2]*alpha[10])+93.91485505499116*(alpha[5]*(favg[9]+favg[7])+favg[5]*alpha[7])+105.0*(alpha[4]*favg[6]+favg[4]*alpha[6]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[3]+favg[1]*alpha[3])); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[18]+fl[18])-1.0*fr[9]+fl[9])*amax+0.001683587574253684*(92.2225568936364*(alpha[14]*favg[18]+favg[14]*alpha[18])+(84.0*alpha[14]+93.91485505499116*alpha[3])*favg[16]+93.91485505499116*(alpha[10]*favg[15]+alpha[2]*favg[14]+favg[2]*alpha[14])+105.0*(alpha[11]*favg[13]+favg[11]*alpha[13])+93.91485505499116*(alpha[10]*favg[12]+favg[10]*alpha[12])+105.0*(alpha[1]*favg[10]+favg[1]*alpha[10])+93.91485505499116*(alpha[6]*(favg[9]+favg[8])+favg[6]*alpha[8])+105.0*(alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[4]*favg[5]+favg[4]*alpha[5]+alpha[2]*favg[3]+favg[2]*alpha[3])); 
  Ghat[7] = 0.07071067811865475*(8.660254037844387*(fr[21]+fl[21])+5.0*(fl[11]-1.0*fr[11]))*amax+0.001683587574253684*((62.60990336999411*alpha[17]+92.22255689363637*alpha[1])*favg[17]+92.22255689363637*favg[1]*alpha[17]+(67.0820393249937*alpha[13]+105.0*alpha[3])*favg[13]+105.0*favg[3]*alpha[13]+93.91485505499116*alpha[12]*favg[12]+(67.0820393249937*alpha[11]+105.0*alpha[2])*favg[11]+105.0*favg[2]*alpha[11]+93.91485505499116*alpha[10]*favg[10]+(67.0820393249937*alpha[7]+105.0*alpha[0])*favg[7]+105.0*favg[0]*alpha[7]+93.91485505499116*(alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[1]*favg[1])); 
  Ghat[8] = 0.07071067811865475*(8.660254037844387*(fr[22]+fl[22])+5.0*(fl[12]-1.0*fr[12]))*amax+0.001683587574253684*((62.60990336999411*alpha[18]+92.22255689363637*alpha[2])*favg[18]+92.22255689363637*favg[2]*alpha[18]+(67.0820393249937*alpha[14]+105.0*alpha[3])*favg[14]+105.0*favg[3]*alpha[14]+(67.0820393249937*alpha[12]+105.0*alpha[1])*favg[12]+105.0*favg[1]*alpha[12]+93.91485505499116*(alpha[11]*favg[11]+alpha[10]*favg[10])+(67.0820393249937*alpha[8]+105.0*alpha[0])*favg[8]+105.0*favg[0]*alpha[8]+93.91485505499116*(alpha[6]*favg[6]+alpha[4]*favg[4]+alpha[2]*favg[2])); 
  Ghat[9] = 0.07071067811865475*(8.660254037844387*(fr[30]+fl[30])+5.0*(fl[14]-1.0*fr[14]))*amax+0.001683587574253684*(92.22255689363637*alpha[3]*favg[19]+105.0*(alpha[2]*favg[16]+alpha[1]*favg[15])+93.91485505499116*(alpha[14]*favg[14]+alpha[13]*favg[13]+alpha[10]*favg[10])+105.0*alpha[0]*favg[9]+93.91485505499116*(alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[3]*favg[3])); 
  Ghat[10] = 0.01178511301977579*(13.41640786499874*(alpha[5]*favg[16]+alpha[6]*favg[15]+alpha[4]*favg[14]+favg[4]*alpha[14]+alpha[4]*favg[13]+favg[4]*alpha[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[5]*favg[11]+favg[5]*alpha[11])+(13.41640786499874*(alpha[8]+alpha[7])+15.0*alpha[0])*favg[10]+(13.41640786499874*(favg[9]+favg[8]+favg[7])+15.0*favg[0])*alpha[10]+15.0*(alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[3]*favg[4]+favg[3]*alpha[4]))-0.3535533905932737*(fr[16]-1.0*fl[16])*amax; 
  Ghat[11] = 0.001683587574253684*(92.2225568936364*(alpha[4]*favg[17]+favg[4]*alpha[17])+105.0*(alpha[6]*favg[13]+favg[6]*alpha[13])+84.0*(alpha[4]*favg[12]+favg[4]*alpha[12])+(93.91485505499116*alpha[8]+67.0820393249937*alpha[7]+105.0*alpha[0])*favg[11]+(93.91485505499116*favg[8]+67.0820393249937*favg[7]+105.0*favg[0])*alpha[11]+93.91485505499116*(alpha[5]*favg[10]+favg[5]*alpha[10])+105.0*(alpha[2]*favg[7]+favg[2]*alpha[7])+93.91485505499116*(alpha[1]*favg[4]+favg[1]*alpha[4]))-0.3535533905932737*(fr[19]-1.0*fl[19])*amax; 
  Ghat[12] = 0.001683587574253684*(92.2225568936364*(alpha[4]*favg[18]+favg[4]*alpha[18])+105.0*(alpha[5]*favg[14]+favg[5]*alpha[14])+(67.0820393249937*alpha[8]+93.91485505499116*alpha[7]+105.0*alpha[0])*favg[12]+(67.0820393249937*favg[8]+93.91485505499116*favg[7]+105.0*favg[0])*alpha[12]+84.0*(alpha[4]*favg[11]+favg[4]*alpha[11])+93.91485505499116*(alpha[6]*favg[10]+favg[6]*alpha[10])+105.0*(alpha[1]*favg[8]+favg[1]*alpha[8])+93.91485505499116*(alpha[2]*favg[4]+favg[2]*alpha[4]))-0.3535533905932737*(fr[20]-1.0*fl[20])*amax; 
  Ghat[13] = 0.001683587574253684*(92.2225568936364*(alpha[5]*favg[17]+favg[5]*alpha[17])+84.0*alpha[5]*favg[15]+(67.0820393249937*alpha[7]+105.0*alpha[0])*favg[13]+(93.91485505499116*favg[9]+67.0820393249937*favg[7]+105.0*favg[0])*alpha[13]+105.0*(alpha[6]*favg[11]+favg[6]*alpha[11])+93.91485505499116*(alpha[4]*favg[10]+favg[4]*alpha[10])+105.0*(alpha[3]*favg[7]+favg[3]*alpha[7])+93.91485505499116*(alpha[1]*favg[5]+favg[1]*alpha[5]))-0.3535533905932737*(fr[25]-1.0*fl[25])*amax; 
  Ghat[14] = 0.001683587574253684*(92.2225568936364*(alpha[6]*favg[18]+favg[6]*alpha[18])+84.0*alpha[6]*favg[16]+(67.0820393249937*alpha[8]+105.0*alpha[0])*favg[14]+(93.91485505499116*favg[9]+67.0820393249937*favg[8]+105.0*favg[0])*alpha[14]+105.0*(alpha[5]*favg[12]+favg[5]*alpha[12])+93.91485505499116*(alpha[4]*favg[10]+favg[4]*alpha[10])+105.0*(alpha[3]*favg[8]+favg[3]*alpha[8])+93.91485505499116*(alpha[2]*favg[6]+favg[2]*alpha[6]))-0.3535533905932737*(fr[26]-1.0*fl[26])*amax; 
  Ghat[15] = 0.001683587574253684*(92.2225568936364*alpha[5]*favg[19]+105.0*alpha[4]*favg[16]+(93.91485505499116*alpha[7]+105.0*alpha[0])*favg[15]+84.0*(alpha[5]*favg[13]+favg[5]*alpha[13])+93.91485505499116*(alpha[6]*favg[10]+favg[6]*alpha[10])+105.0*alpha[1]*favg[9]+93.91485505499116*(alpha[3]*favg[5]+favg[3]*alpha[5]))-0.3535533905932737*(fr[28]-1.0*fl[28])*amax; 
  Ghat[16] = 0.001683587574253684*(92.2225568936364*alpha[6]*favg[19]+(93.91485505499116*alpha[8]+105.0*alpha[0])*favg[16]+105.0*alpha[4]*favg[15]+84.0*(alpha[6]*favg[14]+favg[6]*alpha[14])+93.91485505499116*(alpha[5]*favg[10]+favg[5]*alpha[10])+105.0*alpha[2]*favg[9]+93.91485505499116*(alpha[3]*favg[6]+favg[3]*alpha[6]))-0.3535533905932737*(fr[29]-1.0*fl[29])*amax; 
  Ghat[17] = 0.001683587574253684*((62.60990336999411*alpha[7]+105.0*alpha[0])*favg[17]+(62.60990336999411*favg[7]+105.0*favg[0])*alpha[17]+92.2225568936364*(alpha[5]*favg[13]+favg[5]*alpha[13]+alpha[4]*favg[11]+favg[4]*alpha[11])+92.22255689363637*(alpha[1]*favg[7]+favg[1]*alpha[7]))-0.3535533905932737*(fr[31]-1.0*fl[31])*amax; 
  Ghat[18] = 0.001683587574253684*((62.60990336999411*alpha[8]+105.0*alpha[0])*favg[18]+(62.60990336999411*favg[8]+105.0*favg[0])*alpha[18]+92.2225568936364*(alpha[6]*favg[14]+favg[6]*alpha[14]+alpha[4]*favg[12]+favg[4]*alpha[12])+92.22255689363637*(alpha[2]*favg[8]+favg[2]*alpha[8]))-0.3535533905932737*(fr[32]-1.0*fl[32])*amax; 
  Ghat[19] = 0.005050762722761053*(35.0*alpha[0]*favg[19]+30.7408522978788*(alpha[6]*favg[16]+alpha[5]*favg[15])+30.7408522978788*alpha[3]*favg[9])-0.3535533905932737*(fr[34]-1.0*fl[34])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv10r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv10r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv10r; 
  outr[3] += -1.224744871391589*Ghat[0]*dv10r; 
  outr[4] += 0.7071067811865475*Ghat[3]*dv10r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv10r; 
  outr[6] += -1.224744871391589*Ghat[1]*dv10r; 
  outr[7] += -1.224744871391589*Ghat[2]*dv10r; 
  outr[8] += 0.7071067811865475*Ghat[5]*dv10r; 
  outr[9] += 0.7071067811865475*Ghat[6]*dv10r; 
  outr[10] += -1.224744871391589*Ghat[3]*dv10r; 
  outr[11] += 0.7071067811865475*Ghat[7]*dv10r; 
  outr[12] += 0.7071067811865475*Ghat[8]*dv10r; 
  outr[13] += 1.58113883008419*Ghat[0]*dv10r; 
  outr[14] += 0.7071067811865475*Ghat[9]*dv10r; 
  outr[15] += -1.224744871391589*Ghat[4]*dv10r; 
  outr[16] += 0.7071067811865475*Ghat[10]*dv10r; 
  outr[17] += -1.224744871391589*Ghat[5]*dv10r; 
  outr[18] += -1.224744871391589*Ghat[6]*dv10r; 
  outr[19] += 0.7071067811865475*Ghat[11]*dv10r; 
  outr[20] += 0.7071067811865475*Ghat[12]*dv10r; 
  outr[21] += -1.224744871391589*Ghat[7]*dv10r; 
  outr[22] += -1.224744871391589*Ghat[8]*dv10r; 
  outr[23] += 1.58113883008419*Ghat[1]*dv10r; 
  outr[24] += 1.58113883008419*Ghat[2]*dv10r; 
  outr[25] += 0.7071067811865475*Ghat[13]*dv10r; 
  outr[26] += 0.7071067811865475*Ghat[14]*dv10r; 
  outr[27] += 1.58113883008419*Ghat[3]*dv10r; 
  outr[28] += 0.7071067811865475*Ghat[15]*dv10r; 
  outr[29] += 0.7071067811865475*Ghat[16]*dv10r; 
  outr[30] += -1.224744871391589*Ghat[9]*dv10r; 
  outr[31] += 0.7071067811865475*Ghat[17]*dv10r; 
  outr[32] += 0.7071067811865475*Ghat[18]*dv10r; 
  outr[33] += -1.870828693386971*Ghat[0]*dv10r; 
  outr[34] += 0.7071067811865475*Ghat[19]*dv10r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv10l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv10l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv10l; 
  outl[3] += -1.224744871391589*Ghat[0]*dv10l; 
  outl[4] += -0.7071067811865475*Ghat[3]*dv10l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv10l; 
  outl[6] += -1.224744871391589*Ghat[1]*dv10l; 
  outl[7] += -1.224744871391589*Ghat[2]*dv10l; 
  outl[8] += -0.7071067811865475*Ghat[5]*dv10l; 
  outl[9] += -0.7071067811865475*Ghat[6]*dv10l; 
  outl[10] += -1.224744871391589*Ghat[3]*dv10l; 
  outl[11] += -0.7071067811865475*Ghat[7]*dv10l; 
  outl[12] += -0.7071067811865475*Ghat[8]*dv10l; 
  outl[13] += -1.58113883008419*Ghat[0]*dv10l; 
  outl[14] += -0.7071067811865475*Ghat[9]*dv10l; 
  outl[15] += -1.224744871391589*Ghat[4]*dv10l; 
  outl[16] += -0.7071067811865475*Ghat[10]*dv10l; 
  outl[17] += -1.224744871391589*Ghat[5]*dv10l; 
  outl[18] += -1.224744871391589*Ghat[6]*dv10l; 
  outl[19] += -0.7071067811865475*Ghat[11]*dv10l; 
  outl[20] += -0.7071067811865475*Ghat[12]*dv10l; 
  outl[21] += -1.224744871391589*Ghat[7]*dv10l; 
  outl[22] += -1.224744871391589*Ghat[8]*dv10l; 
  outl[23] += -1.58113883008419*Ghat[1]*dv10l; 
  outl[24] += -1.58113883008419*Ghat[2]*dv10l; 
  outl[25] += -0.7071067811865475*Ghat[13]*dv10l; 
  outl[26] += -0.7071067811865475*Ghat[14]*dv10l; 
  outl[27] += -1.58113883008419*Ghat[3]*dv10l; 
  outl[28] += -0.7071067811865475*Ghat[15]*dv10l; 
  outl[29] += -0.7071067811865475*Ghat[16]*dv10l; 
  outl[30] += -1.224744871391589*Ghat[9]*dv10l; 
  outl[31] += -0.7071067811865475*Ghat[17]*dv10l; 
  outl[32] += -0.7071067811865475*Ghat[18]*dv10l; 
  outl[33] += -1.870828693386971*Ghat[0]*dv10l; 
  outl[34] += -0.7071067811865475*Ghat[19]*dv10l; 

  return std::abs(amid); 
} 
__host__ __device__ double VlasovSurfElcMag2x2vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11l = 2/dxvl[3]; 
  double dv11r = 2/dxvr[3]; 
  const double *E1 = &EM[10]; 
  const double dv1 = dxvr[2], wv1 = wr[2]; 
  const double dv2 = dxvr[3], wv2 = wr[3]; 
  const double *B2 = &EM[50]; 

  double Ghat[20]; 
  double favg[20]; 
  double alpha[20]; 

  favg[0] = (-1.870828693386971*fr[34])+1.870828693386971*fl[34]+1.58113883008419*fr[14]+1.58113883008419*fl[14]-1.224744871391589*fr[4]+1.224744871391589*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0]; 
  favg[1] = 1.58113883008419*fr[28]+1.58113883008419*fl[28]-1.224744871391589*fr[8]+1.224744871391589*fl[8]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1]; 
  favg[2] = 1.58113883008419*fr[29]+1.58113883008419*fl[29]-1.224744871391589*fr[9]+1.224744871391589*fl[9]+0.7071067811865475*fr[2]+0.7071067811865475*fl[2]; 
  favg[3] = 1.58113883008419*fr[30]+1.58113883008419*fl[30]-1.224744871391589*fr[10]+1.224744871391589*fl[10]+0.7071067811865475*fr[3]+0.7071067811865475*fl[3]; 
  favg[4] = (-1.224744871391589*fr[16])+1.224744871391589*fl[16]+0.7071067811865475*fr[5]+0.7071067811865475*fl[5]; 
  favg[5] = (-1.224744871391589*fr[17])+1.224744871391589*fl[17]+0.7071067811865475*fr[6]+0.7071067811865475*fl[6]; 
  favg[6] = (-1.224744871391589*fr[18])+1.224744871391589*fl[18]+0.7071067811865475*fr[7]+0.7071067811865475*fl[7]; 
  favg[7] = (-1.224744871391589*fr[25])+1.224744871391589*fl[25]+0.7071067811865475*fr[11]+0.7071067811865475*fl[11]; 
  favg[8] = (-1.224744871391589*fr[26])+1.224744871391589*fl[26]+0.7071067811865475*fr[12]+0.7071067811865475*fl[12]; 
  favg[9] = (-1.224744871391589*fr[27])+1.224744871391589*fl[27]+0.7071067811865475*fr[13]+0.7071067811865475*fl[13]; 
  favg[10] = 0.7071067811865475*fr[15]+0.7071067811865475*fl[15]; 
  favg[11] = 0.7071067811865475*fr[19]+0.7071067811865475*fl[19]; 
  favg[12] = 0.7071067811865475*fr[20]+0.7071067811865475*fl[20]; 
  favg[13] = 0.7071067811865475*fr[21]+0.7071067811865475*fl[21]; 
  favg[14] = 0.7071067811865475*fr[22]+0.7071067811865475*fl[22]; 
  favg[15] = 0.7071067811865475*fr[23]+0.7071067811865475*fl[23]; 
  favg[16] = 0.7071067811865475*fr[24]+0.7071067811865475*fl[24]; 
  favg[17] = 0.7071067811865475*fr[31]+0.7071067811865475*fl[31]; 
  favg[18] = 0.7071067811865475*fr[32]+0.7071067811865475*fl[32]; 
  favg[19] = 0.7071067811865475*fr[33]+0.7071067811865475*fl[33]; 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = 1.414213562373095*E1[2]-1.414213562373095*B2[2]*wv1; 
  alpha[3] = -0.408248290463863*B2[0]*dv1; 
  alpha[4] = 1.414213562373095*E1[3]-1.414213562373095*B2[3]*wv1; 
  alpha[5] = -0.408248290463863*B2[1]*dv1; 
  alpha[6] = -0.408248290463863*B2[2]*dv1; 
  alpha[7] = 1.414213562373095*E1[4]-1.414213562373095*B2[4]*wv1; 
  alpha[8] = 1.414213562373095*E1[5]-1.414213562373095*B2[5]*wv1; 
  alpha[10] = -0.408248290463863*B2[3]*dv1; 
  alpha[11] = 1.414213562373095*E1[6]-1.414213562373095*B2[6]*wv1; 
  alpha[12] = 1.414213562373095*E1[7]-1.414213562373095*B2[7]*wv1; 
  alpha[13] = -0.408248290463863*B2[4]*dv1; 
  alpha[14] = -0.408248290463863*B2[5]*dv1; 
  alpha[17] = 1.414213562373095*E1[8]-1.414213562373095*B2[8]*wv1; 
  alpha[18] = 1.414213562373095*E1[9]-1.414213562373095*B2[9]*wv1; 

  const double amid = (-0.3952847075210473*alpha[8])-0.3952847075210473*alpha[7]+0.3535533905932737*alpha[0]; 

  Ghat[0] = 0.3535533905932737*(2.645751311064591*(fr[34]+fl[34])+2.23606797749979*(fl[14]-1.0*fr[14])+1.732050807568877*(fr[4]+fl[4])-1.0*fr[0]+fl[0])*amax+0.1767766952966368*(alpha[18]*favg[18]+alpha[17]*favg[17]+alpha[14]*favg[14]+alpha[13]*favg[13]+alpha[12]*favg[12]+alpha[11]*favg[11]+alpha[10]*favg[10]+alpha[8]*favg[8]+alpha[7]*favg[7]+alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[3]*favg[3]+alpha[2]*favg[2]+alpha[1]*favg[1]+alpha[0]*favg[0]); 
  Ghat[1] = 0.001683587574253684*(92.22255689363637*(alpha[7]*favg[17]+favg[7]*alpha[17])+93.91485505499116*(alpha[5]*favg[13]+favg[5]*alpha[13])+105.0*(alpha[8]*favg[12]+favg[8]*alpha[12])+93.91485505499116*(alpha[4]*favg[11]+favg[4]*alpha[11])+105.0*(alpha[6]*favg[10]+favg[6]*alpha[10])+93.91485505499116*(alpha[1]*favg[7]+favg[1]*alpha[7])+105.0*(alpha[3]*favg[5]+favg[3]*alpha[5]+alpha[2]*favg[4]+favg[2]*alpha[4]+alpha[0]*favg[1]+favg[0]*alpha[1]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[28]-1.0*(3.872983346207417*fl[28]+3.0*(fr[8]+fl[8])))+3.0*(fr[1]-1.0*fl[1]))*amax; 
  Ghat[2] = 0.001683587574253684*(92.22255689363637*(alpha[8]*favg[18]+favg[8]*alpha[18])+93.91485505499116*(alpha[6]*favg[14]+favg[6]*alpha[14]+alpha[4]*favg[12]+favg[4]*alpha[12])+105.0*(alpha[7]*favg[11]+favg[7]*alpha[11])+105.0*(alpha[5]*favg[10]+favg[5]*alpha[10])+93.91485505499116*(alpha[2]*favg[8]+favg[2]*alpha[8])+105.0*(alpha[3]*favg[6]+favg[3]*alpha[6]+alpha[1]*favg[4]+favg[1]*alpha[4]+alpha[0]*favg[2]+favg[0]*alpha[2]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[29]-1.0*(3.872983346207417*fl[29]+3.0*(fr[9]+fl[9])))+3.0*(fr[2]-1.0*fl[2]))*amax; 
  Ghat[3] = 0.01178511301977579*(13.41640786499874*(alpha[6]*favg[16]+alpha[5]*favg[15])+15.0*(alpha[8]*favg[14]+favg[8]*alpha[14]+alpha[7]*favg[13]+favg[7]*alpha[13])+15.0*(alpha[4]*favg[10]+favg[4]*alpha[10])+13.41640786499874*alpha[3]*favg[9]+15.0*(alpha[2]*favg[6]+favg[2]*alpha[6]+alpha[1]*favg[5]+favg[1]*alpha[5]+alpha[0]*favg[3]+favg[0]*alpha[3]))-0.1178511301977579*(1.732050807568877*(3.872983346207417*fr[30]-1.0*(3.872983346207417*fl[30]+3.0*(fr[10]+fl[10])))+3.0*(fr[3]-1.0*fl[3]))*amax; 
  Ghat[4] = 0.3535533905932737*(1.732050807568877*(fr[16]+fl[16])-1.0*fr[5]+fl[5])*amax+0.001683587574253684*(92.2225568936364*(alpha[12]*favg[18]+favg[12]*alpha[18]+alpha[11]*favg[17]+favg[11]*alpha[17])+93.91485505499116*(alpha[10]*favg[14]+favg[10]*alpha[14]+alpha[10]*favg[13]+favg[10]*alpha[13])+(84.0*alpha[11]+93.91485505499116*alpha[2])*favg[12]+(84.0*favg[11]+93.91485505499116*favg[2])*alpha[12]+93.91485505499116*(alpha[1]*favg[11]+favg[1]*alpha[11])+105.0*(alpha[3]*favg[10]+favg[3]*alpha[10])+93.91485505499116*(alpha[4]*favg[8]+favg[4]*alpha[8]+alpha[4]*favg[7]+favg[4]*alpha[7])+105.0*(alpha[5]*favg[6]+favg[5]*alpha[6]+alpha[0]*favg[4]+favg[0]*alpha[4]+alpha[1]*favg[2]+favg[1]*alpha[2])); 
  Ghat[5] = 0.3535533905932737*(1.732050807568877*(fr[17]+fl[17])-1.0*fr[6]+fl[6])*amax+0.001683587574253684*(92.2225568936364*(alpha[13]*favg[17]+favg[13]*alpha[17])+93.91485505499116*alpha[10]*favg[16]+(84.0*alpha[13]+93.91485505499116*alpha[3])*favg[15]+105.0*(alpha[12]*favg[14]+favg[12]*alpha[14])+93.91485505499116*(alpha[1]*favg[13]+favg[1]*alpha[13]+alpha[10]*favg[11]+favg[10]*alpha[11])+105.0*(alpha[2]*favg[10]+favg[2]*alpha[10])+93.91485505499116*(alpha[5]*(favg[9]+favg[7])+favg[5]*alpha[7])+105.0*(alpha[4]*favg[6]+favg[4]*alpha[6]+alpha[0]*favg[5]+favg[0]*alpha[5]+alpha[1]*favg[3]+favg[1]*alpha[3])); 
  Ghat[6] = 0.3535533905932737*(1.732050807568877*(fr[18]+fl[18])-1.0*fr[7]+fl[7])*amax+0.001683587574253684*(92.2225568936364*(alpha[14]*favg[18]+favg[14]*alpha[18])+(84.0*alpha[14]+93.91485505499116*alpha[3])*favg[16]+93.91485505499116*(alpha[10]*favg[15]+alpha[2]*favg[14]+favg[2]*alpha[14])+105.0*(alpha[11]*favg[13]+favg[11]*alpha[13])+93.91485505499116*(alpha[10]*favg[12]+favg[10]*alpha[12])+105.0*(alpha[1]*favg[10]+favg[1]*alpha[10])+93.91485505499116*(alpha[6]*(favg[9]+favg[8])+favg[6]*alpha[8])+105.0*(alpha[0]*favg[6]+favg[0]*alpha[6]+alpha[4]*favg[5]+favg[4]*alpha[5]+alpha[2]*favg[3]+favg[2]*alpha[3])); 
  Ghat[7] = 0.07071067811865475*(8.660254037844387*(fr[25]+fl[25])+5.0*(fl[11]-1.0*fr[11]))*amax+0.001683587574253684*((62.60990336999411*alpha[17]+92.22255689363637*alpha[1])*favg[17]+92.22255689363637*favg[1]*alpha[17]+(67.0820393249937*alpha[13]+105.0*alpha[3])*favg[13]+105.0*favg[3]*alpha[13]+93.91485505499116*alpha[12]*favg[12]+(67.0820393249937*alpha[11]+105.0*alpha[2])*favg[11]+105.0*favg[2]*alpha[11]+93.91485505499116*alpha[10]*favg[10]+(67.0820393249937*alpha[7]+105.0*alpha[0])*favg[7]+105.0*favg[0]*alpha[7]+93.91485505499116*(alpha[5]*favg[5]+alpha[4]*favg[4]+alpha[1]*favg[1])); 
  Ghat[8] = 0.07071067811865475*(8.660254037844387*(fr[26]+fl[26])+5.0*(fl[12]-1.0*fr[12]))*amax+0.001683587574253684*((62.60990336999411*alpha[18]+92.22255689363637*alpha[2])*favg[18]+92.22255689363637*favg[2]*alpha[18]+(67.0820393249937*alpha[14]+105.0*alpha[3])*favg[14]+105.0*favg[3]*alpha[14]+(67.0820393249937*alpha[12]+105.0*alpha[1])*favg[12]+105.0*favg[1]*alpha[12]+93.91485505499116*(alpha[11]*favg[11]+alpha[10]*favg[10])+(67.0820393249937*alpha[8]+105.0*alpha[0])*favg[8]+105.0*favg[0]*alpha[8]+93.91485505499116*(alpha[6]*favg[6]+alpha[4]*favg[4]+alpha[2]*favg[2])); 
  Ghat[9] = 0.07071067811865475*(8.660254037844387*(fr[27]+fl[27])+5.0*(fl[13]-1.0*fr[13]))*amax+0.001683587574253684*(92.22255689363637*alpha[3]*favg[19]+105.0*(alpha[2]*favg[16]+alpha[1]*favg[15])+93.91485505499116*(alpha[14]*favg[14]+alpha[13]*favg[13]+alpha[10]*favg[10])+105.0*alpha[0]*favg[9]+93.91485505499116*(alpha[6]*favg[6]+alpha[5]*favg[5]+alpha[3]*favg[3])); 
  Ghat[10] = 0.01178511301977579*(13.41640786499874*(alpha[5]*favg[16]+alpha[6]*favg[15]+alpha[4]*favg[14]+favg[4]*alpha[14]+alpha[4]*favg[13]+favg[4]*alpha[13]+alpha[6]*favg[12]+favg[6]*alpha[12]+alpha[5]*favg[11]+favg[5]*alpha[11])+(13.41640786499874*(alpha[8]+alpha[7])+15.0*alpha[0])*favg[10]+(13.41640786499874*(favg[9]+favg[8]+favg[7])+15.0*favg[0])*alpha[10]+15.0*(alpha[1]*favg[6]+favg[1]*alpha[6]+alpha[2]*favg[5]+favg[2]*alpha[5]+alpha[3]*favg[4]+favg[3]*alpha[4]))-0.3535533905932737*(fr[15]-1.0*fl[15])*amax; 
  Ghat[11] = 0.001683587574253684*(92.2225568936364*(alpha[4]*favg[17]+favg[4]*alpha[17])+105.0*(alpha[6]*favg[13]+favg[6]*alpha[13])+84.0*(alpha[4]*favg[12]+favg[4]*alpha[12])+(93.91485505499116*alpha[8]+67.0820393249937*alpha[7]+105.0*alpha[0])*favg[11]+(93.91485505499116*favg[8]+67.0820393249937*favg[7]+105.0*favg[0])*alpha[11]+93.91485505499116*(alpha[5]*favg[10]+favg[5]*alpha[10])+105.0*(alpha[2]*favg[7]+favg[2]*alpha[7])+93.91485505499116*(alpha[1]*favg[4]+favg[1]*alpha[4]))-0.3535533905932737*(fr[19]-1.0*fl[19])*amax; 
  Ghat[12] = 0.001683587574253684*(92.2225568936364*(alpha[4]*favg[18]+favg[4]*alpha[18])+105.0*(alpha[5]*favg[14]+favg[5]*alpha[14])+(67.0820393249937*alpha[8]+93.91485505499116*alpha[7]+105.0*alpha[0])*favg[12]+(67.0820393249937*favg[8]+93.91485505499116*favg[7]+105.0*favg[0])*alpha[12]+84.0*(alpha[4]*favg[11]+favg[4]*alpha[11])+93.91485505499116*(alpha[6]*favg[10]+favg[6]*alpha[10])+105.0*(alpha[1]*favg[8]+favg[1]*alpha[8])+93.91485505499116*(alpha[2]*favg[4]+favg[2]*alpha[4]))-0.3535533905932737*(fr[20]-1.0*fl[20])*amax; 
  Ghat[13] = 0.001683587574253684*(92.2225568936364*(alpha[5]*favg[17]+favg[5]*alpha[17])+84.0*alpha[5]*favg[15]+(67.0820393249937*alpha[7]+105.0*alpha[0])*favg[13]+(93.91485505499116*favg[9]+67.0820393249937*favg[7]+105.0*favg[0])*alpha[13]+105.0*(alpha[6]*favg[11]+favg[6]*alpha[11])+93.91485505499116*(alpha[4]*favg[10]+favg[4]*alpha[10])+105.0*(alpha[3]*favg[7]+favg[3]*alpha[7])+93.91485505499116*(alpha[1]*favg[5]+favg[1]*alpha[5]))-0.3535533905932737*(fr[21]-1.0*fl[21])*amax; 
  Ghat[14] = 0.001683587574253684*(92.2225568936364*(alpha[6]*favg[18]+favg[6]*alpha[18])+84.0*alpha[6]*favg[16]+(67.0820393249937*alpha[8]+105.0*alpha[0])*favg[14]+(93.91485505499116*favg[9]+67.0820393249937*favg[8]+105.0*favg[0])*alpha[14]+105.0*(alpha[5]*favg[12]+favg[5]*alpha[12])+93.91485505499116*(alpha[4]*favg[10]+favg[4]*alpha[10])+105.0*(alpha[3]*favg[8]+favg[3]*alpha[8])+93.91485505499116*(alpha[2]*favg[6]+favg[2]*alpha[6]))-0.3535533905932737*(fr[22]-1.0*fl[22])*amax; 
  Ghat[15] = 0.001683587574253684*(92.2225568936364*alpha[5]*favg[19]+105.0*alpha[4]*favg[16]+(93.91485505499116*alpha[7]+105.0*alpha[0])*favg[15]+84.0*(alpha[5]*favg[13]+favg[5]*alpha[13])+93.91485505499116*(alpha[6]*favg[10]+favg[6]*alpha[10])+105.0*alpha[1]*favg[9]+93.91485505499116*(alpha[3]*favg[5]+favg[3]*alpha[5]))-0.3535533905932737*(fr[23]-1.0*fl[23])*amax; 
  Ghat[16] = 0.001683587574253684*(92.2225568936364*alpha[6]*favg[19]+(93.91485505499116*alpha[8]+105.0*alpha[0])*favg[16]+105.0*alpha[4]*favg[15]+84.0*(alpha[6]*favg[14]+favg[6]*alpha[14])+93.91485505499116*(alpha[5]*favg[10]+favg[5]*alpha[10])+105.0*alpha[2]*favg[9]+93.91485505499116*(alpha[3]*favg[6]+favg[3]*alpha[6]))-0.3535533905932737*(fr[24]-1.0*fl[24])*amax; 
  Ghat[17] = 0.001683587574253684*((62.60990336999411*alpha[7]+105.0*alpha[0])*favg[17]+(62.60990336999411*favg[7]+105.0*favg[0])*alpha[17]+92.2225568936364*(alpha[5]*favg[13]+favg[5]*alpha[13]+alpha[4]*favg[11]+favg[4]*alpha[11])+92.22255689363637*(alpha[1]*favg[7]+favg[1]*alpha[7]))-0.3535533905932737*(fr[31]-1.0*fl[31])*amax; 
  Ghat[18] = 0.001683587574253684*((62.60990336999411*alpha[8]+105.0*alpha[0])*favg[18]+(62.60990336999411*favg[8]+105.0*favg[0])*alpha[18]+92.2225568936364*(alpha[6]*favg[14]+favg[6]*alpha[14]+alpha[4]*favg[12]+favg[4]*alpha[12])+92.22255689363637*(alpha[2]*favg[8]+favg[2]*alpha[8]))-0.3535533905932737*(fr[32]-1.0*fl[32])*amax; 
  Ghat[19] = 0.005050762722761053*(35.0*alpha[0]*favg[19]+30.7408522978788*(alpha[6]*favg[16]+alpha[5]*favg[15])+30.7408522978788*alpha[3]*favg[9])-0.3535533905932737*(fr[33]-1.0*fl[33])*amax; 

  outr[0] += 0.7071067811865475*Ghat[0]*dv11r; 
  outr[1] += 0.7071067811865475*Ghat[1]*dv11r; 
  outr[2] += 0.7071067811865475*Ghat[2]*dv11r; 
  outr[3] += 0.7071067811865475*Ghat[3]*dv11r; 
  outr[4] += -1.224744871391589*Ghat[0]*dv11r; 
  outr[5] += 0.7071067811865475*Ghat[4]*dv11r; 
  outr[6] += 0.7071067811865475*Ghat[5]*dv11r; 
  outr[7] += 0.7071067811865475*Ghat[6]*dv11r; 
  outr[8] += -1.224744871391589*Ghat[1]*dv11r; 
  outr[9] += -1.224744871391589*Ghat[2]*dv11r; 
  outr[10] += -1.224744871391589*Ghat[3]*dv11r; 
  outr[11] += 0.7071067811865475*Ghat[7]*dv11r; 
  outr[12] += 0.7071067811865475*Ghat[8]*dv11r; 
  outr[13] += 0.7071067811865475*Ghat[9]*dv11r; 
  outr[14] += 1.58113883008419*Ghat[0]*dv11r; 
  outr[15] += 0.7071067811865475*Ghat[10]*dv11r; 
  outr[16] += -1.224744871391589*Ghat[4]*dv11r; 
  outr[17] += -1.224744871391589*Ghat[5]*dv11r; 
  outr[18] += -1.224744871391589*Ghat[6]*dv11r; 
  outr[19] += 0.7071067811865475*Ghat[11]*dv11r; 
  outr[20] += 0.7071067811865475*Ghat[12]*dv11r; 
  outr[21] += 0.7071067811865475*Ghat[13]*dv11r; 
  outr[22] += 0.7071067811865475*Ghat[14]*dv11r; 
  outr[23] += 0.7071067811865475*Ghat[15]*dv11r; 
  outr[24] += 0.7071067811865475*Ghat[16]*dv11r; 
  outr[25] += -1.224744871391589*Ghat[7]*dv11r; 
  outr[26] += -1.224744871391589*Ghat[8]*dv11r; 
  outr[27] += -1.224744871391589*Ghat[9]*dv11r; 
  outr[28] += 1.58113883008419*Ghat[1]*dv11r; 
  outr[29] += 1.58113883008419*Ghat[2]*dv11r; 
  outr[30] += 1.58113883008419*Ghat[3]*dv11r; 
  outr[31] += 0.7071067811865475*Ghat[17]*dv11r; 
  outr[32] += 0.7071067811865475*Ghat[18]*dv11r; 
  outr[33] += 0.7071067811865475*Ghat[19]*dv11r; 
  outr[34] += -1.870828693386971*Ghat[0]*dv11r; 

  outl[0] += -0.7071067811865475*Ghat[0]*dv11l; 
  outl[1] += -0.7071067811865475*Ghat[1]*dv11l; 
  outl[2] += -0.7071067811865475*Ghat[2]*dv11l; 
  outl[3] += -0.7071067811865475*Ghat[3]*dv11l; 
  outl[4] += -1.224744871391589*Ghat[0]*dv11l; 
  outl[5] += -0.7071067811865475*Ghat[4]*dv11l; 
  outl[6] += -0.7071067811865475*Ghat[5]*dv11l; 
  outl[7] += -0.7071067811865475*Ghat[6]*dv11l; 
  outl[8] += -1.224744871391589*Ghat[1]*dv11l; 
  outl[9] += -1.224744871391589*Ghat[2]*dv11l; 
  outl[10] += -1.224744871391589*Ghat[3]*dv11l; 
  outl[11] += -0.7071067811865475*Ghat[7]*dv11l; 
  outl[12] += -0.7071067811865475*Ghat[8]*dv11l; 
  outl[13] += -0.7071067811865475*Ghat[9]*dv11l; 
  outl[14] += -1.58113883008419*Ghat[0]*dv11l; 
  outl[15] += -0.7071067811865475*Ghat[10]*dv11l; 
  outl[16] += -1.224744871391589*Ghat[4]*dv11l; 
  outl[17] += -1.224744871391589*Ghat[5]*dv11l; 
  outl[18] += -1.224744871391589*Ghat[6]*dv11l; 
  outl[19] += -0.7071067811865475*Ghat[11]*dv11l; 
  outl[20] += -0.7071067811865475*Ghat[12]*dv11l; 
  outl[21] += -0.7071067811865475*Ghat[13]*dv11l; 
  outl[22] += -0.7071067811865475*Ghat[14]*dv11l; 
  outl[23] += -0.7071067811865475*Ghat[15]*dv11l; 
  outl[24] += -0.7071067811865475*Ghat[16]*dv11l; 
  outl[25] += -1.224744871391589*Ghat[7]*dv11l; 
  outl[26] += -1.224744871391589*Ghat[8]*dv11l; 
  outl[27] += -1.224744871391589*Ghat[9]*dv11l; 
  outl[28] += -1.58113883008419*Ghat[1]*dv11l; 
  outl[29] += -1.58113883008419*Ghat[2]*dv11l; 
  outl[30] += -1.58113883008419*Ghat[3]*dv11l; 
  outl[31] += -0.7071067811865475*Ghat[17]*dv11l; 
  outl[32] += -0.7071067811865475*Ghat[18]*dv11l; 
  outl[33] += -0.7071067811865475*Ghat[19]*dv11l; 
  outl[34] += -1.870828693386971*Ghat[0]*dv11l; 

  return std::abs(amid); 
} 