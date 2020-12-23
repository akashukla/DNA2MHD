#include <VlasovModDecl.h> 
__host__ __device__ void VlasovSurfStream2x3vSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double incr[32]; 

  if (wr[2]>0) { 
  incr[0] = dxvl[2]*(0.25*fl[7]+0.1443375672974065*fl[3])+(0.8660254037844386*fl[1]+0.5*fl[0])*wl[2]; 
  incr[1] = dxvl[2]*((-0.4330127018922193*fl[7])-0.25*fl[3])+((-1.5*fl[1])-0.8660254037844386*fl[0])*wl[2]; 
  incr[2] = dxvl[2]*(0.25*fl[16]+0.1443375672974065*fl[8])+wl[2]*(0.8660254037844386*fl[6]+0.5*fl[2]); 
  incr[3] = wl[2]*(0.8660254037844386*fl[7]+0.5*fl[3])+(0.25*fl[1]+0.1443375672974065*fl[0])*dxvl[2]; 
  incr[4] = dxvl[2]*(0.25*fl[18]+0.1443375672974065*fl[11])+wl[2]*(0.8660254037844386*fl[9]+0.5*fl[4]); 
  incr[5] = dxvl[2]*(0.25*fl[21]+0.1443375672974065*fl[14])+wl[2]*(0.8660254037844386*fl[12]+0.5*fl[5]); 
  incr[6] = dxvl[2]*((-0.4330127018922193*fl[16])-0.25*fl[8])+wl[2]*((-1.5*fl[6])-0.8660254037844386*fl[2]); 
  incr[7] = wl[2]*((-1.5*fl[7])-0.8660254037844386*fl[3])+((-0.4330127018922193*fl[1])-0.25*fl[0])*dxvl[2]; 
  incr[8] = wl[2]*(0.8660254037844386*fl[16]+0.5*fl[8])+dxvl[2]*(0.25*fl[6]+0.1443375672974065*fl[2]); 
  incr[9] = dxvl[2]*((-0.4330127018922193*fl[18])-0.25*fl[11])+wl[2]*((-1.5*fl[9])-0.8660254037844386*fl[4]); 
  incr[10] = dxvl[2]*(0.25*fl[26]+0.1443375672974065*fl[19])+wl[2]*(0.8660254037844386*fl[17]+0.5*fl[10]); 
  incr[11] = wl[2]*(0.8660254037844386*fl[18]+0.5*fl[11])+dxvl[2]*(0.25*fl[9]+0.1443375672974065*fl[4]); 
  incr[12] = dxvl[2]*((-0.4330127018922193*fl[21])-0.25*fl[14])+wl[2]*((-1.5*fl[12])-0.8660254037844386*fl[5]); 
  incr[13] = dxvl[2]*(0.25*fl[27]+0.1443375672974065*fl[22])+wl[2]*(0.8660254037844386*fl[20]+0.5*fl[13]); 
  incr[14] = wl[2]*(0.8660254037844386*fl[21]+0.5*fl[14])+dxvl[2]*(0.25*fl[12]+0.1443375672974065*fl[5]); 
  incr[15] = dxvl[2]*(0.25*fl[29]+0.1443375672974065*fl[25])+wl[2]*(0.8660254037844386*fl[23]+0.5*fl[15]); 
  incr[16] = wl[2]*((-1.5*fl[16])-0.8660254037844386*fl[8])+dxvl[2]*((-0.4330127018922193*fl[6])-0.25*fl[2]); 
  incr[17] = dxvl[2]*((-0.4330127018922193*fl[26])-0.25*fl[19])+wl[2]*((-1.5*fl[17])-0.8660254037844386*fl[10]); 
  incr[18] = wl[2]*((-1.5*fl[18])-0.8660254037844386*fl[11])+dxvl[2]*((-0.4330127018922193*fl[9])-0.25*fl[4]); 
  incr[19] = wl[2]*(0.8660254037844386*fl[26]+0.5*fl[19])+dxvl[2]*(0.25*fl[17]+0.1443375672974065*fl[10]); 
  incr[20] = dxvl[2]*((-0.4330127018922193*fl[27])-0.25*fl[22])+wl[2]*((-1.5*fl[20])-0.8660254037844386*fl[13]); 
  incr[21] = wl[2]*((-1.5*fl[21])-0.8660254037844386*fl[14])+dxvl[2]*((-0.4330127018922193*fl[12])-0.25*fl[5]); 
  incr[22] = wl[2]*(0.8660254037844386*fl[27]+0.5*fl[22])+dxvl[2]*(0.25*fl[20]+0.1443375672974065*fl[13]); 
  incr[23] = dxvl[2]*((-0.4330127018922193*fl[29])-0.25*fl[25])+wl[2]*((-1.5*fl[23])-0.8660254037844386*fl[15]); 
  incr[24] = dxvl[2]*(0.25*fl[31]+0.1443375672974065*fl[30])+wl[2]*(0.8660254037844386*fl[28]+0.5*fl[24]); 
  incr[25] = wl[2]*(0.8660254037844386*fl[29]+0.5*fl[25])+dxvl[2]*(0.25*fl[23]+0.1443375672974065*fl[15]); 
  incr[26] = wl[2]*((-1.5*fl[26])-0.8660254037844386*fl[19])+dxvl[2]*((-0.4330127018922193*fl[17])-0.25*fl[10]); 
  incr[27] = wl[2]*((-1.5*fl[27])-0.8660254037844386*fl[22])+dxvl[2]*((-0.4330127018922193*fl[20])-0.25*fl[13]); 
  incr[28] = dxvl[2]*((-0.4330127018922193*fl[31])-0.25*fl[30])+wl[2]*((-1.5*fl[28])-0.8660254037844386*fl[24]); 
  incr[29] = wl[2]*((-1.5*fl[29])-0.8660254037844386*fl[25])+dxvl[2]*((-0.4330127018922193*fl[23])-0.25*fl[15]); 
  incr[30] = wl[2]*(0.8660254037844386*fl[31]+0.5*fl[30])+dxvl[2]*(0.25*fl[28]+0.1443375672974065*fl[24]); 
  incr[31] = wl[2]*((-1.5*fl[31])-0.8660254037844386*fl[30])+dxvl[2]*((-0.4330127018922193*fl[28])-0.25*fl[24]); 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 
  outr[8] += incr[8]*rdxr2; 
  outr[9] += incr[9]*rdxr2; 
  outr[10] += incr[10]*rdxr2; 
  outr[11] += incr[11]*rdxr2; 
  outr[12] += incr[12]*rdxr2; 
  outr[13] += incr[13]*rdxr2; 
  outr[14] += incr[14]*rdxr2; 
  outr[15] += incr[15]*rdxr2; 
  outr[16] += incr[16]*rdxr2; 
  outr[17] += incr[17]*rdxr2; 
  outr[18] += incr[18]*rdxr2; 
  outr[19] += incr[19]*rdxr2; 
  outr[20] += incr[20]*rdxr2; 
  outr[21] += incr[21]*rdxr2; 
  outr[22] += incr[22]*rdxr2; 
  outr[23] += incr[23]*rdxr2; 
  outr[24] += incr[24]*rdxr2; 
  outr[25] += incr[25]*rdxr2; 
  outr[26] += incr[26]*rdxr2; 
  outr[27] += incr[27]*rdxr2; 
  outr[28] += incr[28]*rdxr2; 
  outr[29] += incr[29]*rdxr2; 
  outr[30] += incr[30]*rdxr2; 
  outr[31] += incr[31]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  outl[8] += -1.0*incr[8]*rdxl2; 
  outl[9] += incr[9]*rdxl2; 
  outl[10] += -1.0*incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += incr[12]*rdxl2; 
  outl[13] += -1.0*incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += -1.0*incr[15]*rdxl2; 
  outl[16] += incr[16]*rdxl2; 
  outl[17] += incr[17]*rdxl2; 
  outl[18] += incr[18]*rdxl2; 
  outl[19] += -1.0*incr[19]*rdxl2; 
  outl[20] += incr[20]*rdxl2; 
  outl[21] += incr[21]*rdxl2; 
  outl[22] += -1.0*incr[22]*rdxl2; 
  outl[23] += incr[23]*rdxl2; 
  outl[24] += -1.0*incr[24]*rdxl2; 
  outl[25] += -1.0*incr[25]*rdxl2; 
  outl[26] += incr[26]*rdxl2; 
  outl[27] += incr[27]*rdxl2; 
  outl[28] += incr[28]*rdxl2; 
  outl[29] += incr[29]*rdxl2; 
  outl[30] += -1.0*incr[30]*rdxl2; 
  outl[31] += incr[31]*rdxl2; 
  } else { 
  incr[0] = dxvr[2]*(0.1443375672974065*fr[3]-0.25*fr[7])+(0.5*fr[0]-0.8660254037844386*fr[1])*wr[2]; 
  incr[1] = dxvr[2]*(0.4330127018922193*fr[7]-0.25*fr[3])+(1.5*fr[1]-0.8660254037844386*fr[0])*wr[2]; 
  incr[2] = dxvr[2]*(0.1443375672974065*fr[8]-0.25*fr[16])+wr[2]*(0.5*fr[2]-0.8660254037844386*fr[6]); 
  incr[3] = wr[2]*(0.5*fr[3]-0.8660254037844386*fr[7])+(0.1443375672974065*fr[0]-0.25*fr[1])*dxvr[2]; 
  incr[4] = dxvr[2]*(0.1443375672974065*fr[11]-0.25*fr[18])+wr[2]*(0.5*fr[4]-0.8660254037844386*fr[9]); 
  incr[5] = dxvr[2]*(0.1443375672974065*fr[14]-0.25*fr[21])+wr[2]*(0.5*fr[5]-0.8660254037844386*fr[12]); 
  incr[6] = dxvr[2]*(0.4330127018922193*fr[16]-0.25*fr[8])+wr[2]*(1.5*fr[6]-0.8660254037844386*fr[2]); 
  incr[7] = wr[2]*(1.5*fr[7]-0.8660254037844386*fr[3])+(0.4330127018922193*fr[1]-0.25*fr[0])*dxvr[2]; 
  incr[8] = wr[2]*(0.5*fr[8]-0.8660254037844386*fr[16])+dxvr[2]*(0.1443375672974065*fr[2]-0.25*fr[6]); 
  incr[9] = dxvr[2]*(0.4330127018922193*fr[18]-0.25*fr[11])+wr[2]*(1.5*fr[9]-0.8660254037844386*fr[4]); 
  incr[10] = dxvr[2]*(0.1443375672974065*fr[19]-0.25*fr[26])+wr[2]*(0.5*fr[10]-0.8660254037844386*fr[17]); 
  incr[11] = wr[2]*(0.5*fr[11]-0.8660254037844386*fr[18])+dxvr[2]*(0.1443375672974065*fr[4]-0.25*fr[9]); 
  incr[12] = dxvr[2]*(0.4330127018922193*fr[21]-0.25*fr[14])+wr[2]*(1.5*fr[12]-0.8660254037844386*fr[5]); 
  incr[13] = dxvr[2]*(0.1443375672974065*fr[22]-0.25*fr[27])+wr[2]*(0.5*fr[13]-0.8660254037844386*fr[20]); 
  incr[14] = wr[2]*(0.5*fr[14]-0.8660254037844386*fr[21])+dxvr[2]*(0.1443375672974065*fr[5]-0.25*fr[12]); 
  incr[15] = dxvr[2]*(0.1443375672974065*fr[25]-0.25*fr[29])+wr[2]*(0.5*fr[15]-0.8660254037844386*fr[23]); 
  incr[16] = wr[2]*(1.5*fr[16]-0.8660254037844386*fr[8])+dxvr[2]*(0.4330127018922193*fr[6]-0.25*fr[2]); 
  incr[17] = dxvr[2]*(0.4330127018922193*fr[26]-0.25*fr[19])+wr[2]*(1.5*fr[17]-0.8660254037844386*fr[10]); 
  incr[18] = wr[2]*(1.5*fr[18]-0.8660254037844386*fr[11])+dxvr[2]*(0.4330127018922193*fr[9]-0.25*fr[4]); 
  incr[19] = wr[2]*(0.5*fr[19]-0.8660254037844386*fr[26])+dxvr[2]*(0.1443375672974065*fr[10]-0.25*fr[17]); 
  incr[20] = dxvr[2]*(0.4330127018922193*fr[27]-0.25*fr[22])+wr[2]*(1.5*fr[20]-0.8660254037844386*fr[13]); 
  incr[21] = wr[2]*(1.5*fr[21]-0.8660254037844386*fr[14])+dxvr[2]*(0.4330127018922193*fr[12]-0.25*fr[5]); 
  incr[22] = wr[2]*(0.5*fr[22]-0.8660254037844386*fr[27])+dxvr[2]*(0.1443375672974065*fr[13]-0.25*fr[20]); 
  incr[23] = dxvr[2]*(0.4330127018922193*fr[29]-0.25*fr[25])+wr[2]*(1.5*fr[23]-0.8660254037844386*fr[15]); 
  incr[24] = dxvr[2]*(0.1443375672974065*fr[30]-0.25*fr[31])+wr[2]*(0.5*fr[24]-0.8660254037844386*fr[28]); 
  incr[25] = wr[2]*(0.5*fr[25]-0.8660254037844386*fr[29])+dxvr[2]*(0.1443375672974065*fr[15]-0.25*fr[23]); 
  incr[26] = wr[2]*(1.5*fr[26]-0.8660254037844386*fr[19])+dxvr[2]*(0.4330127018922193*fr[17]-0.25*fr[10]); 
  incr[27] = wr[2]*(1.5*fr[27]-0.8660254037844386*fr[22])+dxvr[2]*(0.4330127018922193*fr[20]-0.25*fr[13]); 
  incr[28] = dxvr[2]*(0.4330127018922193*fr[31]-0.25*fr[30])+wr[2]*(1.5*fr[28]-0.8660254037844386*fr[24]); 
  incr[29] = wr[2]*(1.5*fr[29]-0.8660254037844386*fr[25])+dxvr[2]*(0.4330127018922193*fr[23]-0.25*fr[15]); 
  incr[30] = wr[2]*(0.5*fr[30]-0.8660254037844386*fr[31])+dxvr[2]*(0.1443375672974065*fr[24]-0.25*fr[28]); 
  incr[31] = wr[2]*(1.5*fr[31]-0.8660254037844386*fr[30])+dxvr[2]*(0.4330127018922193*fr[28]-0.25*fr[24]); 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 
  outr[8] += incr[8]*rdxr2; 
  outr[9] += incr[9]*rdxr2; 
  outr[10] += incr[10]*rdxr2; 
  outr[11] += incr[11]*rdxr2; 
  outr[12] += incr[12]*rdxr2; 
  outr[13] += incr[13]*rdxr2; 
  outr[14] += incr[14]*rdxr2; 
  outr[15] += incr[15]*rdxr2; 
  outr[16] += incr[16]*rdxr2; 
  outr[17] += incr[17]*rdxr2; 
  outr[18] += incr[18]*rdxr2; 
  outr[19] += incr[19]*rdxr2; 
  outr[20] += incr[20]*rdxr2; 
  outr[21] += incr[21]*rdxr2; 
  outr[22] += incr[22]*rdxr2; 
  outr[23] += incr[23]*rdxr2; 
  outr[24] += incr[24]*rdxr2; 
  outr[25] += incr[25]*rdxr2; 
  outr[26] += incr[26]*rdxr2; 
  outr[27] += incr[27]*rdxr2; 
  outr[28] += incr[28]*rdxr2; 
  outr[29] += incr[29]*rdxr2; 
  outr[30] += incr[30]*rdxr2; 
  outr[31] += incr[31]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  outl[8] += -1.0*incr[8]*rdxl2; 
  outl[9] += incr[9]*rdxl2; 
  outl[10] += -1.0*incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += incr[12]*rdxl2; 
  outl[13] += -1.0*incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += -1.0*incr[15]*rdxl2; 
  outl[16] += incr[16]*rdxl2; 
  outl[17] += incr[17]*rdxl2; 
  outl[18] += incr[18]*rdxl2; 
  outl[19] += -1.0*incr[19]*rdxl2; 
  outl[20] += incr[20]*rdxl2; 
  outl[21] += incr[21]*rdxl2; 
  outl[22] += -1.0*incr[22]*rdxl2; 
  outl[23] += incr[23]*rdxl2; 
  outl[24] += -1.0*incr[24]*rdxl2; 
  outl[25] += -1.0*incr[25]*rdxl2; 
  outl[26] += incr[26]*rdxl2; 
  outl[27] += incr[27]*rdxl2; 
  outl[28] += incr[28]*rdxl2; 
  outl[29] += incr[29]*rdxl2; 
  outl[30] += -1.0*incr[30]*rdxl2; 
  outl[31] += incr[31]*rdxl2; 
  } 
} 
__host__ __device__ void VlasovSurfStream2x3vSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[1]; 
  double rdxr2 = 2.0/dxvr[1]; 

  double incr[32]; 

  if (wr[3]>0) { 
  incr[0] = dxvl[3]*(0.25*fl[10]+0.1443375672974065*fl[4])+(0.8660254037844386*fl[2]+0.5*fl[0])*wl[3]; 
  incr[1] = dxvl[3]*(0.25*fl[17]+0.1443375672974065*fl[9])+wl[3]*(0.8660254037844386*fl[6]+0.5*fl[1]); 
  incr[2] = dxvl[3]*((-0.4330127018922193*fl[10])-0.25*fl[4])+((-1.5*fl[2])-0.8660254037844386*fl[0])*wl[3]; 
  incr[3] = dxvl[3]*(0.25*fl[19]+0.1443375672974065*fl[11])+wl[3]*(0.8660254037844386*fl[8]+0.5*fl[3]); 
  incr[4] = wl[3]*(0.8660254037844386*fl[10]+0.5*fl[4])+(0.25*fl[2]+0.1443375672974065*fl[0])*dxvl[3]; 
  incr[5] = dxvl[3]*(0.25*fl[24]+0.1443375672974065*fl[15])+wl[3]*(0.8660254037844386*fl[13]+0.5*fl[5]); 
  incr[6] = dxvl[3]*((-0.4330127018922193*fl[17])-0.25*fl[9])+wl[3]*((-1.5*fl[6])-0.8660254037844386*fl[1]); 
  incr[7] = dxvl[3]*(0.25*fl[26]+0.1443375672974065*fl[18])+wl[3]*(0.8660254037844386*fl[16]+0.5*fl[7]); 
  incr[8] = dxvl[3]*((-0.4330127018922193*fl[19])-0.25*fl[11])+wl[3]*((-1.5*fl[8])-0.8660254037844386*fl[3]); 
  incr[9] = wl[3]*(0.8660254037844386*fl[17]+0.5*fl[9])+dxvl[3]*(0.25*fl[6]+0.1443375672974065*fl[1]); 
  incr[10] = wl[3]*((-1.5*fl[10])-0.8660254037844386*fl[4])+((-0.4330127018922193*fl[2])-0.25*fl[0])*dxvl[3]; 
  incr[11] = wl[3]*(0.8660254037844386*fl[19]+0.5*fl[11])+dxvl[3]*(0.25*fl[8]+0.1443375672974065*fl[3]); 
  incr[12] = dxvl[3]*(0.25*fl[28]+0.1443375672974065*fl[23])+wl[3]*(0.8660254037844386*fl[20]+0.5*fl[12]); 
  incr[13] = dxvl[3]*((-0.4330127018922193*fl[24])-0.25*fl[15])+wl[3]*((-1.5*fl[13])-0.8660254037844386*fl[5]); 
  incr[14] = dxvl[3]*(0.25*fl[30]+0.1443375672974065*fl[25])+wl[3]*(0.8660254037844386*fl[22]+0.5*fl[14]); 
  incr[15] = wl[3]*(0.8660254037844386*fl[24]+0.5*fl[15])+dxvl[3]*(0.25*fl[13]+0.1443375672974065*fl[5]); 
  incr[16] = dxvl[3]*((-0.4330127018922193*fl[26])-0.25*fl[18])+wl[3]*((-1.5*fl[16])-0.8660254037844386*fl[7]); 
  incr[17] = wl[3]*((-1.5*fl[17])-0.8660254037844386*fl[9])+dxvl[3]*((-0.4330127018922193*fl[6])-0.25*fl[1]); 
  incr[18] = wl[3]*(0.8660254037844386*fl[26]+0.5*fl[18])+dxvl[3]*(0.25*fl[16]+0.1443375672974065*fl[7]); 
  incr[19] = wl[3]*((-1.5*fl[19])-0.8660254037844386*fl[11])+dxvl[3]*((-0.4330127018922193*fl[8])-0.25*fl[3]); 
  incr[20] = dxvl[3]*((-0.4330127018922193*fl[28])-0.25*fl[23])+wl[3]*((-1.5*fl[20])-0.8660254037844386*fl[12]); 
  incr[21] = dxvl[3]*(0.25*fl[31]+0.1443375672974065*fl[29])+wl[3]*(0.8660254037844386*fl[27]+0.5*fl[21]); 
  incr[22] = dxvl[3]*((-0.4330127018922193*fl[30])-0.25*fl[25])+wl[3]*((-1.5*fl[22])-0.8660254037844386*fl[14]); 
  incr[23] = wl[3]*(0.8660254037844386*fl[28]+0.5*fl[23])+dxvl[3]*(0.25*fl[20]+0.1443375672974065*fl[12]); 
  incr[24] = wl[3]*((-1.5*fl[24])-0.8660254037844386*fl[15])+dxvl[3]*((-0.4330127018922193*fl[13])-0.25*fl[5]); 
  incr[25] = wl[3]*(0.8660254037844386*fl[30]+0.5*fl[25])+dxvl[3]*(0.25*fl[22]+0.1443375672974065*fl[14]); 
  incr[26] = wl[3]*((-1.5*fl[26])-0.8660254037844386*fl[18])+dxvl[3]*((-0.4330127018922193*fl[16])-0.25*fl[7]); 
  incr[27] = dxvl[3]*((-0.4330127018922193*fl[31])-0.25*fl[29])+wl[3]*((-1.5*fl[27])-0.8660254037844386*fl[21]); 
  incr[28] = wl[3]*((-1.5*fl[28])-0.8660254037844386*fl[23])+dxvl[3]*((-0.4330127018922193*fl[20])-0.25*fl[12]); 
  incr[29] = wl[3]*(0.8660254037844386*fl[31]+0.5*fl[29])+dxvl[3]*(0.25*fl[27]+0.1443375672974065*fl[21]); 
  incr[30] = wl[3]*((-1.5*fl[30])-0.8660254037844386*fl[25])+dxvl[3]*((-0.4330127018922193*fl[22])-0.25*fl[14]); 
  incr[31] = wl[3]*((-1.5*fl[31])-0.8660254037844386*fl[29])+dxvl[3]*((-0.4330127018922193*fl[27])-0.25*fl[21]); 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 
  outr[8] += incr[8]*rdxr2; 
  outr[9] += incr[9]*rdxr2; 
  outr[10] += incr[10]*rdxr2; 
  outr[11] += incr[11]*rdxr2; 
  outr[12] += incr[12]*rdxr2; 
  outr[13] += incr[13]*rdxr2; 
  outr[14] += incr[14]*rdxr2; 
  outr[15] += incr[15]*rdxr2; 
  outr[16] += incr[16]*rdxr2; 
  outr[17] += incr[17]*rdxr2; 
  outr[18] += incr[18]*rdxr2; 
  outr[19] += incr[19]*rdxr2; 
  outr[20] += incr[20]*rdxr2; 
  outr[21] += incr[21]*rdxr2; 
  outr[22] += incr[22]*rdxr2; 
  outr[23] += incr[23]*rdxr2; 
  outr[24] += incr[24]*rdxr2; 
  outr[25] += incr[25]*rdxr2; 
  outr[26] += incr[26]*rdxr2; 
  outr[27] += incr[27]*rdxr2; 
  outr[28] += incr[28]*rdxr2; 
  outr[29] += incr[29]*rdxr2; 
  outr[30] += incr[30]*rdxr2; 
  outr[31] += incr[31]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += -1.0*incr[1]*rdxl2; 
  outl[2] += incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += -1.0*incr[7]*rdxl2; 
  outl[8] += incr[8]*rdxl2; 
  outl[9] += -1.0*incr[9]*rdxl2; 
  outl[10] += incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += -1.0*incr[12]*rdxl2; 
  outl[13] += incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += -1.0*incr[15]*rdxl2; 
  outl[16] += incr[16]*rdxl2; 
  outl[17] += incr[17]*rdxl2; 
  outl[18] += -1.0*incr[18]*rdxl2; 
  outl[19] += incr[19]*rdxl2; 
  outl[20] += incr[20]*rdxl2; 
  outl[21] += -1.0*incr[21]*rdxl2; 
  outl[22] += incr[22]*rdxl2; 
  outl[23] += -1.0*incr[23]*rdxl2; 
  outl[24] += incr[24]*rdxl2; 
  outl[25] += -1.0*incr[25]*rdxl2; 
  outl[26] += incr[26]*rdxl2; 
  outl[27] += incr[27]*rdxl2; 
  outl[28] += incr[28]*rdxl2; 
  outl[29] += -1.0*incr[29]*rdxl2; 
  outl[30] += incr[30]*rdxl2; 
  outl[31] += incr[31]*rdxl2; 
  } else { 
  incr[0] = dxvr[3]*(0.1443375672974065*fr[4]-0.25*fr[10])+(0.5*fr[0]-0.8660254037844386*fr[2])*wr[3]; 
  incr[1] = dxvr[3]*(0.1443375672974065*fr[9]-0.25*fr[17])+wr[3]*(0.5*fr[1]-0.8660254037844386*fr[6]); 
  incr[2] = dxvr[3]*(0.4330127018922193*fr[10]-0.25*fr[4])+(1.5*fr[2]-0.8660254037844386*fr[0])*wr[3]; 
  incr[3] = dxvr[3]*(0.1443375672974065*fr[11]-0.25*fr[19])+wr[3]*(0.5*fr[3]-0.8660254037844386*fr[8]); 
  incr[4] = wr[3]*(0.5*fr[4]-0.8660254037844386*fr[10])+(0.1443375672974065*fr[0]-0.25*fr[2])*dxvr[3]; 
  incr[5] = dxvr[3]*(0.1443375672974065*fr[15]-0.25*fr[24])+wr[3]*(0.5*fr[5]-0.8660254037844386*fr[13]); 
  incr[6] = dxvr[3]*(0.4330127018922193*fr[17]-0.25*fr[9])+wr[3]*(1.5*fr[6]-0.8660254037844386*fr[1]); 
  incr[7] = dxvr[3]*(0.1443375672974065*fr[18]-0.25*fr[26])+wr[3]*(0.5*fr[7]-0.8660254037844386*fr[16]); 
  incr[8] = dxvr[3]*(0.4330127018922193*fr[19]-0.25*fr[11])+wr[3]*(1.5*fr[8]-0.8660254037844386*fr[3]); 
  incr[9] = wr[3]*(0.5*fr[9]-0.8660254037844386*fr[17])+dxvr[3]*(0.1443375672974065*fr[1]-0.25*fr[6]); 
  incr[10] = wr[3]*(1.5*fr[10]-0.8660254037844386*fr[4])+(0.4330127018922193*fr[2]-0.25*fr[0])*dxvr[3]; 
  incr[11] = wr[3]*(0.5*fr[11]-0.8660254037844386*fr[19])+dxvr[3]*(0.1443375672974065*fr[3]-0.25*fr[8]); 
  incr[12] = dxvr[3]*(0.1443375672974065*fr[23]-0.25*fr[28])+wr[3]*(0.5*fr[12]-0.8660254037844386*fr[20]); 
  incr[13] = dxvr[3]*(0.4330127018922193*fr[24]-0.25*fr[15])+wr[3]*(1.5*fr[13]-0.8660254037844386*fr[5]); 
  incr[14] = dxvr[3]*(0.1443375672974065*fr[25]-0.25*fr[30])+wr[3]*(0.5*fr[14]-0.8660254037844386*fr[22]); 
  incr[15] = wr[3]*(0.5*fr[15]-0.8660254037844386*fr[24])+dxvr[3]*(0.1443375672974065*fr[5]-0.25*fr[13]); 
  incr[16] = dxvr[3]*(0.4330127018922193*fr[26]-0.25*fr[18])+wr[3]*(1.5*fr[16]-0.8660254037844386*fr[7]); 
  incr[17] = wr[3]*(1.5*fr[17]-0.8660254037844386*fr[9])+dxvr[3]*(0.4330127018922193*fr[6]-0.25*fr[1]); 
  incr[18] = wr[3]*(0.5*fr[18]-0.8660254037844386*fr[26])+dxvr[3]*(0.1443375672974065*fr[7]-0.25*fr[16]); 
  incr[19] = wr[3]*(1.5*fr[19]-0.8660254037844386*fr[11])+dxvr[3]*(0.4330127018922193*fr[8]-0.25*fr[3]); 
  incr[20] = dxvr[3]*(0.4330127018922193*fr[28]-0.25*fr[23])+wr[3]*(1.5*fr[20]-0.8660254037844386*fr[12]); 
  incr[21] = dxvr[3]*(0.1443375672974065*fr[29]-0.25*fr[31])+wr[3]*(0.5*fr[21]-0.8660254037844386*fr[27]); 
  incr[22] = dxvr[3]*(0.4330127018922193*fr[30]-0.25*fr[25])+wr[3]*(1.5*fr[22]-0.8660254037844386*fr[14]); 
  incr[23] = wr[3]*(0.5*fr[23]-0.8660254037844386*fr[28])+dxvr[3]*(0.1443375672974065*fr[12]-0.25*fr[20]); 
  incr[24] = wr[3]*(1.5*fr[24]-0.8660254037844386*fr[15])+dxvr[3]*(0.4330127018922193*fr[13]-0.25*fr[5]); 
  incr[25] = wr[3]*(0.5*fr[25]-0.8660254037844386*fr[30])+dxvr[3]*(0.1443375672974065*fr[14]-0.25*fr[22]); 
  incr[26] = wr[3]*(1.5*fr[26]-0.8660254037844386*fr[18])+dxvr[3]*(0.4330127018922193*fr[16]-0.25*fr[7]); 
  incr[27] = dxvr[3]*(0.4330127018922193*fr[31]-0.25*fr[29])+wr[3]*(1.5*fr[27]-0.8660254037844386*fr[21]); 
  incr[28] = wr[3]*(1.5*fr[28]-0.8660254037844386*fr[23])+dxvr[3]*(0.4330127018922193*fr[20]-0.25*fr[12]); 
  incr[29] = wr[3]*(0.5*fr[29]-0.8660254037844386*fr[31])+dxvr[3]*(0.1443375672974065*fr[21]-0.25*fr[27]); 
  incr[30] = wr[3]*(1.5*fr[30]-0.8660254037844386*fr[25])+dxvr[3]*(0.4330127018922193*fr[22]-0.25*fr[14]); 
  incr[31] = wr[3]*(1.5*fr[31]-0.8660254037844386*fr[29])+dxvr[3]*(0.4330127018922193*fr[27]-0.25*fr[21]); 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 
  outr[8] += incr[8]*rdxr2; 
  outr[9] += incr[9]*rdxr2; 
  outr[10] += incr[10]*rdxr2; 
  outr[11] += incr[11]*rdxr2; 
  outr[12] += incr[12]*rdxr2; 
  outr[13] += incr[13]*rdxr2; 
  outr[14] += incr[14]*rdxr2; 
  outr[15] += incr[15]*rdxr2; 
  outr[16] += incr[16]*rdxr2; 
  outr[17] += incr[17]*rdxr2; 
  outr[18] += incr[18]*rdxr2; 
  outr[19] += incr[19]*rdxr2; 
  outr[20] += incr[20]*rdxr2; 
  outr[21] += incr[21]*rdxr2; 
  outr[22] += incr[22]*rdxr2; 
  outr[23] += incr[23]*rdxr2; 
  outr[24] += incr[24]*rdxr2; 
  outr[25] += incr[25]*rdxr2; 
  outr[26] += incr[26]*rdxr2; 
  outr[27] += incr[27]*rdxr2; 
  outr[28] += incr[28]*rdxr2; 
  outr[29] += incr[29]*rdxr2; 
  outr[30] += incr[30]*rdxr2; 
  outr[31] += incr[31]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += -1.0*incr[1]*rdxl2; 
  outl[2] += incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += -1.0*incr[7]*rdxl2; 
  outl[8] += incr[8]*rdxl2; 
  outl[9] += -1.0*incr[9]*rdxl2; 
  outl[10] += incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += -1.0*incr[12]*rdxl2; 
  outl[13] += incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += -1.0*incr[15]*rdxl2; 
  outl[16] += incr[16]*rdxl2; 
  outl[17] += incr[17]*rdxl2; 
  outl[18] += -1.0*incr[18]*rdxl2; 
  outl[19] += incr[19]*rdxl2; 
  outl[20] += incr[20]*rdxl2; 
  outl[21] += -1.0*incr[21]*rdxl2; 
  outl[22] += incr[22]*rdxl2; 
  outl[23] += -1.0*incr[23]*rdxl2; 
  outl[24] += incr[24]*rdxl2; 
  outl[25] += -1.0*incr[25]*rdxl2; 
  outl[26] += incr[26]*rdxl2; 
  outl[27] += incr[27]*rdxl2; 
  outl[28] += incr[28]*rdxl2; 
  outl[29] += -1.0*incr[29]*rdxl2; 
  outl[30] += incr[30]*rdxl2; 
  outl[31] += incr[31]*rdxl2; 
  } 
} 
