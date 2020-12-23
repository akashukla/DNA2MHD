#include <VlasovModDecl.h> 
__host__ __device__ void VlasovSurfStream1x3vMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double incr[35]; 

  if (wr[1]>0) { 
  incr[0] = wl[1]*(1.322875655532295*fl[31]+1.118033988749895*fl[11]+0.8660254037844386*fl[1]+0.5*fl[0])+dxvl[1]*(0.3227486121839514*fl[19]+0.25*fl[5]+0.1443375672974065*fl[2]); 
  incr[1] = wl[1]*((-2.29128784747792*fl[31])-1.936491673103709*fl[11]-1.5*fl[1]-0.8660254037844386*fl[0])+dxvl[1]*((-0.5590169943749476*fl[19])-0.4330127018922193*fl[5]-0.25*fl[2]); 
  incr[2] = dxvl[1]*(0.3818813079129867*fl[31]+0.223606797749979*fl[20]+0.1290994448735806*fl[12]+0.3227486121839515*fl[11]+0.25*fl[1]+0.1443375672974065*fl[0])+wl[1]*(1.118033988749895*fl[19]+0.8660254037844386*fl[5]+0.5*fl[2]); 
  incr[3] = wl[1]*(1.118033988749895*fl[21]+0.8660254037844386*fl[6]+0.5*fl[3])+dxvl[1]*(0.25*fl[15]+0.1443375672974065*fl[7]); 
  incr[4] = wl[1]*(1.118033988749895*fl[25]+0.8660254037844386*fl[8]+0.5*fl[4])+dxvl[1]*(0.25*fl[16]+0.1443375672974065*fl[9]); 
  incr[5] = dxvl[1]*((-0.6614378277661477*fl[31])-0.3872983346207417*fl[20]-0.223606797749979*fl[12]-0.5590169943749475*fl[11]-0.4330127018922193*fl[1]-0.25*fl[0])+wl[1]*((-1.936491673103709*fl[19])-1.5*fl[5]-0.8660254037844386*fl[2]); 
  incr[6] = wl[1]*((-1.936491673103709*fl[21])-1.5*fl[6]-0.8660254037844386*fl[3])+dxvl[1]*((-0.4330127018922193*fl[15])-0.25*fl[7]); 
  incr[7] = dxvl[1]*(0.1290994448735805*fl[22]+0.3227486121839514*fl[21]+0.25*fl[6]+0.1443375672974065*fl[3])+wl[1]*(0.8660254037844386*fl[15]+0.5*fl[7]); 
  incr[8] = wl[1]*((-1.936491673103709*fl[25])-1.5*fl[8]-0.8660254037844386*fl[4])+dxvl[1]*((-0.4330127018922193*fl[16])-0.25*fl[9]); 
  incr[9] = dxvl[1]*(0.1290994448735805*fl[26]+0.3227486121839514*fl[25]+0.25*fl[8]+0.1443375672974065*fl[4])+wl[1]*(0.8660254037844386*fl[16]+0.5*fl[9]); 
  incr[10] = 0.1443375672974065*dxvl[1]*fl[18]+wl[1]*(0.8660254037844386*fl[17]+0.5*fl[10]); 
  incr[11] = wl[1]*(2.958039891549809*fl[31]+2.5*fl[11]+1.936491673103709*fl[1]+1.118033988749895*fl[0])+dxvl[1]*(0.7216878364870323*fl[19]+0.5590169943749475*fl[5]+0.3227486121839515*fl[2]); 
  incr[12] = dxvl[1]*(0.1267731382092775*fl[32]+0.2886751345948129*fl[19]+0.223606797749979*fl[5]+0.1290994448735806*fl[2])+wl[1]*(0.8660254037844387*fl[20]+0.5*fl[12]); 
  incr[13] = 0.1443375672974064*dxvl[1]*fl[24]+wl[1]*(0.8660254037844387*fl[23]+0.5*fl[13]); 
  incr[14] = 0.1443375672974064*dxvl[1]*fl[29]+wl[1]*(0.8660254037844387*fl[28]+0.5*fl[14]); 
  incr[15] = dxvl[1]*((-0.223606797749979*fl[22])-0.5590169943749476*fl[21]-0.4330127018922193*fl[6]-0.25*fl[3])+wl[1]*((-1.5*fl[15])-0.8660254037844386*fl[7]); 
  incr[16] = dxvl[1]*((-0.223606797749979*fl[26])-0.5590169943749476*fl[25]-0.4330127018922193*fl[8]-0.25*fl[4])+wl[1]*((-1.5*fl[16])-0.8660254037844386*fl[9]); 
  incr[17] = wl[1]*((-1.5*fl[17])-0.8660254037844386*fl[10])-0.25*dxvl[1]*fl[18]; 
  incr[18] = 0.5*wl[1]*fl[18]+dxvl[1]*(0.25*fl[17]+0.1443375672974065*fl[10]); 
  incr[19] = dxvl[1]*(0.8539125638299666*fl[31]+0.5*fl[20]+0.2886751345948129*fl[12]+0.7216878364870323*fl[11]+0.5590169943749476*fl[1]+0.3227486121839514*fl[0])+wl[1]*(2.5*fl[19]+1.936491673103709*fl[5]+1.118033988749895*fl[2]); 
  incr[20] = dxvl[1]*((-0.2195775164134199*fl[32])-0.5*fl[19]-0.3872983346207417*fl[5]-0.223606797749979*fl[2])+wl[1]*((-1.5*fl[20])-0.8660254037844387*fl[12]); 
  incr[21] = wl[1]*(2.5*fl[21]+1.936491673103709*fl[6]+1.118033988749895*fl[3])+dxvl[1]*(0.5590169943749476*fl[15]+0.3227486121839514*fl[7]); 
  incr[22] = 0.5*wl[1]*fl[22]+dxvl[1]*(0.223606797749979*fl[15]+0.1290994448735805*fl[7]); 
  incr[23] = wl[1]*((-1.5*fl[23])-0.8660254037844387*fl[13])-0.25*dxvl[1]*fl[24]; 
  incr[24] = 0.5*wl[1]*fl[24]+dxvl[1]*(0.25*fl[23]+0.1443375672974064*fl[13]); 
  incr[25] = wl[1]*(2.5*fl[25]+1.936491673103709*fl[8]+1.118033988749895*fl[4])+dxvl[1]*(0.5590169943749476*fl[16]+0.3227486121839514*fl[9]); 
  incr[26] = 0.5*wl[1]*fl[26]+dxvl[1]*(0.223606797749979*fl[16]+0.1290994448735805*fl[9]); 
  incr[27] = 0.5*wl[1]*fl[27]; 
  incr[28] = wl[1]*((-1.5*fl[28])-0.8660254037844387*fl[14])-0.25*dxvl[1]*fl[29]; 
  incr[29] = 0.5*wl[1]*fl[29]+dxvl[1]*(0.25*fl[28]+0.1443375672974064*fl[14]); 
  incr[30] = 0.5*wl[1]*fl[30]; 
  incr[31] = wl[1]*((-3.5*fl[31])-2.958039891549809*fl[11]-2.29128784747792*fl[1]-1.322875655532295*fl[0])+dxvl[1]*((-0.8539125638299666*fl[19])-0.6614378277661477*fl[5]-0.3818813079129867*fl[2]); 
  incr[32] = 0.5*wl[1]*fl[32]+dxvl[1]*(0.2195775164134199*fl[20]+0.1267731382092775*fl[12]); 
  incr[33] = 0.5*wl[1]*fl[33]; 
  incr[34] = 0.5*wl[1]*fl[34]; 

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
  outr[32] += incr[32]*rdxr2; 
  outr[33] += incr[33]*rdxr2; 
  outr[34] += incr[34]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += -1.0*incr[7]*rdxl2; 
  outl[8] += incr[8]*rdxl2; 
  outl[9] += -1.0*incr[9]*rdxl2; 
  outl[10] += -1.0*incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += -1.0*incr[12]*rdxl2; 
  outl[13] += -1.0*incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += incr[15]*rdxl2; 
  outl[16] += incr[16]*rdxl2; 
  outl[17] += incr[17]*rdxl2; 
  outl[18] += -1.0*incr[18]*rdxl2; 
  outl[19] += -1.0*incr[19]*rdxl2; 
  outl[20] += incr[20]*rdxl2; 
  outl[21] += -1.0*incr[21]*rdxl2; 
  outl[22] += -1.0*incr[22]*rdxl2; 
  outl[23] += incr[23]*rdxl2; 
  outl[24] += -1.0*incr[24]*rdxl2; 
  outl[25] += -1.0*incr[25]*rdxl2; 
  outl[26] += -1.0*incr[26]*rdxl2; 
  outl[27] += -1.0*incr[27]*rdxl2; 
  outl[28] += incr[28]*rdxl2; 
  outl[29] += -1.0*incr[29]*rdxl2; 
  outl[30] += -1.0*incr[30]*rdxl2; 
  outl[31] += incr[31]*rdxl2; 
  outl[32] += -1.0*incr[32]*rdxl2; 
  outl[33] += -1.0*incr[33]*rdxl2; 
  outl[34] += -1.0*incr[34]*rdxl2; 
  } else { 
  incr[0] = wr[1]*((-1.322875655532295*fr[31])+1.118033988749895*fr[11]-0.8660254037844386*fr[1]+0.5*fr[0])+dxvr[1]*(0.3227486121839514*fr[19]-0.25*fr[5]+0.1443375672974065*fr[2]); 
  incr[1] = wr[1]*(2.29128784747792*fr[31]-1.936491673103709*fr[11]+1.5*fr[1]-0.8660254037844386*fr[0])+dxvr[1]*((-0.5590169943749476*fr[19])+0.4330127018922193*fr[5]-0.25*fr[2]); 
  incr[2] = dxvr[1]*((-0.3818813079129867*fr[31])-0.223606797749979*fr[20]+0.1290994448735806*fr[12]+0.3227486121839515*fr[11]-0.25*fr[1]+0.1443375672974065*fr[0])+wr[1]*(1.118033988749895*fr[19]-0.8660254037844386*fr[5]+0.5*fr[2]); 
  incr[3] = wr[1]*(1.118033988749895*fr[21]-0.8660254037844386*fr[6]+0.5*fr[3])+dxvr[1]*(0.1443375672974065*fr[7]-0.25*fr[15]); 
  incr[4] = wr[1]*(1.118033988749895*fr[25]-0.8660254037844386*fr[8]+0.5*fr[4])+dxvr[1]*(0.1443375672974065*fr[9]-0.25*fr[16]); 
  incr[5] = dxvr[1]*(0.6614378277661477*fr[31]+0.3872983346207417*fr[20]-0.223606797749979*fr[12]-0.5590169943749475*fr[11]+0.4330127018922193*fr[1]-0.25*fr[0])+wr[1]*((-1.936491673103709*fr[19])+1.5*fr[5]-0.8660254037844386*fr[2]); 
  incr[6] = wr[1]*((-1.936491673103709*fr[21])+1.5*fr[6]-0.8660254037844386*fr[3])+dxvr[1]*(0.4330127018922193*fr[15]-0.25*fr[7]); 
  incr[7] = dxvr[1]*(0.1290994448735805*fr[22]+0.3227486121839514*fr[21]-0.25*fr[6]+0.1443375672974065*fr[3])+wr[1]*(0.5*fr[7]-0.8660254037844386*fr[15]); 
  incr[8] = wr[1]*((-1.936491673103709*fr[25])+1.5*fr[8]-0.8660254037844386*fr[4])+dxvr[1]*(0.4330127018922193*fr[16]-0.25*fr[9]); 
  incr[9] = dxvr[1]*(0.1290994448735805*fr[26]+0.3227486121839514*fr[25]-0.25*fr[8]+0.1443375672974065*fr[4])+wr[1]*(0.5*fr[9]-0.8660254037844386*fr[16]); 
  incr[10] = 0.1443375672974065*dxvr[1]*fr[18]+wr[1]*(0.5*fr[10]-0.8660254037844386*fr[17]); 
  incr[11] = wr[1]*((-2.958039891549809*fr[31])+2.5*fr[11]-1.936491673103709*fr[1]+1.118033988749895*fr[0])+dxvr[1]*(0.7216878364870323*fr[19]-0.5590169943749475*fr[5]+0.3227486121839515*fr[2]); 
  incr[12] = dxvr[1]*(0.1267731382092775*fr[32]+0.2886751345948129*fr[19]-0.223606797749979*fr[5]+0.1290994448735806*fr[2])+wr[1]*(0.5*fr[12]-0.8660254037844387*fr[20]); 
  incr[13] = 0.1443375672974064*dxvr[1]*fr[24]+wr[1]*(0.5*fr[13]-0.8660254037844387*fr[23]); 
  incr[14] = 0.1443375672974064*dxvr[1]*fr[29]+wr[1]*(0.5*fr[14]-0.8660254037844387*fr[28]); 
  incr[15] = dxvr[1]*((-0.223606797749979*fr[22])-0.5590169943749476*fr[21]+0.4330127018922193*fr[6]-0.25*fr[3])+wr[1]*(1.5*fr[15]-0.8660254037844386*fr[7]); 
  incr[16] = dxvr[1]*((-0.223606797749979*fr[26])-0.5590169943749476*fr[25]+0.4330127018922193*fr[8]-0.25*fr[4])+wr[1]*(1.5*fr[16]-0.8660254037844386*fr[9]); 
  incr[17] = wr[1]*(1.5*fr[17]-0.8660254037844386*fr[10])-0.25*dxvr[1]*fr[18]; 
  incr[18] = 0.5*wr[1]*fr[18]+dxvr[1]*(0.1443375672974065*fr[10]-0.25*fr[17]); 
  incr[19] = dxvr[1]*((-0.8539125638299666*fr[31])-0.5*fr[20]+0.2886751345948129*fr[12]+0.7216878364870323*fr[11]-0.5590169943749476*fr[1]+0.3227486121839514*fr[0])+wr[1]*(2.5*fr[19]-1.936491673103709*fr[5]+1.118033988749895*fr[2]); 
  incr[20] = dxvr[1]*((-0.2195775164134199*fr[32])-0.5*fr[19]+0.3872983346207417*fr[5]-0.223606797749979*fr[2])+wr[1]*(1.5*fr[20]-0.8660254037844387*fr[12]); 
  incr[21] = wr[1]*(2.5*fr[21]-1.936491673103709*fr[6]+1.118033988749895*fr[3])+dxvr[1]*(0.3227486121839514*fr[7]-0.5590169943749476*fr[15]); 
  incr[22] = 0.5*wr[1]*fr[22]+dxvr[1]*(0.1290994448735805*fr[7]-0.223606797749979*fr[15]); 
  incr[23] = wr[1]*(1.5*fr[23]-0.8660254037844387*fr[13])-0.25*dxvr[1]*fr[24]; 
  incr[24] = 0.5*wr[1]*fr[24]+dxvr[1]*(0.1443375672974064*fr[13]-0.25*fr[23]); 
  incr[25] = wr[1]*(2.5*fr[25]-1.936491673103709*fr[8]+1.118033988749895*fr[4])+dxvr[1]*(0.3227486121839514*fr[9]-0.5590169943749476*fr[16]); 
  incr[26] = 0.5*wr[1]*fr[26]+dxvr[1]*(0.1290994448735805*fr[9]-0.223606797749979*fr[16]); 
  incr[27] = 0.5*wr[1]*fr[27]; 
  incr[28] = wr[1]*(1.5*fr[28]-0.8660254037844387*fr[14])-0.25*dxvr[1]*fr[29]; 
  incr[29] = 0.5*wr[1]*fr[29]+dxvr[1]*(0.1443375672974064*fr[14]-0.25*fr[28]); 
  incr[30] = 0.5*wr[1]*fr[30]; 
  incr[31] = wr[1]*(3.5*fr[31]-2.958039891549809*fr[11]+2.29128784747792*fr[1]-1.322875655532295*fr[0])+dxvr[1]*((-0.8539125638299666*fr[19])+0.6614378277661477*fr[5]-0.3818813079129867*fr[2]); 
  incr[32] = 0.5*wr[1]*fr[32]+dxvr[1]*(0.1267731382092775*fr[12]-0.2195775164134199*fr[20]); 
  incr[33] = 0.5*wr[1]*fr[33]; 
  incr[34] = 0.5*wr[1]*fr[34]; 

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
  outr[32] += incr[32]*rdxr2; 
  outr[33] += incr[33]*rdxr2; 
  outr[34] += incr[34]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += -1.0*incr[7]*rdxl2; 
  outl[8] += incr[8]*rdxl2; 
  outl[9] += -1.0*incr[9]*rdxl2; 
  outl[10] += -1.0*incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += -1.0*incr[12]*rdxl2; 
  outl[13] += -1.0*incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += incr[15]*rdxl2; 
  outl[16] += incr[16]*rdxl2; 
  outl[17] += incr[17]*rdxl2; 
  outl[18] += -1.0*incr[18]*rdxl2; 
  outl[19] += -1.0*incr[19]*rdxl2; 
  outl[20] += incr[20]*rdxl2; 
  outl[21] += -1.0*incr[21]*rdxl2; 
  outl[22] += -1.0*incr[22]*rdxl2; 
  outl[23] += incr[23]*rdxl2; 
  outl[24] += -1.0*incr[24]*rdxl2; 
  outl[25] += -1.0*incr[25]*rdxl2; 
  outl[26] += -1.0*incr[26]*rdxl2; 
  outl[27] += -1.0*incr[27]*rdxl2; 
  outl[28] += incr[28]*rdxl2; 
  outl[29] += -1.0*incr[29]*rdxl2; 
  outl[30] += -1.0*incr[30]*rdxl2; 
  outl[31] += incr[31]*rdxl2; 
  outl[32] += -1.0*incr[32]*rdxl2; 
  outl[33] += -1.0*incr[33]*rdxl2; 
  outl[34] += -1.0*incr[34]*rdxl2; 
  } 
} 
