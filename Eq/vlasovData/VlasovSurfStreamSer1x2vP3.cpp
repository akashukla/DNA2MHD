#include <VlasovModDecl.h> 
__host__ __device__ void VlasovSurfStream1x2vSer_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double incr[32]; 

  if (wr[1]>0) { 
  incr[0] = dxvl[1]*(0.3818813079129866*fl[23]+0.3227486121839514*fl[11]+0.25*fl[4]+0.1443375672974065*fl[2])+wl[1]*(1.322875655532295*fl[17]+1.118033988749895*fl[7]+0.8660254037844386*fl[1]+0.5*fl[0]); 
  incr[1] = dxvl[1]*((-0.6614378277661477*fl[23])-0.5590169943749476*fl[11]-0.4330127018922193*fl[4]-0.25*fl[2])+wl[1]*((-2.29128784747792*fl[17])-1.936491673103709*fl[7]-1.5*fl[1]-0.8660254037844386*fl[0]); 
  incr[2] = wl[1]*(1.322875655532295*fl[23]+1.118033988749895*fl[11]+0.8660254037844386*fl[4]+0.5*fl[2])+dxvl[1]*(0.3818813079129867*fl[17]+0.223606797749979*fl[12]+0.1290994448735806*fl[8]+0.3227486121839515*fl[7]+0.25*fl[1]+0.1443375672974065*fl[0]); 
  incr[3] = dxvl[1]*(0.3818813079129867*fl[29]+0.3227486121839515*fl[20]+0.25*fl[10]+0.1443375672974065*fl[6])+wl[1]*(1.322875655532295*fl[25]+1.118033988749895*fl[13]+0.8660254037844386*fl[5]+0.5*fl[3]); 
  incr[4] = wl[1]*((-2.29128784747792*fl[23])-1.936491673103709*fl[11]-1.5*fl[4]-0.8660254037844386*fl[2])+dxvl[1]*((-0.6614378277661477*fl[17])-0.3872983346207417*fl[12]-0.223606797749979*fl[8]-0.5590169943749475*fl[7]-0.4330127018922193*fl[1]-0.25*fl[0]); 
  incr[5] = dxvl[1]*((-0.6614378277661477*fl[29])-0.5590169943749475*fl[20]-0.4330127018922193*fl[10]-0.25*fl[6])+wl[1]*((-2.29128784747792*fl[25])-1.936491673103709*fl[13]-1.5*fl[5]-0.8660254037844386*fl[3]); 
  incr[6] = wl[1]*(1.322875655532295*fl[29]+1.118033988749895*fl[20]+0.8660254037844386*fl[10]+0.5*fl[6])+dxvl[1]*(0.3818813079129866*fl[25]+0.223606797749979*fl[21]+0.1290994448735805*fl[14]+0.3227486121839514*fl[13]+0.25*fl[5]+0.1443375672974065*fl[3]); 
  incr[7] = dxvl[1]*(0.8539125638299665*fl[23]+0.7216878364870323*fl[11]+0.5590169943749475*fl[4]+0.3227486121839515*fl[2])+wl[1]*(2.958039891549809*fl[17]+2.5*fl[7]+1.936491673103709*fl[1]+1.118033988749895*fl[0]); 
  incr[8] = dxvl[1]*(0.2195775164134199*fl[24]+0.3415650255319865*fl[23]+0.1267731382092775*fl[18]+0.2886751345948129*fl[11]+0.223606797749979*fl[4]+0.1290994448735806*fl[2])+wl[1]*(0.8660254037844387*fl[12]+0.5*fl[8]); 
  incr[9] = dxvl[1]*(0.25*fl[22]+0.1443375672974064*fl[16])+wl[1]*(0.8660254037844387*fl[15]+0.5*fl[9]); 
  incr[10] = wl[1]*((-2.29128784747792*fl[29])-1.936491673103709*fl[20]-1.5*fl[10]-0.8660254037844386*fl[6])+dxvl[1]*((-0.6614378277661477*fl[25])-0.3872983346207416*fl[21]-0.223606797749979*fl[14]-0.5590169943749476*fl[13]-0.4330127018922193*fl[5]-0.25*fl[3]); 
  incr[11] = wl[1]*(2.958039891549808*fl[23]+2.5*fl[11]+1.936491673103709*fl[4]+1.118033988749895*fl[2])+dxvl[1]*(0.8539125638299666*fl[17]+0.5*fl[12]+0.2886751345948129*fl[8]+0.7216878364870323*fl[7]+0.5590169943749476*fl[1]+0.3227486121839514*fl[0]); 
  incr[12] = dxvl[1]*((-0.3803194146278325*fl[24])-0.5916079783099615*fl[23]-0.2195775164134199*fl[18]-0.5*fl[11]-0.3872983346207417*fl[4]-0.223606797749979*fl[2])+wl[1]*((-1.5*fl[12])-0.8660254037844387*fl[8]); 
  incr[13] = dxvl[1]*(0.8539125638299666*fl[29]+0.7216878364870323*fl[20]+0.5590169943749476*fl[10]+0.3227486121839514*fl[6])+wl[1]*(2.958039891549808*fl[25]+2.5*fl[13]+1.936491673103709*fl[5]+1.118033988749895*fl[3]); 
  incr[14] = dxvl[1]*(0.2195775164134199*fl[30]+0.3415650255319866*fl[29]+0.1267731382092775*fl[26]+0.2886751345948129*fl[20]+0.223606797749979*fl[10]+0.1290994448735805*fl[6])+wl[1]*(0.8660254037844387*fl[21]+0.5*fl[14]); 
  incr[15] = dxvl[1]*((-0.4330127018922194*fl[22])-0.25*fl[16])+wl[1]*((-1.5*fl[15])-0.8660254037844387*fl[9]); 
  incr[16] = wl[1]*(0.8660254037844387*fl[22]+0.5*fl[16])+dxvl[1]*(0.25*fl[15]+0.1443375672974064*fl[9]); 
  incr[17] = dxvl[1]*((-1.010362971081845*fl[23])-0.8539125638299666*fl[11]-0.6614378277661477*fl[4]-0.3818813079129867*fl[2])+wl[1]*((-3.5*fl[17])-2.958039891549809*fl[7]-2.29128784747792*fl[1]-1.322875655532295*fl[0]); 
  incr[18] = wl[1]*(0.8660254037844386*fl[24]+0.5*fl[18])+dxvl[1]*(0.2195775164134199*fl[12]+0.1267731382092775*fl[8]); 
  incr[19] = dxvl[1]*(0.25*fl[31]+0.1443375672974064*fl[28])+wl[1]*(0.8660254037844386*fl[27]+0.5*fl[19]); 
  incr[20] = wl[1]*(2.958039891549809*fl[29]+2.5*fl[20]+1.936491673103709*fl[10]+1.118033988749895*fl[6])+dxvl[1]*(0.8539125638299665*fl[25]+0.5*fl[21]+0.2886751345948129*fl[14]+0.7216878364870323*fl[13]+0.5590169943749475*fl[5]+0.3227486121839515*fl[3]); 
  incr[21] = dxvl[1]*((-0.3803194146278324*fl[30])-0.5916079783099616*fl[29]-0.2195775164134199*fl[26]-0.5*fl[20]-0.3872983346207416*fl[10]-0.223606797749979*fl[6])+wl[1]*((-1.5*fl[21])-0.8660254037844387*fl[14]); 
  incr[22] = wl[1]*((-1.5*fl[22])-0.8660254037844387*fl[16])+dxvl[1]*((-0.4330127018922194*fl[15])-0.25*fl[9]); 
  incr[23] = wl[1]*((-3.5*fl[23])-2.958039891549808*fl[11]-2.29128784747792*fl[4]-1.322875655532295*fl[2])+dxvl[1]*((-1.010362971081845*fl[17])-0.5916079783099615*fl[12]-0.3415650255319865*fl[8]-0.8539125638299665*fl[7]-0.6614378277661477*fl[1]-0.3818813079129866*fl[0]); 
  incr[24] = wl[1]*((-1.5*fl[24])-0.8660254037844386*fl[18])+dxvl[1]*((-0.3803194146278325*fl[12])-0.2195775164134199*fl[8]); 
  incr[25] = dxvl[1]*((-1.010362971081845*fl[29])-0.8539125638299665*fl[20]-0.6614378277661477*fl[10]-0.3818813079129866*fl[6])+wl[1]*((-3.5*fl[25])-2.958039891549808*fl[13]-2.29128784747792*fl[5]-1.322875655532295*fl[3]); 
  incr[26] = wl[1]*(0.8660254037844386*fl[30]+0.5*fl[26])+dxvl[1]*(0.2195775164134199*fl[21]+0.1267731382092775*fl[14]); 
  incr[27] = dxvl[1]*((-0.4330127018922193*fl[31])-0.25*fl[28])+wl[1]*((-1.5*fl[27])-0.8660254037844386*fl[19]); 
  incr[28] = wl[1]*(0.8660254037844386*fl[31]+0.5*fl[28])+dxvl[1]*(0.25*fl[27]+0.1443375672974064*fl[19]); 
  incr[29] = wl[1]*((-3.5*fl[29])-2.958039891549809*fl[20]-2.29128784747792*fl[10]-1.322875655532295*fl[6])+dxvl[1]*((-1.010362971081845*fl[25])-0.5916079783099616*fl[21]-0.3415650255319866*fl[14]-0.8539125638299666*fl[13]-0.6614378277661477*fl[5]-0.3818813079129867*fl[3]); 
  incr[30] = wl[1]*((-1.5*fl[30])-0.8660254037844386*fl[26])+dxvl[1]*((-0.3803194146278324*fl[21])-0.2195775164134199*fl[14]); 
  incr[31] = wl[1]*((-1.5*fl[31])-0.8660254037844386*fl[28])+dxvl[1]*((-0.4330127018922193*fl[27])-0.25*fl[19]); 

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
  outl[4] += incr[4]*rdxl2; 
  outl[5] += incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  outl[7] += -1.0*incr[7]*rdxl2; 
  outl[8] += -1.0*incr[8]*rdxl2; 
  outl[9] += -1.0*incr[9]*rdxl2; 
  outl[10] += incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += incr[12]*rdxl2; 
  outl[13] += -1.0*incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += incr[15]*rdxl2; 
  outl[16] += -1.0*incr[16]*rdxl2; 
  outl[17] += incr[17]*rdxl2; 
  outl[18] += -1.0*incr[18]*rdxl2; 
  outl[19] += -1.0*incr[19]*rdxl2; 
  outl[20] += -1.0*incr[20]*rdxl2; 
  outl[21] += incr[21]*rdxl2; 
  outl[22] += incr[22]*rdxl2; 
  outl[23] += incr[23]*rdxl2; 
  outl[24] += incr[24]*rdxl2; 
  outl[25] += incr[25]*rdxl2; 
  outl[26] += -1.0*incr[26]*rdxl2; 
  outl[27] += incr[27]*rdxl2; 
  outl[28] += -1.0*incr[28]*rdxl2; 
  outl[29] += incr[29]*rdxl2; 
  outl[30] += incr[30]*rdxl2; 
  outl[31] += incr[31]*rdxl2; 
  } else { 
  incr[0] = dxvr[1]*((-0.3818813079129866*fr[23])+0.3227486121839514*fr[11]-0.25*fr[4]+0.1443375672974065*fr[2])+wr[1]*((-1.322875655532295*fr[17])+1.118033988749895*fr[7]-0.8660254037844386*fr[1]+0.5*fr[0]); 
  incr[1] = dxvr[1]*(0.6614378277661477*fr[23]-0.5590169943749476*fr[11]+0.4330127018922193*fr[4]-0.25*fr[2])+wr[1]*(2.29128784747792*fr[17]-1.936491673103709*fr[7]+1.5*fr[1]-0.8660254037844386*fr[0]); 
  incr[2] = wr[1]*((-1.322875655532295*fr[23])+1.118033988749895*fr[11]-0.8660254037844386*fr[4]+0.5*fr[2])+dxvr[1]*((-0.3818813079129867*fr[17])-0.223606797749979*fr[12]+0.1290994448735806*fr[8]+0.3227486121839515*fr[7]-0.25*fr[1]+0.1443375672974065*fr[0]); 
  incr[3] = dxvr[1]*((-0.3818813079129867*fr[29])+0.3227486121839515*fr[20]-0.25*fr[10]+0.1443375672974065*fr[6])+wr[1]*((-1.322875655532295*fr[25])+1.118033988749895*fr[13]-0.8660254037844386*fr[5]+0.5*fr[3]); 
  incr[4] = wr[1]*(2.29128784747792*fr[23]-1.936491673103709*fr[11]+1.5*fr[4]-0.8660254037844386*fr[2])+dxvr[1]*(0.6614378277661477*fr[17]+0.3872983346207417*fr[12]-0.223606797749979*fr[8]-0.5590169943749475*fr[7]+0.4330127018922193*fr[1]-0.25*fr[0]); 
  incr[5] = dxvr[1]*(0.6614378277661477*fr[29]-0.5590169943749475*fr[20]+0.4330127018922193*fr[10]-0.25*fr[6])+wr[1]*(2.29128784747792*fr[25]-1.936491673103709*fr[13]+1.5*fr[5]-0.8660254037844386*fr[3]); 
  incr[6] = wr[1]*((-1.322875655532295*fr[29])+1.118033988749895*fr[20]-0.8660254037844386*fr[10]+0.5*fr[6])+dxvr[1]*((-0.3818813079129866*fr[25])-0.223606797749979*fr[21]+0.1290994448735805*fr[14]+0.3227486121839514*fr[13]-0.25*fr[5]+0.1443375672974065*fr[3]); 
  incr[7] = dxvr[1]*((-0.8539125638299665*fr[23])+0.7216878364870323*fr[11]-0.5590169943749475*fr[4]+0.3227486121839515*fr[2])+wr[1]*((-2.958039891549809*fr[17])+2.5*fr[7]-1.936491673103709*fr[1]+1.118033988749895*fr[0]); 
  incr[8] = dxvr[1]*((-0.2195775164134199*fr[24])-0.3415650255319865*fr[23]+0.1267731382092775*fr[18]+0.2886751345948129*fr[11]-0.223606797749979*fr[4]+0.1290994448735806*fr[2])+wr[1]*(0.5*fr[8]-0.8660254037844387*fr[12]); 
  incr[9] = dxvr[1]*(0.1443375672974064*fr[16]-0.25*fr[22])+wr[1]*(0.5*fr[9]-0.8660254037844387*fr[15]); 
  incr[10] = wr[1]*(2.29128784747792*fr[29]-1.936491673103709*fr[20]+1.5*fr[10]-0.8660254037844386*fr[6])+dxvr[1]*(0.6614378277661477*fr[25]+0.3872983346207416*fr[21]-0.223606797749979*fr[14]-0.5590169943749476*fr[13]+0.4330127018922193*fr[5]-0.25*fr[3]); 
  incr[11] = wr[1]*((-2.958039891549808*fr[23])+2.5*fr[11]-1.936491673103709*fr[4]+1.118033988749895*fr[2])+dxvr[1]*((-0.8539125638299666*fr[17])-0.5*fr[12]+0.2886751345948129*fr[8]+0.7216878364870323*fr[7]-0.5590169943749476*fr[1]+0.3227486121839514*fr[0]); 
  incr[12] = dxvr[1]*(0.3803194146278325*fr[24]+0.5916079783099615*fr[23]-0.2195775164134199*fr[18]-0.5*fr[11]+0.3872983346207417*fr[4]-0.223606797749979*fr[2])+wr[1]*(1.5*fr[12]-0.8660254037844387*fr[8]); 
  incr[13] = dxvr[1]*((-0.8539125638299666*fr[29])+0.7216878364870323*fr[20]-0.5590169943749476*fr[10]+0.3227486121839514*fr[6])+wr[1]*((-2.958039891549808*fr[25])+2.5*fr[13]-1.936491673103709*fr[5]+1.118033988749895*fr[3]); 
  incr[14] = dxvr[1]*((-0.2195775164134199*fr[30])-0.3415650255319866*fr[29]+0.1267731382092775*fr[26]+0.2886751345948129*fr[20]-0.223606797749979*fr[10]+0.1290994448735805*fr[6])+wr[1]*(0.5*fr[14]-0.8660254037844387*fr[21]); 
  incr[15] = dxvr[1]*(0.4330127018922194*fr[22]-0.25*fr[16])+wr[1]*(1.5*fr[15]-0.8660254037844387*fr[9]); 
  incr[16] = wr[1]*(0.5*fr[16]-0.8660254037844387*fr[22])+dxvr[1]*(0.1443375672974064*fr[9]-0.25*fr[15]); 
  incr[17] = dxvr[1]*(1.010362971081845*fr[23]-0.8539125638299666*fr[11]+0.6614378277661477*fr[4]-0.3818813079129867*fr[2])+wr[1]*(3.5*fr[17]-2.958039891549809*fr[7]+2.29128784747792*fr[1]-1.322875655532295*fr[0]); 
  incr[18] = wr[1]*(0.5*fr[18]-0.8660254037844386*fr[24])+dxvr[1]*(0.1267731382092775*fr[8]-0.2195775164134199*fr[12]); 
  incr[19] = dxvr[1]*(0.1443375672974064*fr[28]-0.25*fr[31])+wr[1]*(0.5*fr[19]-0.8660254037844386*fr[27]); 
  incr[20] = wr[1]*((-2.958039891549809*fr[29])+2.5*fr[20]-1.936491673103709*fr[10]+1.118033988749895*fr[6])+dxvr[1]*((-0.8539125638299665*fr[25])-0.5*fr[21]+0.2886751345948129*fr[14]+0.7216878364870323*fr[13]-0.5590169943749475*fr[5]+0.3227486121839515*fr[3]); 
  incr[21] = dxvr[1]*(0.3803194146278324*fr[30]+0.5916079783099616*fr[29]-0.2195775164134199*fr[26]-0.5*fr[20]+0.3872983346207416*fr[10]-0.223606797749979*fr[6])+wr[1]*(1.5*fr[21]-0.8660254037844387*fr[14]); 
  incr[22] = wr[1]*(1.5*fr[22]-0.8660254037844387*fr[16])+dxvr[1]*(0.4330127018922194*fr[15]-0.25*fr[9]); 
  incr[23] = wr[1]*(3.5*fr[23]-2.958039891549808*fr[11]+2.29128784747792*fr[4]-1.322875655532295*fr[2])+dxvr[1]*(1.010362971081845*fr[17]+0.5916079783099615*fr[12]-0.3415650255319865*fr[8]-0.8539125638299665*fr[7]+0.6614378277661477*fr[1]-0.3818813079129866*fr[0]); 
  incr[24] = wr[1]*(1.5*fr[24]-0.8660254037844386*fr[18])+dxvr[1]*(0.3803194146278325*fr[12]-0.2195775164134199*fr[8]); 
  incr[25] = dxvr[1]*(1.010362971081845*fr[29]-0.8539125638299665*fr[20]+0.6614378277661477*fr[10]-0.3818813079129866*fr[6])+wr[1]*(3.5*fr[25]-2.958039891549808*fr[13]+2.29128784747792*fr[5]-1.322875655532295*fr[3]); 
  incr[26] = wr[1]*(0.5*fr[26]-0.8660254037844386*fr[30])+dxvr[1]*(0.1267731382092775*fr[14]-0.2195775164134199*fr[21]); 
  incr[27] = dxvr[1]*(0.4330127018922193*fr[31]-0.25*fr[28])+wr[1]*(1.5*fr[27]-0.8660254037844386*fr[19]); 
  incr[28] = wr[1]*(0.5*fr[28]-0.8660254037844386*fr[31])+dxvr[1]*(0.1443375672974064*fr[19]-0.25*fr[27]); 
  incr[29] = wr[1]*(3.5*fr[29]-2.958039891549809*fr[20]+2.29128784747792*fr[10]-1.322875655532295*fr[6])+dxvr[1]*(1.010362971081845*fr[25]+0.5916079783099616*fr[21]-0.3415650255319866*fr[14]-0.8539125638299666*fr[13]+0.6614378277661477*fr[5]-0.3818813079129867*fr[3]); 
  incr[30] = wr[1]*(1.5*fr[30]-0.8660254037844386*fr[26])+dxvr[1]*(0.3803194146278324*fr[21]-0.2195775164134199*fr[14]); 
  incr[31] = wr[1]*(1.5*fr[31]-0.8660254037844386*fr[28])+dxvr[1]*(0.4330127018922193*fr[27]-0.25*fr[19]); 

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
  outl[4] += incr[4]*rdxl2; 
  outl[5] += incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  outl[7] += -1.0*incr[7]*rdxl2; 
  outl[8] += -1.0*incr[8]*rdxl2; 
  outl[9] += -1.0*incr[9]*rdxl2; 
  outl[10] += incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += incr[12]*rdxl2; 
  outl[13] += -1.0*incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += incr[15]*rdxl2; 
  outl[16] += -1.0*incr[16]*rdxl2; 
  outl[17] += incr[17]*rdxl2; 
  outl[18] += -1.0*incr[18]*rdxl2; 
  outl[19] += -1.0*incr[19]*rdxl2; 
  outl[20] += -1.0*incr[20]*rdxl2; 
  outl[21] += incr[21]*rdxl2; 
  outl[22] += incr[22]*rdxl2; 
  outl[23] += incr[23]*rdxl2; 
  outl[24] += incr[24]*rdxl2; 
  outl[25] += incr[25]*rdxl2; 
  outl[26] += -1.0*incr[26]*rdxl2; 
  outl[27] += incr[27]*rdxl2; 
  outl[28] += -1.0*incr[28]*rdxl2; 
  outl[29] += incr[29]*rdxl2; 
  outl[30] += incr[30]*rdxl2; 
  outl[31] += incr[31]*rdxl2; 
  } 
} 