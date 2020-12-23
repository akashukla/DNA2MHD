#include <VlasovModDecl.h> 
__host__ __device__ void VlasovSurfStream1x2vMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double incr[10]; 

  if (wr[1]>0) { 
  incr[0] = wl[1]*(1.118033988749895*fl[7]+0.8660254037844386*fl[1]+0.5*fl[0])+dxvl[1]*(0.25*fl[4]+0.1443375672974065*fl[2]); 
  incr[1] = wl[1]*((-1.936491673103709*fl[7])-1.5*fl[1]-0.8660254037844386*fl[0])+dxvl[1]*((-0.4330127018922193*fl[4])-0.25*fl[2]); 
  incr[2] = dxvl[1]*(0.1290994448735806*fl[8]+0.3227486121839515*fl[7]+0.25*fl[1]+0.1443375672974065*fl[0])+wl[1]*(0.8660254037844386*fl[4]+0.5*fl[2]); 
  incr[3] = 0.1443375672974065*dxvl[1]*fl[6]+wl[1]*(0.8660254037844386*fl[5]+0.5*fl[3]); 
  incr[4] = dxvl[1]*((-0.223606797749979*fl[8])-0.5590169943749475*fl[7]-0.4330127018922193*fl[1]-0.25*fl[0])+wl[1]*((-1.5*fl[4])-0.8660254037844386*fl[2]); 
  incr[5] = wl[1]*((-1.5*fl[5])-0.8660254037844386*fl[3])-0.25*dxvl[1]*fl[6]; 
  incr[6] = 0.5*wl[1]*fl[6]+dxvl[1]*(0.25*fl[5]+0.1443375672974065*fl[3]); 
  incr[7] = wl[1]*(2.5*fl[7]+1.936491673103709*fl[1]+1.118033988749895*fl[0])+dxvl[1]*(0.5590169943749475*fl[4]+0.3227486121839515*fl[2]); 
  incr[8] = 0.5*wl[1]*fl[8]+dxvl[1]*(0.223606797749979*fl[4]+0.1290994448735806*fl[2]); 
  incr[9] = 0.5*wl[1]*fl[9]; 

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
  } else { 
  incr[0] = wr[1]*(1.118033988749895*fr[7]-0.8660254037844386*fr[1]+0.5*fr[0])+dxvr[1]*(0.1443375672974065*fr[2]-0.25*fr[4]); 
  incr[1] = wr[1]*((-1.936491673103709*fr[7])+1.5*fr[1]-0.8660254037844386*fr[0])+dxvr[1]*(0.4330127018922193*fr[4]-0.25*fr[2]); 
  incr[2] = dxvr[1]*(0.1290994448735806*fr[8]+0.3227486121839515*fr[7]-0.25*fr[1]+0.1443375672974065*fr[0])+wr[1]*(0.5*fr[2]-0.8660254037844386*fr[4]); 
  incr[3] = 0.1443375672974065*dxvr[1]*fr[6]+wr[1]*(0.5*fr[3]-0.8660254037844386*fr[5]); 
  incr[4] = dxvr[1]*((-0.223606797749979*fr[8])-0.5590169943749475*fr[7]+0.4330127018922193*fr[1]-0.25*fr[0])+wr[1]*(1.5*fr[4]-0.8660254037844386*fr[2]); 
  incr[5] = wr[1]*(1.5*fr[5]-0.8660254037844386*fr[3])-0.25*dxvr[1]*fr[6]; 
  incr[6] = 0.5*wr[1]*fr[6]+dxvr[1]*(0.1443375672974065*fr[3]-0.25*fr[5]); 
  incr[7] = wr[1]*(2.5*fr[7]-1.936491673103709*fr[1]+1.118033988749895*fr[0])+dxvr[1]*(0.3227486121839515*fr[2]-0.5590169943749475*fr[4]); 
  incr[8] = 0.5*wr[1]*fr[8]+dxvr[1]*(0.1290994448735806*fr[2]-0.223606797749979*fr[4]); 
  incr[9] = 0.5*wr[1]*fr[9]; 

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
  } 
} 
