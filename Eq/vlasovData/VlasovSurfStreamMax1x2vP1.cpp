#include <VlasovModDecl.h> 
__host__ __device__ void VlasovSurfStream1x2vMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double incr[4]; 

  if (wr[1]>0) { 
  incr[0] = 0.1443375672974065*dxvl[1]*fl[2]+(0.8660254037844386*fl[1]+0.5*fl[0])*wl[1]; 
  incr[1] = ((-1.5*fl[1])-0.8660254037844386*fl[0])*wl[1]-0.25*dxvl[1]*fl[2]; 
  incr[2] = 0.5*wl[1]*fl[2]+dxvl[1]*(0.25*fl[1]+0.1443375672974065*fl[0]); 
  incr[3] = 0.5*wl[1]*fl[3]; 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  } else { 
  incr[0] = 0.1443375672974065*dxvr[1]*fr[2]+(0.5*fr[0]-0.8660254037844386*fr[1])*wr[1]; 
  incr[1] = (1.5*fr[1]-0.8660254037844386*fr[0])*wr[1]-0.25*dxvr[1]*fr[2]; 
  incr[2] = 0.5*wr[1]*fr[2]+dxvr[1]*(0.1443375672974065*fr[0]-0.25*fr[1]); 
  incr[3] = 0.5*wr[1]*fr[3]; 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  } 
} 
