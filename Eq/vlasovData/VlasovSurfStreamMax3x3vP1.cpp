#include <VlasovModDecl.h> 
__host__ __device__ void VlasovSurfStream3x3vMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double incr[7]; 

  if (wr[3]>0) { 
  incr[0] = 0.1443375672974065*dxvl[3]*fl[4]+(0.8660254037844386*fl[1]+0.5*fl[0])*wl[3]; 
  incr[1] = ((-1.5*fl[1])-0.8660254037844386*fl[0])*wl[3]-0.25*dxvl[3]*fl[4]; 
  incr[2] = 0.5*fl[2]*wl[3]; 
  incr[3] = 0.5*fl[3]*wl[3]; 
  incr[4] = 0.5*wl[3]*fl[4]+(0.25*fl[1]+0.1443375672974065*fl[0])*dxvl[3]; 
  incr[5] = 0.5*wl[3]*fl[5]; 
  incr[6] = 0.5*wl[3]*fl[6]; 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  } else { 
  incr[0] = 0.1443375672974065*dxvr[3]*fr[4]+(0.5*fr[0]-0.8660254037844386*fr[1])*wr[3]; 
  incr[1] = (1.5*fr[1]-0.8660254037844386*fr[0])*wr[3]-0.25*dxvr[3]*fr[4]; 
  incr[2] = 0.5*fr[2]*wr[3]; 
  incr[3] = 0.5*fr[3]*wr[3]; 
  incr[4] = 0.5*wr[3]*fr[4]+(0.1443375672974065*fr[0]-0.25*fr[1])*dxvr[3]; 
  incr[5] = 0.5*wr[3]*fr[5]; 
  incr[6] = 0.5*wr[3]*fr[6]; 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  } 
} 
__host__ __device__ void VlasovSurfStream3x3vMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[1]; 
  double rdxr2 = 2.0/dxvr[1]; 

  double incr[7]; 

  if (wr[4]>0) { 
  incr[0] = 0.1443375672974065*dxvl[4]*fl[5]+(0.8660254037844386*fl[2]+0.5*fl[0])*wl[4]; 
  incr[1] = 0.5*fl[1]*wl[4]; 
  incr[2] = ((-1.5*fl[2])-0.8660254037844386*fl[0])*wl[4]-0.25*dxvl[4]*fl[5]; 
  incr[3] = 0.5*fl[3]*wl[4]; 
  incr[4] = 0.5*fl[4]*wl[4]; 
  incr[5] = 0.5*wl[4]*fl[5]+(0.25*fl[2]+0.1443375672974065*fl[0])*dxvl[4]; 
  incr[6] = 0.5*wl[4]*fl[6]; 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += -1.0*incr[1]*rdxl2; 
  outl[2] += incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  } else { 
  incr[0] = 0.1443375672974065*dxvr[4]*fr[5]+(0.5*fr[0]-0.8660254037844386*fr[2])*wr[4]; 
  incr[1] = 0.5*fr[1]*wr[4]; 
  incr[2] = (1.5*fr[2]-0.8660254037844386*fr[0])*wr[4]-0.25*dxvr[4]*fr[5]; 
  incr[3] = 0.5*fr[3]*wr[4]; 
  incr[4] = 0.5*fr[4]*wr[4]; 
  incr[5] = 0.5*wr[4]*fr[5]+(0.1443375672974065*fr[0]-0.25*fr[2])*dxvr[4]; 
  incr[6] = 0.5*wr[4]*fr[6]; 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += -1.0*incr[1]*rdxl2; 
  outl[2] += incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  } 
} 
__host__ __device__ void VlasovSurfStream3x3vMax_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[2]; 
  double rdxr2 = 2.0/dxvr[2]; 

  double incr[7]; 

  if (wr[5]>0) { 
  incr[0] = 0.1443375672974065*dxvl[5]*fl[6]+(0.8660254037844386*fl[3]+0.5*fl[0])*wl[5]; 
  incr[1] = 0.5*fl[1]*wl[5]; 
  incr[2] = 0.5*fl[2]*wl[5]; 
  incr[3] = ((-1.5*fl[3])-0.8660254037844386*fl[0])*wl[5]-0.25*dxvl[5]*fl[6]; 
  incr[4] = 0.5*fl[4]*wl[5]; 
  incr[5] = 0.5*fl[5]*wl[5]; 
  incr[6] = 0.5*wl[5]*fl[6]+(0.25*fl[3]+0.1443375672974065*fl[0])*dxvl[5]; 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += -1.0*incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  } else { 
  incr[0] = 0.1443375672974065*dxvr[5]*fr[6]+(0.5*fr[0]-0.8660254037844386*fr[3])*wr[5]; 
  incr[1] = 0.5*fr[1]*wr[5]; 
  incr[2] = 0.5*fr[2]*wr[5]; 
  incr[3] = (1.5*fr[3]-0.8660254037844386*fr[0])*wr[5]-0.25*dxvr[5]*fr[6]; 
  incr[4] = 0.5*fr[4]*wr[5]; 
  incr[5] = 0.5*fr[5]*wr[5]; 
  incr[6] = 0.5*wr[5]*fr[6]+(0.1443375672974065*fr[0]-0.25*fr[3])*dxvr[5]; 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += -1.0*incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += -1.0*incr[6]*rdxl2; 
  } 
} 
