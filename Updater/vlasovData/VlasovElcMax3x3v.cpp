#include <VlasovModDecl.h> 
void VlasovVolElc3x3vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. E/f: Input electric-field/distribution function. out: Incremented output 
  double dv10 = 2/dxv[3]; 
  const double *E0 = &E[0]; 
  double dv11 = 2/dxv[4]; 
  const double *E1 = &E[4]; 
  double dv12 = 2/dxv[5]; 
  const double *E2 = &E[8]; 
  out[4] += 0.6123724356957944*E0[3]*f[3]*dv10+0.6123724356957944*E0[2]*f[2]*dv10+0.6123724356957944*E0[1]*f[1]*dv10+0.6123724356957944*E0[0]*f[0]*dv10; 
  out[5] += 0.6123724356957944*E1[3]*f[3]*dv11+0.6123724356957944*E1[2]*f[2]*dv11+0.6123724356957944*E1[1]*f[1]*dv11+0.6123724356957944*E1[0]*f[0]*dv11; 
  out[6] += 0.6123724356957944*E2[3]*f[3]*dv12+0.6123724356957944*E2[2]*f[2]*dv12+0.6123724356957944*E2[1]*f[1]*dv12+0.6123724356957944*E2[0]*f[0]*dv12; 
} 
void VlasovVolElc3x3vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. E/f: Input electric-field/distribution function. out: Incremented output 
  double dv10 = 2/dxv[3]; 
  const double *E0 = &E[0]; 
  double dv11 = 2/dxv[4]; 
  const double *E1 = &E[10]; 
  double dv12 = 2/dxv[5]; 
  const double *E2 = &E[20]; 
  out[4] += 0.6123724356957944*E0[9]*f[24]*dv10+0.6123724356957944*E0[8]*f[23]*dv10+0.6123724356957944*E0[7]*f[22]*dv10+0.6123724356957944*E0[6]*f[9]*dv10+0.6123724356957944*E0[5]*f[8]*dv10+0.6123724356957944*E0[4]*f[7]*dv10+0.6123724356957944*E0[3]*f[3]*dv10+0.6123724356957944*E0[2]*f[2]*dv10+0.6123724356957944*E0[1]*f[1]*dv10+0.6123724356957944*E0[0]*f[0]*dv10; 
  out[5] += 0.6123724356957944*E1[9]*f[24]*dv11+0.6123724356957944*E1[8]*f[23]*dv11+0.6123724356957944*E1[7]*f[22]*dv11+0.6123724356957944*E1[6]*f[9]*dv11+0.6123724356957944*E1[5]*f[8]*dv11+0.6123724356957944*E1[4]*f[7]*dv11+0.6123724356957944*E1[3]*f[3]*dv11+0.6123724356957944*E1[2]*f[2]*dv11+0.6123724356957944*E1[1]*f[1]*dv11+0.6123724356957944*E1[0]*f[0]*dv11; 
  out[6] += 0.6123724356957944*E2[9]*f[24]*dv12+0.6123724356957944*E2[8]*f[23]*dv12+0.6123724356957944*E2[7]*f[22]*dv12+0.6123724356957944*E2[6]*f[9]*dv12+0.6123724356957944*E2[5]*f[8]*dv12+0.6123724356957944*E2[4]*f[7]*dv12+0.6123724356957944*E2[3]*f[3]*dv12+0.6123724356957944*E2[2]*f[2]*dv12+0.6123724356957944*E2[1]*f[1]*dv12+0.6123724356957944*E2[0]*f[0]*dv12; 
  out[10] += 0.5477225575051661*E0[1]*f[22]*dv10+0.6123724356957944*E0[3]*f[8]*dv10+0.6123724356957944*E0[2]*f[7]*dv10+0.5477225575051661*f[1]*E0[7]*dv10+0.6123724356957944*f[3]*E0[5]*dv10+0.6123724356957944*f[2]*E0[4]*dv10+0.6123724356957944*E0[0]*f[1]*dv10+0.6123724356957944*f[0]*E0[1]*dv10; 
  out[11] += 0.5477225575051661*E0[2]*f[23]*dv10+0.6123724356957944*E0[3]*f[9]*dv10+0.5477225575051661*f[2]*E0[8]*dv10+0.6123724356957944*E0[1]*f[7]*dv10+0.6123724356957944*f[3]*E0[6]*dv10+0.6123724356957944*f[1]*E0[4]*dv10+0.6123724356957944*E0[0]*f[2]*dv10+0.6123724356957944*f[0]*E0[2]*dv10; 
  out[12] += 0.5477225575051661*E0[3]*f[24]*dv10+0.6123724356957944*E0[2]*f[9]*dv10+0.5477225575051661*f[3]*E0[9]*dv10+0.6123724356957944*E0[1]*f[8]*dv10+0.6123724356957944*f[2]*E0[6]*dv10+0.6123724356957944*f[1]*E0[5]*dv10+0.6123724356957944*E0[0]*f[3]*dv10+0.6123724356957944*f[0]*E0[3]*dv10; 
  out[13] += 0.5477225575051661*E1[1]*f[22]*dv11+0.6123724356957944*E1[3]*f[8]*dv11+0.6123724356957944*E1[2]*f[7]*dv11+0.5477225575051661*f[1]*E1[7]*dv11+0.6123724356957944*f[3]*E1[5]*dv11+0.6123724356957944*f[2]*E1[4]*dv11+0.6123724356957944*E1[0]*f[1]*dv11+0.6123724356957944*f[0]*E1[1]*dv11; 
  out[14] += 0.5477225575051661*E1[2]*f[23]*dv11+0.6123724356957944*E1[3]*f[9]*dv11+0.5477225575051661*f[2]*E1[8]*dv11+0.6123724356957944*E1[1]*f[7]*dv11+0.6123724356957944*f[3]*E1[6]*dv11+0.6123724356957944*f[1]*E1[4]*dv11+0.6123724356957944*E1[0]*f[2]*dv11+0.6123724356957944*f[0]*E1[2]*dv11; 
  out[15] += 0.5477225575051661*E1[3]*f[24]*dv11+0.6123724356957944*E1[2]*f[9]*dv11+0.5477225575051661*f[3]*E1[9]*dv11+0.6123724356957944*E1[1]*f[8]*dv11+0.6123724356957944*f[2]*E1[6]*dv11+0.6123724356957944*f[1]*E1[5]*dv11+0.6123724356957944*E1[0]*f[3]*dv11+0.6123724356957944*f[0]*E1[3]*dv11; 
  out[16] += 0.6123724356957944*E1[3]*f[12]*dv11+0.6123724356957944*E1[2]*f[11]*dv11+0.6123724356957944*E1[1]*f[10]*dv11+0.6123724356957944*E1[0]*f[4]*dv11+0.6123724356957944*E0[3]*f[15]*dv10+0.6123724356957944*E0[2]*f[14]*dv10+0.6123724356957944*E0[1]*f[13]*dv10+0.6123724356957944*E0[0]*f[5]*dv10; 
  out[17] += 0.5477225575051661*E2[1]*f[22]*dv12+0.6123724356957944*E2[3]*f[8]*dv12+0.6123724356957944*E2[2]*f[7]*dv12+0.5477225575051661*f[1]*E2[7]*dv12+0.6123724356957944*f[3]*E2[5]*dv12+0.6123724356957944*f[2]*E2[4]*dv12+0.6123724356957944*E2[0]*f[1]*dv12+0.6123724356957944*f[0]*E2[1]*dv12; 
  out[18] += 0.5477225575051661*E2[2]*f[23]*dv12+0.6123724356957944*E2[3]*f[9]*dv12+0.5477225575051661*f[2]*E2[8]*dv12+0.6123724356957944*E2[1]*f[7]*dv12+0.6123724356957944*f[3]*E2[6]*dv12+0.6123724356957944*f[1]*E2[4]*dv12+0.6123724356957944*E2[0]*f[2]*dv12+0.6123724356957944*f[0]*E2[2]*dv12; 
  out[19] += 0.5477225575051661*E2[3]*f[24]*dv12+0.6123724356957944*E2[2]*f[9]*dv12+0.5477225575051661*f[3]*E2[9]*dv12+0.6123724356957944*E2[1]*f[8]*dv12+0.6123724356957944*f[2]*E2[6]*dv12+0.6123724356957944*f[1]*E2[5]*dv12+0.6123724356957944*E2[0]*f[3]*dv12+0.6123724356957944*f[0]*E2[3]*dv12; 
  out[20] += 0.6123724356957944*E2[3]*f[12]*dv12+0.6123724356957944*E2[2]*f[11]*dv12+0.6123724356957944*E2[1]*f[10]*dv12+0.6123724356957944*E2[0]*f[4]*dv12+0.6123724356957944*E0[3]*f[19]*dv10+0.6123724356957944*E0[2]*f[18]*dv10+0.6123724356957944*E0[1]*f[17]*dv10+0.6123724356957944*E0[0]*f[6]*dv10; 
  out[21] += 0.6123724356957944*E2[3]*f[15]*dv12+0.6123724356957944*E2[2]*f[14]*dv12+0.6123724356957944*E2[1]*f[13]*dv12+0.6123724356957944*E2[0]*f[5]*dv12+0.6123724356957944*E1[3]*f[19]*dv11+0.6123724356957944*E1[2]*f[18]*dv11+0.6123724356957944*E1[1]*f[17]*dv11+0.6123724356957944*E1[0]*f[6]*dv11; 
  out[25] += 1.369306393762915*E0[3]*f[12]*dv10+1.369306393762915*E0[2]*f[11]*dv10+1.369306393762915*E0[1]*f[10]*dv10+1.369306393762915*E0[0]*f[4]*dv10; 
  out[26] += 1.369306393762915*E1[3]*f[15]*dv11+1.369306393762915*E1[2]*f[14]*dv11+1.369306393762915*E1[1]*f[13]*dv11+1.369306393762915*E1[0]*f[5]*dv11; 
  out[27] += 1.369306393762915*E2[3]*f[19]*dv12+1.369306393762915*E2[2]*f[18]*dv12+1.369306393762915*E2[1]*f[17]*dv12+1.369306393762915*E2[0]*f[6]*dv12; 
} 