#include <PassiveAdvectionModDecl.h> 
double PassiveAdvectionVol4xSerP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing.
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  double dfac2 = 2.0/dxv[1]; 
  double w2 = w[1]; 
  double dfac3 = 2.0/dxv[2]; 
  double w3 = w[2]; 
  double dfac4 = 2.0/dxv[3]; 
  double w4 = w[3]; 
  const double *v1 = &f[16]; 
  const double *v2 = &f[32]; 
  const double *v3 = &f[48]; 
  const double *v4 = &f[64]; 
  double cflRate = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*v1[1]-1.0*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*v1[1]+v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.4330127018922193*v1[15]-0.25*v1[14]-0.4330127018922193*(v1[13]+v1[12]+v1[11])+0.25*(v1[10]+v1[9])+0.4330127018922193*v1[8]+0.25*v1[7]+0.4330127018922193*(v1[6]+v1[5])-0.25*(v1[4]+v1[3]+v1[2])-0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*v1[15])+0.25*v1[14]-0.4330127018922193*v1[13]+0.4330127018922193*(v1[12]+v1[11])+0.25*v1[10]-0.25*v1[9]+0.4330127018922193*v1[8]-0.25*v1[7]+0.4330127018922193*v1[6]-0.4330127018922193*v1[5]-0.25*(v1[4]+v1[3])+0.25*v1[2]-0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*v1[15])+0.25*v1[14]+0.4330127018922193*v1[13]-0.4330127018922193*v1[12]+0.4330127018922193*v1[11]-0.25*v1[10]+0.25*v1[9]+0.4330127018922193*v1[8]-0.25*v1[7]-0.4330127018922193*v1[6]+0.4330127018922193*v1[5]-0.25*v1[4]+0.25*v1[3]-0.25*v1[2]-0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*v1[15]-0.25*v1[14]+0.4330127018922193*(v1[13]+v1[12])-0.4330127018922193*v1[11]-0.25*(v1[10]+v1[9])+0.4330127018922193*v1[8]+0.25*v1[7]-0.4330127018922193*(v1[6]+v1[5])-0.25*v1[4]+0.25*(v1[3]+v1[2])-0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*v1[15])+0.25*v1[14]+0.4330127018922193*(v1[13]+v1[12])-0.4330127018922193*v1[11]-0.25*(v1[10]+v1[9])-0.4330127018922193*v1[8]+0.25*v1[7]+0.4330127018922193*(v1[6]+v1[5])+0.25*v1[4]-0.25*(v1[3]+v1[2])-0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*v1[15]-0.25*v1[14]+0.4330127018922193*v1[13]-0.4330127018922193*v1[12]+0.4330127018922193*v1[11]-0.25*v1[10]+0.25*v1[9]-0.4330127018922193*v1[8]-0.25*v1[7]+0.4330127018922193*v1[6]-0.4330127018922193*v1[5]+0.25*v1[4]-0.25*v1[3]+0.25*v1[2]-0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*v1[15]-0.25*v1[14]-0.4330127018922193*v1[13]+0.4330127018922193*(v1[12]+v1[11])+0.25*v1[10]-0.25*v1[9]-0.4330127018922193*v1[8]-0.25*v1[7]-0.4330127018922193*v1[6]+0.4330127018922193*v1[5]+0.25*(v1[4]+v1[3])-0.25*v1[2]-0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*v1[15])+0.25*v1[14]-0.4330127018922193*(v1[13]+v1[12]+v1[11])+0.25*(v1[10]+v1[9])-0.4330127018922193*v1[8]+0.25*v1[7]-0.4330127018922193*(v1[6]+v1[5])+0.25*(v1[4]+v1[3]+v1[2])-0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*((-0.4330127018922193*v1[15])-0.25*v1[14]+0.4330127018922193*(v1[13]+v1[12]+v1[11])+0.25*(v1[10]+v1[9])-0.4330127018922193*v1[8]+0.25*v1[7]-0.4330127018922193*(v1[6]+v1[5])-0.25*(v1[4]+v1[3]+v1[2])+0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*v1[15]+0.25*v1[14]+0.4330127018922193*v1[13]-0.4330127018922193*(v1[12]+v1[11])+0.25*v1[10]-0.25*v1[9]-0.4330127018922193*v1[8]-0.25*v1[7]-0.4330127018922193*v1[6]+0.4330127018922193*v1[5]-0.25*(v1[4]+v1[3])+0.25*v1[2]+0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*v1[15]+0.25*v1[14]-0.4330127018922193*v1[13]+0.4330127018922193*v1[12]-0.4330127018922193*v1[11]-0.25*v1[10]+0.25*v1[9]-0.4330127018922193*v1[8]-0.25*v1[7]+0.4330127018922193*v1[6]-0.4330127018922193*v1[5]-0.25*v1[4]+0.25*v1[3]-0.25*v1[2]+0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*v1[15])-0.25*v1[14]-0.4330127018922193*(v1[13]+v1[12])+0.4330127018922193*v1[11]-0.25*(v1[10]+v1[9])-0.4330127018922193*v1[8]+0.25*v1[7]+0.4330127018922193*(v1[6]+v1[5])-0.25*v1[4]+0.25*(v1[3]+v1[2])+0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*v1[15]+0.25*v1[14]-0.4330127018922193*(v1[13]+v1[12])+0.4330127018922193*v1[11]-0.25*(v1[10]+v1[9])+0.4330127018922193*v1[8]+0.25*v1[7]-0.4330127018922193*(v1[6]+v1[5])+0.25*v1[4]-0.25*(v1[3]+v1[2])+0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*v1[15])-0.25*v1[14]-0.4330127018922193*v1[13]+0.4330127018922193*v1[12]-0.4330127018922193*v1[11]-0.25*v1[10]+0.25*v1[9]+0.4330127018922193*v1[8]-0.25*v1[7]-0.4330127018922193*v1[6]+0.4330127018922193*v1[5]+0.25*v1[4]-0.25*v1[3]+0.25*v1[2]+0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*v1[15])-0.25*v1[14]+0.4330127018922193*v1[13]-0.4330127018922193*(v1[12]+v1[11])+0.25*v1[10]-0.25*v1[9]+0.4330127018922193*v1[8]-0.25*v1[7]+0.4330127018922193*v1[6]-0.4330127018922193*v1[5]+0.25*(v1[4]+v1[3])-0.25*v1[2]+0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*v1[15]+0.25*v1[14]+0.4330127018922193*(v1[13]+v1[12]+v1[11])+0.25*(v1[10]+v1[9])+0.4330127018922193*v1[8]+0.25*v1[7]+0.4330127018922193*(v1[6]+v1[5])+0.25*(v1[4]+v1[3]+v1[2])+0.4330127018922193*v1[1]+0.25*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*v2[2]-1.0*v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*v2[2]+v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.4330127018922193*v2[15]-0.4330127018922193*v2[14]-0.25*v2[13]-0.4330127018922193*(v2[12]+v2[11])+0.25*v2[10]+0.4330127018922193*v2[9]+0.25*v2[8]+0.4330127018922193*v2[7]+0.25*v2[6]+0.4330127018922193*v2[5]-0.25*(v2[4]+v2[3])-0.4330127018922193*v2[2]-0.25*v2[1]+0.25*v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(v2[15]+v2[14]))+0.25*v2[13]+0.4330127018922193*(v2[12]+v2[11])+0.25*v2[10]+0.4330127018922193*v2[9]-0.25*v2[8]+0.4330127018922193*v2[7]-0.25*v2[6]-0.4330127018922193*v2[5]-0.25*(v2[4]+v2[3])-0.4330127018922193*v2[2]+0.25*(v2[1]+v2[0]))*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*v2[15])+0.4330127018922193*v2[14]+0.25*v2[13]-0.4330127018922193*v2[12]+0.4330127018922193*v2[11]-0.25*v2[10]+0.4330127018922193*v2[9]+0.25*v2[8]-0.4330127018922193*v2[7]-0.25*v2[6]+0.4330127018922193*v2[5]-0.25*v2[4]+0.25*v2[3]-0.4330127018922193*v2[2]-0.25*v2[1]+0.25*v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*(v2[15]+v2[14])-0.25*v2[13]+0.4330127018922193*v2[12]-0.4330127018922193*v2[11]-0.25*v2[10]+0.4330127018922193*v2[9]-0.25*v2[8]-0.4330127018922193*v2[7]+0.25*v2[6]-0.4330127018922193*v2[5]-0.25*v2[4]+0.25*v2[3]-0.4330127018922193*v2[2]+0.25*(v2[1]+v2[0]))*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*v2[15])+0.4330127018922193*v2[14]+0.25*v2[13]+0.4330127018922193*v2[12]-0.4330127018922193*v2[11]-0.25*v2[10]-0.4330127018922193*v2[9]-0.25*v2[8]+0.4330127018922193*v2[7]+0.25*v2[6]+0.4330127018922193*v2[5]+0.25*v2[4]-0.25*v2[3]-0.4330127018922193*v2[2]-0.25*v2[1]+0.25*v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*(v2[15]+v2[14])-0.25*v2[13]-0.4330127018922193*v2[12]+0.4330127018922193*v2[11]-0.25*v2[10]-0.4330127018922193*v2[9]+0.25*v2[8]+0.4330127018922193*v2[7]-0.25*v2[6]-0.4330127018922193*v2[5]+0.25*v2[4]-0.25*v2[3]-0.4330127018922193*v2[2]+0.25*(v2[1]+v2[0]))*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*v2[15]-0.4330127018922193*v2[14]-0.25*v2[13]+0.4330127018922193*(v2[12]+v2[11])+0.25*v2[10]-0.4330127018922193*v2[9]-0.25*v2[8]-0.4330127018922193*v2[7]-0.25*v2[6]+0.4330127018922193*v2[5]+0.25*(v2[4]+v2[3])-0.4330127018922193*v2[2]-0.25*v2[1]+0.25*v2[0])*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(v2[15]+v2[14]))+0.25*v2[13]-0.4330127018922193*(v2[12]+v2[11])+0.25*v2[10]-0.4330127018922193*v2[9]+0.25*v2[8]-0.4330127018922193*v2[7]+0.25*v2[6]-0.4330127018922193*v2[5]+0.25*(v2[4]+v2[3])-0.4330127018922193*v2[2]+0.25*(v2[1]+v2[0]))*dfac2; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*((-0.4330127018922193*v2[15])+0.4330127018922193*v2[14]-0.25*v2[13]+0.4330127018922193*(v2[12]+v2[11])+0.25*v2[10]-0.4330127018922193*v2[9]+0.25*v2[8]-0.4330127018922193*v2[7]+0.25*v2[6]-0.4330127018922193*v2[5]-0.25*(v2[4]+v2[3])+0.4330127018922193*v2[2]-0.25*v2[1]+0.25*v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(v2[15]+v2[14])+0.25*v2[13]-0.4330127018922193*(v2[12]+v2[11])+0.25*v2[10]-0.4330127018922193*v2[9]-0.25*v2[8]-0.4330127018922193*v2[7]-0.25*v2[6]+0.4330127018922193*v2[5]-0.25*(v2[4]+v2[3])+0.4330127018922193*v2[2]+0.25*(v2[1]+v2[0]))*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*v2[15]-0.4330127018922193*v2[14]+0.25*v2[13]+0.4330127018922193*v2[12]-0.4330127018922193*v2[11]-0.25*v2[10]-0.4330127018922193*v2[9]+0.25*v2[8]+0.4330127018922193*v2[7]-0.25*v2[6]-0.4330127018922193*v2[5]-0.25*v2[4]+0.25*v2[3]+0.4330127018922193*v2[2]-0.25*v2[1]+0.25*v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*(v2[15]+v2[14]))-0.25*v2[13]-0.4330127018922193*v2[12]+0.4330127018922193*v2[11]-0.25*v2[10]-0.4330127018922193*v2[9]-0.25*v2[8]+0.4330127018922193*v2[7]+0.25*v2[6]+0.4330127018922193*v2[5]-0.25*v2[4]+0.25*v2[3]+0.4330127018922193*v2[2]+0.25*(v2[1]+v2[0]))*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*v2[15]-0.4330127018922193*v2[14]+0.25*v2[13]-0.4330127018922193*v2[12]+0.4330127018922193*v2[11]-0.25*v2[10]+0.4330127018922193*v2[9]-0.25*v2[8]-0.4330127018922193*v2[7]+0.25*v2[6]-0.4330127018922193*v2[5]+0.25*v2[4]-0.25*v2[3]+0.4330127018922193*v2[2]-0.25*v2[1]+0.25*v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*(v2[15]+v2[14]))-0.25*v2[13]+0.4330127018922193*v2[12]-0.4330127018922193*v2[11]-0.25*v2[10]+0.4330127018922193*v2[9]+0.25*v2[8]-0.4330127018922193*v2[7]-0.25*v2[6]+0.4330127018922193*v2[5]+0.25*v2[4]-0.25*v2[3]+0.4330127018922193*v2[2]+0.25*(v2[1]+v2[0]))*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*v2[15])+0.4330127018922193*v2[14]-0.25*v2[13]-0.4330127018922193*(v2[12]+v2[11])+0.25*v2[10]+0.4330127018922193*v2[9]-0.25*v2[8]+0.4330127018922193*v2[7]-0.25*v2[6]-0.4330127018922193*v2[5]+0.25*(v2[4]+v2[3])+0.4330127018922193*v2[2]-0.25*v2[1]+0.25*v2[0])*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(v2[15]+v2[14])+0.25*v2[13]+0.4330127018922193*(v2[12]+v2[11])+0.25*v2[10]+0.4330127018922193*v2[9]+0.25*v2[8]+0.4330127018922193*v2[7]+0.25*v2[6]+0.4330127018922193*v2[5]+0.25*(v2[4]+v2[3])+0.4330127018922193*v2[2]+0.25*(v2[1]+v2[0]))*dfac2; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*v3[3]-1.0*v3[0])*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*v3[3]+v3[0])*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.4330127018922193*v3[15]-0.4330127018922193*(v3[14]+v3[13])-0.25*v3[12]-0.4330127018922193*v3[11]+0.4330127018922193*v3[10]+0.25*(v3[9]+v3[8])+0.4330127018922193*(v3[7]+v3[6])+0.25*v3[5]-0.25*v3[4]-0.4330127018922193*v3[3]-0.25*(v3[2]+v3[1])+0.25*v3[0])*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(v3[15]+v3[14]))+0.4330127018922193*v3[13]+0.25*v3[12]+0.4330127018922193*(v3[11]+v3[10])+0.25*v3[9]-0.25*v3[8]+0.4330127018922193*v3[7]-0.4330127018922193*v3[6]-0.25*(v3[5]+v3[4])-0.4330127018922193*v3[3]-0.25*v3[2]+0.25*(v3[1]+v3[0]))*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*v3[15])+0.4330127018922193*v3[14]-0.4330127018922193*v3[13]+0.25*v3[12]+0.4330127018922193*(v3[11]+v3[10])-0.25*v3[9]+0.25*v3[8]-0.4330127018922193*v3[7]+0.4330127018922193*v3[6]-0.25*(v3[5]+v3[4])-0.4330127018922193*v3[3]+0.25*v3[2]-0.25*v3[1]+0.25*v3[0])*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*(v3[15]+v3[14]+v3[13])-0.25*v3[12]-0.4330127018922193*v3[11]+0.4330127018922193*v3[10]-0.25*(v3[9]+v3[8])-0.4330127018922193*(v3[7]+v3[6])+0.25*v3[5]-0.25*v3[4]-0.4330127018922193*v3[3]+0.25*(v3[2]+v3[1]+v3[0]))*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*v3[15])+0.4330127018922193*(v3[14]+v3[13])+0.25*v3[12]-0.4330127018922193*(v3[11]+v3[10])-0.25*(v3[9]+v3[8])+0.4330127018922193*(v3[7]+v3[6])+0.25*(v3[5]+v3[4])-0.4330127018922193*v3[3]-0.25*(v3[2]+v3[1])+0.25*v3[0])*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*(v3[15]+v3[14])-0.4330127018922193*v3[13]-0.25*v3[12]+0.4330127018922193*v3[11]-0.4330127018922193*v3[10]-0.25*v3[9]+0.25*v3[8]+0.4330127018922193*v3[7]-0.4330127018922193*v3[6]-0.25*v3[5]+0.25*v3[4]-0.4330127018922193*v3[3]-0.25*v3[2]+0.25*(v3[1]+v3[0]))*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*v3[15]-0.4330127018922193*v3[14]+0.4330127018922193*v3[13]-0.25*v3[12]+0.4330127018922193*v3[11]-0.4330127018922193*v3[10]+0.25*v3[9]-0.25*v3[8]-0.4330127018922193*v3[7]+0.4330127018922193*v3[6]-0.25*v3[5]+0.25*v3[4]-0.4330127018922193*v3[3]+0.25*v3[2]-0.25*v3[1]+0.25*v3[0])*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(v3[15]+v3[14]+v3[13]))+0.25*v3[12]-0.4330127018922193*(v3[11]+v3[10])+0.25*(v3[9]+v3[8])-0.4330127018922193*(v3[7]+v3[6])+0.25*(v3[5]+v3[4])-0.4330127018922193*v3[3]+0.25*(v3[2]+v3[1]+v3[0]))*dfac3; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*((-0.4330127018922193*v3[15])+0.4330127018922193*(v3[14]+v3[13])-0.25*v3[12]+0.4330127018922193*v3[11]-0.4330127018922193*v3[10]+0.25*(v3[9]+v3[8])-0.4330127018922193*(v3[7]+v3[6])+0.25*v3[5]-0.25*v3[4]+0.4330127018922193*v3[3]-0.25*(v3[2]+v3[1])+0.25*v3[0])*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(v3[15]+v3[14])-0.4330127018922193*v3[13]+0.25*v3[12]-0.4330127018922193*(v3[11]+v3[10])+0.25*v3[9]-0.25*v3[8]-0.4330127018922193*v3[7]+0.4330127018922193*v3[6]-0.25*(v3[5]+v3[4])+0.4330127018922193*v3[3]-0.25*v3[2]+0.25*(v3[1]+v3[0]))*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*v3[15]-0.4330127018922193*v3[14]+0.4330127018922193*v3[13]+0.25*v3[12]-0.4330127018922193*(v3[11]+v3[10])-0.25*v3[9]+0.25*v3[8]+0.4330127018922193*v3[7]-0.4330127018922193*v3[6]-0.25*(v3[5]+v3[4])+0.4330127018922193*v3[3]+0.25*v3[2]-0.25*v3[1]+0.25*v3[0])*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*(v3[15]+v3[14]+v3[13]))-0.25*v3[12]+0.4330127018922193*v3[11]-0.4330127018922193*v3[10]-0.25*(v3[9]+v3[8])+0.4330127018922193*(v3[7]+v3[6])+0.25*v3[5]-0.25*v3[4]+0.4330127018922193*v3[3]+0.25*(v3[2]+v3[1]+v3[0]))*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*v3[15]-0.4330127018922193*(v3[14]+v3[13])+0.25*v3[12]+0.4330127018922193*(v3[11]+v3[10])-0.25*(v3[9]+v3[8])-0.4330127018922193*(v3[7]+v3[6])+0.25*(v3[5]+v3[4])+0.4330127018922193*v3[3]-0.25*(v3[2]+v3[1])+0.25*v3[0])*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*(v3[15]+v3[14]))+0.4330127018922193*v3[13]-0.25*v3[12]-0.4330127018922193*v3[11]+0.4330127018922193*v3[10]-0.25*v3[9]+0.25*v3[8]-0.4330127018922193*v3[7]+0.4330127018922193*v3[6]-0.25*v3[5]+0.25*v3[4]+0.4330127018922193*v3[3]-0.25*v3[2]+0.25*(v3[1]+v3[0]))*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*v3[15])+0.4330127018922193*v3[14]-0.4330127018922193*v3[13]-0.25*v3[12]-0.4330127018922193*v3[11]+0.4330127018922193*v3[10]+0.25*v3[9]-0.25*v3[8]+0.4330127018922193*v3[7]-0.4330127018922193*v3[6]-0.25*v3[5]+0.25*v3[4]+0.4330127018922193*v3[3]+0.25*v3[2]-0.25*v3[1]+0.25*v3[0])*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(v3[15]+v3[14]+v3[13])+0.25*v3[12]+0.4330127018922193*(v3[11]+v3[10])+0.25*(v3[9]+v3[8])+0.4330127018922193*(v3[7]+v3[6])+0.25*(v3[5]+v3[4])+0.4330127018922193*v3[3]+0.25*(v3[2]+v3[1]+v3[0]))*dfac3; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*v4[4]-1.0*v4[0])*dfac4; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*v4[4]+v4[0])*dfac4; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.4330127018922193*v4[15]-0.4330127018922193*(v4[14]+v4[13]+v4[12])-0.25*v4[11]+0.4330127018922193*(v4[10]+v4[9]+v4[8])+0.25*(v4[7]+v4[6]+v4[5])-0.4330127018922193*v4[4]-0.25*(v4[3]+v4[2]+v4[1])+0.25*v4[0])*dfac4; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(v4[15]+v4[14]))+0.4330127018922193*(v4[13]+v4[12])+0.25*v4[11]+0.4330127018922193*(v4[10]+v4[9])-0.4330127018922193*v4[8]+0.25*v4[7]-0.25*(v4[6]+v4[5])-0.4330127018922193*v4[4]-0.25*(v4[3]+v4[2])+0.25*(v4[1]+v4[0]))*dfac4; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*v4[15])+0.4330127018922193*v4[14]-0.4330127018922193*v4[13]+0.4330127018922193*v4[12]+0.25*v4[11]+0.4330127018922193*v4[10]-0.4330127018922193*v4[9]+0.4330127018922193*v4[8]-0.25*v4[7]+0.25*v4[6]-0.25*v4[5]-0.4330127018922193*v4[4]-0.25*v4[3]+0.25*v4[2]-0.25*v4[1]+0.25*v4[0])*dfac4; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*(v4[15]+v4[14]+v4[13])-0.4330127018922193*v4[12]-0.25*v4[11]+0.4330127018922193*v4[10]-0.4330127018922193*(v4[9]+v4[8])-0.25*(v4[7]+v4[6])+0.25*v4[5]-0.4330127018922193*v4[4]-0.25*v4[3]+0.25*(v4[2]+v4[1]+v4[0]))*dfac4; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*v4[15])+0.4330127018922193*(v4[14]+v4[13])-0.4330127018922193*v4[12]+0.25*v4[11]-0.4330127018922193*v4[10]+0.4330127018922193*(v4[9]+v4[8])-0.25*(v4[7]+v4[6])+0.25*v4[5]-0.4330127018922193*v4[4]+0.25*v4[3]-0.25*(v4[2]+v4[1])+0.25*v4[0])*dfac4; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*(v4[15]+v4[14])-0.4330127018922193*v4[13]+0.4330127018922193*v4[12]-0.25*v4[11]-0.4330127018922193*v4[10]+0.4330127018922193*v4[9]-0.4330127018922193*v4[8]-0.25*v4[7]+0.25*v4[6]-0.25*v4[5]-0.4330127018922193*v4[4]+0.25*v4[3]-0.25*v4[2]+0.25*(v4[1]+v4[0]))*dfac4; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*v4[15]-0.4330127018922193*v4[14]+0.4330127018922193*(v4[13]+v4[12])-0.25*v4[11]-0.4330127018922193*(v4[10]+v4[9])+0.4330127018922193*v4[8]+0.25*v4[7]-0.25*(v4[6]+v4[5])-0.4330127018922193*v4[4]+0.25*(v4[3]+v4[2])-0.25*v4[1]+0.25*v4[0])*dfac4; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(v4[15]+v4[14]+v4[13]+v4[12]))+0.25*v4[11]-0.4330127018922193*(v4[10]+v4[9]+v4[8])+0.25*(v4[7]+v4[6]+v4[5])-0.4330127018922193*v4[4]+0.25*(v4[3]+v4[2]+v4[1]+v4[0]))*dfac4; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*((-0.4330127018922193*v4[15])+0.4330127018922193*(v4[14]+v4[13]+v4[12])-0.25*v4[11]-0.4330127018922193*(v4[10]+v4[9]+v4[8])+0.25*(v4[7]+v4[6]+v4[5])+0.4330127018922193*v4[4]-0.25*(v4[3]+v4[2]+v4[1])+0.25*v4[0])*dfac4; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(v4[15]+v4[14])-0.4330127018922193*(v4[13]+v4[12])+0.25*v4[11]-0.4330127018922193*(v4[10]+v4[9])+0.4330127018922193*v4[8]+0.25*v4[7]-0.25*(v4[6]+v4[5])+0.4330127018922193*v4[4]-0.25*(v4[3]+v4[2])+0.25*(v4[1]+v4[0]))*dfac4; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*v4[15]-0.4330127018922193*v4[14]+0.4330127018922193*v4[13]-0.4330127018922193*v4[12]+0.25*v4[11]-0.4330127018922193*v4[10]+0.4330127018922193*v4[9]-0.4330127018922193*v4[8]-0.25*v4[7]+0.25*v4[6]-0.25*v4[5]+0.4330127018922193*v4[4]-0.25*v4[3]+0.25*v4[2]-0.25*v4[1]+0.25*v4[0])*dfac4; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*(v4[15]+v4[14]+v4[13]))+0.4330127018922193*v4[12]-0.25*v4[11]-0.4330127018922193*v4[10]+0.4330127018922193*(v4[9]+v4[8])-0.25*(v4[7]+v4[6])+0.25*v4[5]+0.4330127018922193*v4[4]-0.25*v4[3]+0.25*(v4[2]+v4[1]+v4[0]))*dfac4; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*v4[15]-0.4330127018922193*(v4[14]+v4[13])+0.4330127018922193*v4[12]+0.25*v4[11]+0.4330127018922193*v4[10]-0.4330127018922193*(v4[9]+v4[8])-0.25*(v4[7]+v4[6])+0.25*v4[5]+0.4330127018922193*v4[4]+0.25*v4[3]-0.25*(v4[2]+v4[1])+0.25*v4[0])*dfac4; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*(v4[15]+v4[14]))+0.4330127018922193*v4[13]-0.4330127018922193*v4[12]-0.25*v4[11]+0.4330127018922193*v4[10]-0.4330127018922193*v4[9]+0.4330127018922193*v4[8]-0.25*v4[7]+0.25*v4[6]-0.25*v4[5]+0.4330127018922193*v4[4]+0.25*v4[3]-0.25*v4[2]+0.25*(v4[1]+v4[0]))*dfac4; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*v4[15])+0.4330127018922193*v4[14]-0.4330127018922193*(v4[13]+v4[12])-0.25*v4[11]+0.4330127018922193*(v4[10]+v4[9])-0.4330127018922193*v4[8]+0.25*v4[7]-0.25*(v4[6]+v4[5])+0.4330127018922193*v4[4]+0.25*(v4[3]+v4[2])-0.25*v4[1]+0.25*v4[0])*dfac4; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(v4[15]+v4[14]+v4[13]+v4[12])+0.25*v4[11]+0.4330127018922193*(v4[10]+v4[9]+v4[8])+0.25*(v4[7]+v4[6]+v4[5])+0.4330127018922193*v4[4]+0.25*(v4[3]+v4[2]+v4[1]+v4[0]))*dfac4; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.4330127018922193*(f[15]*v1[15]+f[14]*v1[14]+f[13]*v1[13]+f[12]*v1[12]+f[11]*v1[11]+f[10]*v1[10]+f[9]*v1[9]+f[8]*v1[8]+f[7]*v1[7]+f[6]*v1[6]+f[5]*v1[5]+f[4]*v1[4]+f[3]*v1[3]+f[2]*v1[2]+f[1]*v1[1]+f[0]*v1[0])*dfac1; 
  out[2] += 0.4330127018922193*(f[15]*v2[15]+f[14]*v2[14]+f[13]*v2[13]+f[12]*v2[12]+f[11]*v2[11]+f[10]*v2[10]+f[9]*v2[9]+f[8]*v2[8]+f[7]*v2[7]+f[6]*v2[6]+f[5]*v2[5]+f[4]*v2[4]+f[3]*v2[3]+f[2]*v2[2]+f[1]*v2[1]+f[0]*v2[0])*dfac2; 
  out[3] += 0.4330127018922193*(f[15]*v3[15]+f[14]*v3[14]+f[13]*v3[13]+f[12]*v3[12]+f[11]*v3[11]+f[10]*v3[10]+f[9]*v3[9]+f[8]*v3[8]+f[7]*v3[7]+f[6]*v3[6]+f[5]*v3[5]+f[4]*v3[4]+f[3]*v3[3]+f[2]*v3[2]+f[1]*v3[1]+f[0]*v3[0])*dfac3; 
  out[4] += 0.4330127018922193*(f[15]*v4[15]+f[14]*v4[14]+f[13]*v4[13]+f[12]*v4[12]+f[11]*v4[11]+f[10]*v4[10]+f[9]*v4[9]+f[8]*v4[8]+f[7]*v4[7]+f[6]*v4[6]+f[5]*v4[5]+f[4]*v4[4]+f[3]*v4[3]+f[2]*v4[2]+f[1]*v4[1]+f[0]*v4[0])*dfac4; 
  out[5] += 0.4330127018922193*((f[14]*v2[15]+v2[14]*f[15]+f[10]*v2[13]+v2[10]*f[13]+f[9]*v2[12]+v2[9]*f[12]+f[7]*v2[11]+v2[7]*f[11]+f[4]*v2[8]+v2[4]*f[8]+f[3]*v2[6]+v2[3]*f[6]+f[2]*v2[5]+v2[2]*f[5]+f[0]*v2[1]+v2[0]*f[1])*dfac2+(f[13]*v1[15]+v1[13]*f[15]+f[10]*v1[14]+v1[10]*f[14]+f[8]*v1[12]+v1[8]*f[12]+f[6]*v1[11]+v1[6]*f[11]+f[4]*v1[9]+v1[4]*f[9]+f[3]*v1[7]+v1[3]*f[7]+f[1]*v1[5]+v1[1]*f[5]+f[0]*v1[2]+v1[0]*f[2])*dfac1); 
  out[6] += 0.4330127018922193*((f[14]*v3[15]+v3[14]*f[15]+f[10]*v3[13]+v3[10]*f[13]+f[9]*v3[12]+v3[9]*f[12]+f[7]*v3[11]+v3[7]*f[11]+f[4]*v3[8]+v3[4]*f[8]+f[3]*v3[6]+v3[3]*f[6]+f[2]*v3[5]+v3[2]*f[5]+f[0]*v3[1]+v3[0]*f[1])*dfac3+(f[12]*v1[15]+v1[12]*f[15]+f[9]*v1[14]+v1[9]*f[14]+f[8]*v1[13]+v1[8]*f[13]+f[5]*v1[11]+v1[5]*f[11]+f[4]*v1[10]+v1[4]*f[10]+f[2]*v1[7]+v1[2]*f[7]+f[1]*v1[6]+v1[1]*f[6]+f[0]*v1[3]+v1[0]*f[3])*dfac1); 
  out[7] += 0.4330127018922193*((f[13]*v3[15]+v3[13]*f[15]+f[10]*v3[14]+v3[10]*f[14]+f[8]*v3[12]+v3[8]*f[12]+f[6]*v3[11]+v3[6]*f[11]+f[4]*v3[9]+v3[4]*f[9]+f[3]*v3[7]+v3[3]*f[7]+f[1]*v3[5]+v3[1]*f[5]+f[0]*v3[2]+v3[0]*f[2])*dfac3+(f[12]*v2[15]+v2[12]*f[15]+f[9]*v2[14]+v2[9]*f[14]+f[8]*v2[13]+v2[8]*f[13]+f[5]*v2[11]+v2[5]*f[11]+f[4]*v2[10]+v2[4]*f[10]+f[2]*v2[7]+v2[2]*f[7]+f[1]*v2[6]+v2[1]*f[6]+f[0]*v2[3]+v2[0]*f[3])*dfac2); 
  out[8] += 0.4330127018922193*((f[14]*v4[15]+v4[14]*f[15]+f[10]*v4[13]+v4[10]*f[13]+f[9]*v4[12]+v4[9]*f[12]+f[7]*v4[11]+v4[7]*f[11]+f[4]*v4[8]+v4[4]*f[8]+f[3]*v4[6]+v4[3]*f[6]+f[2]*v4[5]+v4[2]*f[5]+f[0]*v4[1]+v4[0]*f[1])*dfac4+(f[11]*v1[15]+v1[11]*f[15]+f[7]*v1[14]+v1[7]*f[14]+f[6]*v1[13]+v1[6]*f[13]+f[5]*v1[12]+v1[5]*f[12]+f[3]*v1[10]+v1[3]*f[10]+f[2]*v1[9]+v1[2]*f[9]+f[1]*v1[8]+v1[1]*f[8]+f[0]*v1[4]+v1[0]*f[4])*dfac1); 
  out[9] += 0.4330127018922193*((f[13]*v4[15]+v4[13]*f[15]+f[10]*v4[14]+v4[10]*f[14]+f[8]*v4[12]+v4[8]*f[12]+f[6]*v4[11]+v4[6]*f[11]+f[4]*v4[9]+v4[4]*f[9]+f[3]*v4[7]+v4[3]*f[7]+f[1]*v4[5]+v4[1]*f[5]+f[0]*v4[2]+v4[0]*f[2])*dfac4+(f[11]*v2[15]+v2[11]*f[15]+f[7]*v2[14]+v2[7]*f[14]+f[6]*v2[13]+v2[6]*f[13]+f[5]*v2[12]+v2[5]*f[12]+f[3]*v2[10]+v2[3]*f[10]+f[2]*v2[9]+v2[2]*f[9]+f[1]*v2[8]+v2[1]*f[8]+f[0]*v2[4]+v2[0]*f[4])*dfac2); 
  out[10] += 0.4330127018922193*((f[12]*v4[15]+v4[12]*f[15]+f[9]*v4[14]+v4[9]*f[14]+f[8]*v4[13]+v4[8]*f[13]+f[5]*v4[11]+v4[5]*f[11]+f[4]*v4[10]+v4[4]*f[10]+f[2]*v4[7]+v4[2]*f[7]+f[1]*v4[6]+v4[1]*f[6]+f[0]*v4[3]+v4[0]*f[3])*dfac4+(f[11]*v3[15]+v3[11]*f[15]+f[7]*v3[14]+v3[7]*f[14]+f[6]*v3[13]+v3[6]*f[13]+f[5]*v3[12]+v3[5]*f[12]+f[3]*v3[10]+v3[3]*f[10]+f[2]*v3[9]+v3[2]*f[9]+f[1]*v3[8]+v3[1]*f[8]+f[0]*v3[4]+v3[0]*f[4])*dfac3); 
  out[11] += 0.4330127018922193*((f[10]*v3[15]+v3[10]*f[15]+f[13]*v3[14]+v3[13]*f[14]+f[4]*v3[12]+v3[4]*f[12]+f[3]*v3[11]+v3[3]*f[11]+f[8]*v3[9]+v3[8]*f[9]+f[6]*v3[7]+v3[6]*f[7]+f[0]*v3[5]+v3[0]*f[5]+f[1]*v3[2]+v3[1]*f[2])*dfac3+(f[9]*v2[15]+v2[9]*f[15]+f[12]*v2[14]+v2[12]*f[14]+f[4]*v2[13]+v2[4]*f[13]+f[2]*v2[11]+v2[2]*f[11]+f[8]*v2[10]+v2[8]*f[10]+f[5]*v2[7]+v2[5]*f[7]+f[0]*v2[6]+v2[0]*f[6]+f[1]*v2[3]+v2[1]*f[3])*dfac2+(f[8]*v1[15]+v1[8]*f[15]+f[4]*v1[14]+v1[4]*f[14]+f[12]*v1[13]+v1[12]*f[13]+f[1]*v1[11]+v1[1]*f[11]+f[9]*v1[10]+v1[9]*f[10]+f[0]*v1[7]+v1[0]*f[7]+f[5]*v1[6]+v1[5]*f[6]+f[2]*v1[3]+v1[2]*f[3])*dfac1); 
  out[12] += 0.4330127018922193*((f[10]*v4[15]+v4[10]*f[15]+f[13]*v4[14]+v4[13]*f[14]+f[4]*v4[12]+v4[4]*f[12]+f[3]*v4[11]+v4[3]*f[11]+f[8]*v4[9]+v4[8]*f[9]+f[6]*v4[7]+v4[6]*f[7]+f[0]*v4[5]+v4[0]*f[5]+f[1]*v4[2]+v4[1]*f[2])*dfac4+(f[7]*v2[15]+v2[7]*f[15]+f[11]*v2[14]+v2[11]*f[14]+f[3]*v2[13]+v2[3]*f[13]+f[2]*v2[12]+v2[2]*f[12]+f[6]*v2[10]+v2[6]*f[10]+f[5]*v2[9]+v2[5]*f[9]+f[0]*v2[8]+v2[0]*f[8]+f[1]*v2[4]+v2[1]*f[4])*dfac2+(f[6]*v1[15]+v1[6]*f[15]+f[3]*v1[14]+v1[3]*f[14]+f[11]*v1[13]+v1[11]*f[13]+f[1]*v1[12]+v1[1]*f[12]+f[7]*v1[10]+v1[7]*f[10]+f[0]*v1[9]+v1[0]*f[9]+f[5]*v1[8]+v1[5]*f[8]+f[2]*v1[4]+v1[2]*f[4])*dfac1); 
  out[13] += 0.4330127018922193*((f[9]*v4[15]+v4[9]*f[15]+f[12]*v4[14]+v4[12]*f[14]+f[4]*v4[13]+v4[4]*f[13]+f[2]*v4[11]+v4[2]*f[11]+f[8]*v4[10]+v4[8]*f[10]+f[5]*v4[7]+v4[5]*f[7]+f[0]*v4[6]+v4[0]*f[6]+f[1]*v4[3]+v4[1]*f[3])*dfac4+(f[7]*v3[15]+v3[7]*f[15]+f[11]*v3[14]+v3[11]*f[14]+f[3]*v3[13]+v3[3]*f[13]+f[2]*v3[12]+v3[2]*f[12]+f[6]*v3[10]+v3[6]*f[10]+f[5]*v3[9]+v3[5]*f[9]+f[0]*v3[8]+v3[0]*f[8]+f[1]*v3[4]+v3[1]*f[4])*dfac3+(f[5]*v1[15]+v1[5]*f[15]+f[2]*v1[14]+v1[2]*f[14]+f[1]*v1[13]+v1[1]*f[13]+f[11]*v1[12]+v1[11]*f[12]+f[0]*v1[10]+v1[0]*f[10]+f[7]*v1[9]+v1[7]*f[9]+f[6]*v1[8]+v1[6]*f[8]+f[3]*v1[4]+v1[3]*f[4])*dfac1); 
  out[14] += 0.4330127018922193*((f[8]*v4[15]+v4[8]*f[15]+f[4]*v4[14]+v4[4]*f[14]+f[12]*v4[13]+v4[12]*f[13]+f[1]*v4[11]+v4[1]*f[11]+f[9]*v4[10]+v4[9]*f[10]+f[0]*v4[7]+v4[0]*f[7]+f[5]*v4[6]+v4[5]*f[6]+f[2]*v4[3]+v4[2]*f[3])*dfac4+(f[6]*v3[15]+v3[6]*f[15]+f[3]*v3[14]+v3[3]*f[14]+f[11]*v3[13]+v3[11]*f[13]+f[1]*v3[12]+v3[1]*f[12]+f[7]*v3[10]+v3[7]*f[10]+f[0]*v3[9]+v3[0]*f[9]+f[5]*v3[8]+v3[5]*f[8]+f[2]*v3[4]+v3[2]*f[4])*dfac3+(f[5]*v2[15]+v2[5]*f[15]+f[2]*v2[14]+v2[2]*f[14]+f[1]*v2[13]+v2[1]*f[13]+f[11]*v2[12]+v2[11]*f[12]+f[0]*v2[10]+v2[0]*f[10]+f[7]*v2[9]+v2[7]*f[9]+f[6]*v2[8]+v2[6]*f[8]+f[3]*v2[4]+v2[3]*f[4])*dfac2); 
  out[15] += 0.4330127018922193*((f[4]*v4[15]+v4[4]*f[15]+f[8]*v4[14]+v4[8]*f[14]+f[9]*v4[13]+v4[9]*f[13]+f[10]*v4[12]+v4[10]*f[12]+f[0]*v4[11]+v4[0]*f[11]+f[1]*v4[7]+v4[1]*f[7]+f[2]*v4[6]+v4[2]*f[6]+f[3]*v4[5]+v4[3]*f[5])*dfac4+(f[3]*v3[15]+v3[3]*f[15]+f[6]*v3[14]+v3[6]*f[14]+f[7]*v3[13]+v3[7]*f[13]+f[0]*v3[12]+v3[0]*f[12]+f[10]*v3[11]+v3[10]*f[11]+f[1]*v3[9]+v3[1]*f[9]+f[2]*v3[8]+v3[2]*f[8]+f[4]*v3[5]+v3[4]*f[5])*dfac3+(f[2]*v2[15]+v2[2]*f[15]+f[5]*v2[14]+v2[5]*f[14]+f[0]*v2[13]+v2[0]*f[13]+f[7]*v2[12]+v2[7]*f[12]+f[9]*v2[11]+v2[9]*f[11]+f[1]*v2[10]+v2[1]*f[10]+f[3]*v2[8]+v2[3]*f[8]+f[4]*v2[6]+v2[4]*f[6])*dfac2+(f[1]*v1[15]+v1[1]*f[15]+f[0]*v1[14]+v1[0]*f[14]+f[5]*v1[13]+v1[5]*f[13]+f[6]*v1[12]+v1[6]*f[12]+f[8]*v1[11]+v1[8]*f[11]+f[2]*v1[10]+v1[2]*f[10]+f[3]*v1[9]+v1[3]*f[9]+f[4]*v1[7]+v1[4]*f[7])*dfac1); 
  return cflRate; 
} 
