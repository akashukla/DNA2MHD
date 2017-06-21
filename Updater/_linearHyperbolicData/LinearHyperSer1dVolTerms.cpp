#include <LinearHyperModDecl.h>

void LinearHyperSer1DP1_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[1] += 1.732050807568877*f[0]*vfact; 
  }
} 
void LinearHyperSer1DP2_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[1] += 1.732050807568877*f[0]*vfact; 
  out[2] += 3.872983346207417*f[1]*vfact; 
  }
} 
void LinearHyperSer1DP3_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[1] += 1.732050807568877*f[0]*vfact; 
  out[2] += 3.872983346207417*f[1]*vfact; 
  out[3] += 5.916079783099617*f[2]*vfact+2.645751311064591*f[0]*vfact; 
  }
} 
void LinearHyperSer1DP4_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[1] += 1.732050807568877*f[0]*vfact; 
  out[2] += 3.872983346207417*f[1]*vfact; 
  out[3] += 5.916079783099617*f[2]*vfact+2.645751311064591*f[0]*vfact; 
  out[4] += 7.937253933193772*f[3]*vfact+5.196152422706631*f[1]*vfact; 
  }
} 