#include <sqrt_on_basis_mod_decl.h>

void sqrt_on_basis_gauss_3x_p1_max(const double qExp, const double *fIn, double *out) 
{ 
  // qExp: exponent in sqrt(f)^q.
  // fIn:  input field.
  // out:  output field.
double sqrtfRq[8];
  sqrtfRq[0] = pow(sqrt((-0.3535533905932736*fIn[3])-0.3535533905932736*fIn[2]-0.3535533905932736*fIn[1]+0.3535533905932737*fIn[0]),qExp); 
  sqrtfRq[1] = pow(sqrt(0.3535533905932736*fIn[3]-0.3535533905932736*fIn[2]-0.3535533905932736*fIn[1]+0.3535533905932737*fIn[0]),qExp); 
  sqrtfRq[2] = pow(sqrt((-0.3535533905932736*fIn[3])+0.3535533905932736*fIn[2]-0.3535533905932736*fIn[1]+0.3535533905932737*fIn[0]),qExp); 
  sqrtfRq[3] = pow(sqrt(0.3535533905932736*fIn[3]+0.3535533905932736*fIn[2]-0.3535533905932736*fIn[1]+0.3535533905932737*fIn[0]),qExp); 
  sqrtfRq[4] = pow(sqrt((-0.3535533905932736*fIn[3])-0.3535533905932736*fIn[2]+0.3535533905932736*fIn[1]+0.3535533905932737*fIn[0]),qExp); 
  sqrtfRq[5] = pow(sqrt(0.3535533905932736*fIn[3]-0.3535533905932736*fIn[2]+0.3535533905932736*fIn[1]+0.3535533905932737*fIn[0]),qExp); 
  sqrtfRq[6] = pow(sqrt((-0.3535533905932736*fIn[3])+0.3535533905932736*fIn[2]+0.3535533905932736*fIn[1]+0.3535533905932737*fIn[0]),qExp); 
  sqrtfRq[7] = pow(sqrt(0.3535533905932736*fIn[3]+0.3535533905932736*fIn[2]+0.3535533905932736*fIn[1]+0.3535533905932737*fIn[0]),qExp); 

  out[0] = 0.3535533905932737*(sqrtfRq[7]+sqrtfRq[6]+sqrtfRq[5]+sqrtfRq[4]+sqrtfRq[3]+sqrtfRq[2]+sqrtfRq[1]+sqrtfRq[0]); 
  out[1] = 0.3535533905932736*(sqrtfRq[7]+sqrtfRq[6]+sqrtfRq[5]+sqrtfRq[4])-0.3535533905932736*(sqrtfRq[3]+sqrtfRq[2]+sqrtfRq[1]+sqrtfRq[0]); 
  out[2] = 0.3535533905932736*(sqrtfRq[7]+sqrtfRq[6])-0.3535533905932736*(sqrtfRq[5]+sqrtfRq[4])+0.3535533905932736*(sqrtfRq[3]+sqrtfRq[2])-0.3535533905932736*(sqrtfRq[1]+sqrtfRq[0]); 
  out[3] = 0.3535533905932736*sqrtfRq[7]-0.3535533905932736*sqrtfRq[6]+0.3535533905932736*sqrtfRq[5]-0.3535533905932736*sqrtfRq[4]+0.3535533905932736*sqrtfRq[3]-0.3535533905932736*sqrtfRq[2]+0.3535533905932736*sqrtfRq[1]-0.3535533905932736*sqrtfRq[0]; 
}
