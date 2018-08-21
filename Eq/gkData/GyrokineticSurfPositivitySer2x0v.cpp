#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity2x0vSer_X_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.25*(3.0*BmagInv[0]*Phi[3]-1.732050807568877*BmagInv[0]*Phi[2])*dfac_y; 

  double alpha[4]; 
  alpha[0] = -0.8660254037844386*BmagInv[0]*Phi[2]*dfac_y; 
  alpha[1] = -0.8660254037844386*BmagInv[0]*Phi[3]*dfac_y; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[2] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cfl); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*(3.0*alpha[1]*fhatALVal[1]-1.732050807568877*alpha[0]*fhatALVal[1]+1.732050807568877*fhatALVal[0]*alpha[1]-1.0*alpha[0]*fhatALVal[0])*dfac_x; 
  incr[1] = 0.25*(5.196152422706631*alpha[1]*fhatALVal[1]-3.0*alpha[0]*fhatALVal[1]+3.0*fhatALVal[0]*alpha[1]-1.732050807568877*alpha[0]*fhatALVal[0])*dfac_x; 
  incr[2] = -0.25*(3.0*alpha[1]*fhatALVal[3]-1.732050807568877*alpha[0]*fhatALVal[3]+1.732050807568877*alpha[1]*fhatALVal[2]-1.0*alpha[0]*fhatALVal[2])*dfac_x; 
  incr[3] = 0.25*(5.196152422706631*alpha[1]*fhatALVal[3]-3.0*alpha[0]*fhatALVal[3]+3.0*alpha[1]*fhatALVal[2]-1.732050807568877*alpha[0]*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[2] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cfl); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*(3.0*alpha[1]*fhatALVal[1]-1.732050807568877*alpha[0]*fhatALVal[1]-1.732050807568877*fhatALVal[0]*alpha[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[1] = -0.25*(5.196152422706631*alpha[1]*fhatALVal[1]-3.0*alpha[0]*fhatALVal[1]-3.0*fhatALVal[0]*alpha[1]+1.732050807568877*alpha[0]*fhatALVal[0])*dfac_x; 
  incr[2] = 0.25*(3.0*alpha[1]*fhatALVal[3]-1.732050807568877*alpha[0]*fhatALVal[3]-1.732050807568877*alpha[1]*fhatALVal[2]+alpha[0]*fhatALVal[2])*dfac_x; 
  incr[3] = -0.25*(5.196152422706631*alpha[1]*fhatALVal[3]-3.0*alpha[0]*fhatALVal[3]-3.0*alpha[1]*fhatALVal[2]+1.732050807568877*alpha[0]*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x0vSer_Y_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.25*(3.0*BmagInv[0]*Phi[3]-1.732050807568877*BmagInv[0]*Phi[1])*dfac_x; 

  double alpha[4]; 
  alpha[0] = 0.8660254037844386*BmagInv[0]*Phi[1]*dfac_x; 
  alpha[2] = 0.8660254037844386*BmagInv[0]*Phi[3]*dfac_x; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cfl); 
  fqVal[2] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cfl); 
  fqVal[3] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cfl); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*(3.0*alpha[2]*fhatALVal[2]-1.732050807568877*alpha[0]*fhatALVal[2]+1.732050807568877*fhatALVal[0]*alpha[2]-1.0*alpha[0]*fhatALVal[0])*dfac_y; 
  incr[1] = -0.25*(3.0*alpha[2]*fhatALVal[3]-1.732050807568877*alpha[0]*fhatALVal[3]+1.732050807568877*fhatALVal[1]*alpha[2]-1.0*alpha[0]*fhatALVal[1])*dfac_y; 
  incr[2] = 0.25*(5.196152422706631*alpha[2]*fhatALVal[2]-3.0*alpha[0]*fhatALVal[2]+3.0*fhatALVal[0]*alpha[2]-1.732050807568877*alpha[0]*fhatALVal[0])*dfac_y; 
  incr[3] = 0.25*(5.196152422706631*alpha[2]*fhatALVal[3]-3.0*alpha[0]*fhatALVal[3]+3.0*fhatALVal[1]*alpha[2]-1.732050807568877*alpha[0]*fhatALVal[1])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cfl); 
  fqVal[2] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cfl); 
  fqVal[3] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cfl); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*(3.0*alpha[2]*fhatALVal[2]-1.732050807568877*alpha[0]*fhatALVal[2]-1.732050807568877*fhatALVal[0]*alpha[2]+alpha[0]*fhatALVal[0])*dfac_y; 
  incr[1] = 0.25*(3.0*alpha[2]*fhatALVal[3]-1.732050807568877*alpha[0]*fhatALVal[3]-1.732050807568877*fhatALVal[1]*alpha[2]+alpha[0]*fhatALVal[1])*dfac_y; 
  incr[2] = -0.25*(5.196152422706631*alpha[2]*fhatALVal[2]-3.0*alpha[0]*fhatALVal[2]-3.0*fhatALVal[0]*alpha[2]+1.732050807568877*alpha[0]*fhatALVal[0])*dfac_y; 
  incr[3] = -0.25*(5.196152422706631*alpha[2]*fhatALVal[3]-3.0*alpha[0]*fhatALVal[3]-3.0*fhatALVal[1]*alpha[2]+1.732050807568877*alpha[0]*fhatALVal[1])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 