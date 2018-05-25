#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x1vSerP1(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). u: flow speed. vtSq: thermal speed squared. f: Input distribution function. out: Incremented output 
  double rdv2nu[1]; 
  double rdvSq4nu[1]; 
  rdv2nu[0] = 2.0*nu/dxv[1]; 
  rdvSq4nu[0] = 4.0*nu/(dxv[1]*dxv[1]); 

  double drBar[2]; 
  drBar[0] = u[0]-1.414213562373095*w[1]; 
  drBar[1] = u[1]; 

  out[2] += rdv2nu[0]*(1.224744871391589*drBar[1]*f[1]+1.224744871391589*drBar[0]*f[0])-1.0*f[2]*nu; 
  out[3] += rdv2nu[0]*(1.224744871391589*drBar[0]*f[1]+1.224744871391589*f[0]*drBar[1])-1.0*f[3]*nu; 

return rdvSq4nu[0]*0.5; 

} 
double VmLBOconstNuVol1x1vSerP2(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). u: flow speed. vtSq: thermal speed squared. f: Input distribution function. out: Incremented output 
  double rdv2nu[1]; 
  double rdvSq4nu[1]; 
  rdv2nu[0] = 2.0*nu/dxv[1]; 
  rdvSq4nu[0] = 4.0*nu/(dxv[1]*dxv[1]); 

  double drBar[3]; 
  drBar[0] = u[0]-1.414213562373095*w[1]; 
  drBar[1] = u[1]; 
  drBar[2] = u[2]; 

  out[2] += rdv2nu[0]*(1.224744871391589*drBar[2]*f[4]+1.224744871391589*drBar[1]*f[1]+1.224744871391589*drBar[0]*f[0])-1.0*f[2]*nu; 
  out[3] += rdv2nu[0]*(1.095445115010332*drBar[1]*f[4]+1.095445115010332*f[1]*drBar[2]+1.224744871391589*drBar[0]*f[1]+1.224744871391589*f[0]*drBar[1])-1.0*f[3]*nu; 
  out[5] += (-2.0*f[5]*nu)-2.23606797749979*f[0]*nu+rdv2nu[0]*(2.738612787525831*drBar[2]*f[6]+2.738612787525831*drBar[1]*f[3]+2.738612787525831*drBar[0]*f[2])+rdvSq4nu[0]*(4.743416490252569*vtSq[2]*f[4]+4.743416490252569*f[1]*vtSq[1]+4.743416490252569*f[0]*vtSq[0]); 
  out[6] += rdv2nu[0]*(0.7824607964359517*drBar[2]*f[4]+1.224744871391589*drBar[0]*f[4]+1.224744871391589*f[0]*drBar[2]+1.095445115010332*drBar[1]*f[1])-1.0*f[6]*nu; 
  out[7] += (-2.0*f[7]*nu)-2.23606797749979*f[1]*nu+rdv2nu[0]*(2.449489742783178*drBar[1]*f[6]+2.449489742783178*drBar[2]*f[3]+2.738612787525831*drBar[0]*f[3]+2.738612787525831*drBar[1]*f[2])+rdvSq4nu[0]*(4.242640687119286*vtSq[1]*f[4]+4.242640687119286*f[1]*vtSq[2]+4.743416490252569*f[0]*vtSq[1]+4.743416490252569*vtSq[0]*f[1]); 

return rdvSq4nu[0]*0.5; 

} 
double VmLBOconstNuVol1x1vSerP3(const double *w, const double *dxv, const double nu, const double *u, const double *vtSq, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). u: flow speed. vtSq: thermal speed squared. f: Input distribution function. out: Incremented output 
  double rdv2nu[1]; 
  double rdvSq4nu[1]; 
  rdv2nu[0] = 2.0*nu/dxv[1]; 
  rdvSq4nu[0] = 4.0*nu/(dxv[1]*dxv[1]); 

  double drBar[4]; 
  drBar[0] = u[0]-1.414213562373095*w[1]; 
  drBar[1] = u[1]; 
  drBar[2] = u[2]; 
  drBar[3] = u[3]; 

  out[2] += rdv2nu[0]*(1.224744871391589*drBar[3]*f[8]+1.224744871391589*drBar[2]*f[4]+1.224744871391589*drBar[1]*f[1]+1.224744871391589*drBar[0]*f[0])-1.0*f[2]*nu; 
  out[3] += rdv2nu[0]*(1.075705748400954*drBar[2]*f[8]+1.075705748400954*drBar[3]*f[4]+1.095445115010332*drBar[1]*f[4]+1.095445115010332*f[1]*drBar[2]+1.224744871391589*drBar[0]*f[1]+1.224744871391589*f[0]*drBar[1])-1.0*f[3]*nu; 
  out[5] += (-2.0*f[5]*nu)-2.23606797749979*f[0]*nu+rdv2nu[0]*(2.73861278752583*drBar[3]*f[10]+2.738612787525831*drBar[2]*f[6]+2.738612787525831*drBar[1]*f[3]+2.738612787525831*drBar[0]*f[2])+rdvSq4nu[0]*(4.743416490252569*vtSq[3]*f[8]+4.743416490252569*vtSq[2]*f[4]+4.743416490252569*f[1]*vtSq[1]+4.743416490252569*f[0]*vtSq[0]); 
  out[6] += rdv2nu[0]*(0.7302967433402215*drBar[3]*f[8]+1.075705748400954*drBar[1]*f[8]+0.7824607964359517*drBar[2]*f[4]+1.224744871391589*drBar[0]*f[4]+1.075705748400954*f[1]*drBar[3]+1.224744871391589*f[0]*drBar[2]+1.095445115010332*drBar[1]*f[1])-1.0*f[6]*nu; 
  out[7] += (-2.0*f[7]*nu)-2.23606797749979*f[1]*nu+rdv2nu[0]*(2.405351177211819*drBar[2]*f[10]+2.405351177211819*drBar[3]*f[6]+2.449489742783178*drBar[1]*f[6]+2.449489742783178*drBar[2]*f[3]+2.738612787525831*drBar[0]*f[3]+2.738612787525831*drBar[1]*f[2])+rdvSq4nu[0]*(4.16619044897648*vtSq[2]*f[8]+4.16619044897648*vtSq[3]*f[4]+4.242640687119286*vtSq[1]*f[4]+4.242640687119286*f[1]*vtSq[2]+4.743416490252569*f[0]*vtSq[1]+4.743416490252569*vtSq[0]*f[1]); 
  out[9] += (-3.0*f[9]*nu)-4.58257569495584*f[2]*nu+rdvSq4nu[0]*(16.20185174601965*vtSq[3]*f[10]+16.20185174601965*vtSq[2]*f[6]+16.20185174601965*vtSq[1]*f[3]+16.20185174601965*vtSq[0]*f[2])+rdv2nu[0]*(1.870828693386971*drBar[3]*f[8]+4.183300132670379*drBar[1]*f[7]+4.183300132670378*drBar[0]*f[5]+1.870828693386971*drBar[2]*f[4]+1.870828693386971*drBar[1]*f[1]+1.870828693386971*drBar[0]*f[0]); 
  out[10] += rdv2nu[0]*(0.7302967433402215*drBar[2]*f[8]+1.224744871391589*drBar[0]*f[8]+0.7302967433402215*drBar[3]*f[4]+1.075705748400954*drBar[1]*f[4]+1.224744871391589*f[0]*drBar[3]+1.075705748400954*f[1]*drBar[2])-1.0*f[10]*nu; 
  out[11] += (-3.0*f[11]*nu)-4.58257569495584*f[3]*nu+rdvSq4nu[0]*(14.23024947075771*vtSq[2]*f[10]+14.2302494707577*vtSq[3]*f[6]+14.49137674618944*vtSq[1]*f[6]+14.49137674618944*vtSq[2]*f[3]+16.20185174601965*vtSq[0]*f[3]+16.20185174601965*vtSq[1]*f[2])+rdv2nu[0]*(1.643167672515498*drBar[2]*f[8]+3.741657386773941*drBar[2]*f[7]+4.183300132670377*drBar[0]*f[7]+4.183300132670378*drBar[1]*f[5]+1.643167672515498*drBar[3]*f[4]+1.673320053068151*drBar[1]*f[4]+1.673320053068151*f[1]*drBar[2]+1.870828693386971*drBar[0]*f[1]+1.870828693386971*f[0]*drBar[1]); 

return rdvSq4nu[0]*0.5; 

} 