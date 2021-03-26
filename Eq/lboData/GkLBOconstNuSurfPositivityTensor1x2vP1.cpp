#include <GkLBOModDecl.h> 
double GkLBOconstNuSurfPositivity1x2vTensor_Vpar_P1(const double m_, const double cflL, const double cflR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[4]; 
  alphaDrSurf[0] = (2.0*wl[1]+dxvl[1])*nuSum-1.414213562373095*nuUSum[0]; 
  alphaDrSurf[1] = -1.414213562373095*nuUSum[1]; 

  double f0Quad[4]; 
  double f1Quad[4]; 
  double limQuad[4]; 
  double alphaQuad; 
  // Determine upwinding at each surface quadrature node.
  alphaQuad = 0.5*alphaDrSurf[1]-0.5*alphaDrSurf[0]; 
  if(alphaQuad > 0) {
     f0Quad[0] = 0.5*(fl[5]-1.0*(fl[3]+fl[1])+fl[0]); 
     f1Quad[0] = -0.5*(fl[7]-1.0*(fl[6]+fl[4])+fl[2]); 
     limQuad[0] = fl[0]/cflL*0.3535533905932737; 
  } else {
     f0Quad[0] = 0.5*(fr[5]-1.0*(fr[3]+fr[1])+fr[0]); 
     f1Quad[0] = 0.5*(fr[7]-1.0*(fr[6]+fr[4])+fr[2]); 
     limQuad[0] = fr[0]/cflR*0.3535533905932737; 
  } 
  alphaQuad = -0.5*(alphaDrSurf[1]+alphaDrSurf[0]); 
  if(alphaQuad > 0) {
     f0Quad[1] = -0.5*(fl[5]+fl[3]-1.0*(fl[1]+fl[0])); 
     f1Quad[1] = 0.5*(fl[7]+fl[6]-1.0*(fl[4]+fl[2])); 
     limQuad[1] = fl[0]/cflL*0.3535533905932737; 
  } else {
     f0Quad[1] = -0.5*(fr[5]+fr[3]-1.0*(fr[1]+fr[0])); 
     f1Quad[1] = -0.5*(fr[7]+fr[6]-1.0*(fr[4]+fr[2])); 
     limQuad[1] = fr[0]/cflR*0.3535533905932737; 
  } 
  alphaQuad = 0.5*alphaDrSurf[1]-0.5*alphaDrSurf[0]; 
  if(alphaQuad > 0) {
     f0Quad[2] = -0.5*(fl[5]-1.0*fl[3]+fl[1]-1.0*fl[0]); 
     f1Quad[2] = 0.5*(fl[7]-1.0*fl[6]+fl[4]-1.0*fl[2]); 
     limQuad[2] = fl[0]/cflL*0.3535533905932737; 
  } else {
     f0Quad[2] = -0.5*(fr[5]-1.0*fr[3]+fr[1]-1.0*fr[0]); 
     f1Quad[2] = -0.5*(fr[7]-1.0*fr[6]+fr[4]-1.0*fr[2]); 
     limQuad[2] = fr[0]/cflR*0.3535533905932737; 
  } 
  alphaQuad = -0.5*(alphaDrSurf[1]+alphaDrSurf[0]); 
  if(alphaQuad > 0) {
     f0Quad[3] = 0.5*(fl[5]+fl[3]+fl[1]+fl[0]); 
     f1Quad[3] = -0.5*(fl[7]+fl[6]+fl[4]+fl[2]); 
     limQuad[3] = fl[0]/cflL*0.3535533905932737; 
  } else {
     f0Quad[3] = 0.5*(fr[5]+fr[3]+fr[1]+fr[0]); 
     f1Quad[3] = 0.5*(fr[7]+fr[6]+fr[4]+fr[2]); 
     limQuad[3] = fr[0]/cflR*0.3535533905932737; 
  } 

  double fhat[8]; // (Volume) mode coefficients of fhat.
  fhat[0] = 0.5*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.5*(f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.5*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[3] = 0.5*(f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.5*(f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[5] = 0.5*(f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[6] = 0.5*(f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[7] = 0.5*(f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 

  double rCtrl[4];  // rCtrl=f1/f0 at each control node in dimensions other than vx.
  rCtrl[0] = (1.414213562373095*(3.0*fhat[7]-5.196152422706631*(fhat[6]+fhat[4])+9.0*fhat[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fhat[5]-3.0*(fhat[3]+fhat[1])+5.196152422706631*fhat[0])); 
  rCtrl[1] = -(1.414213562373095*(3.0*fhat[7]+5.196152422706631*fhat[6]-1.0*(5.196152422706631*fhat[4]+9.0*fhat[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fhat[5])+3.0*(fhat[1]-1.0*fhat[3])+5.196152422706631*fhat[0])); 
  rCtrl[2] = -(1.414213562373095*(3.0*fhat[7]+5.196152422706631*(fhat[4]-1.0*fhat[6])-9.0*fhat[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fhat[5])+3.0*(fhat[3]-1.0*fhat[1])+5.196152422706631*fhat[0])); 
  rCtrl[3] = (1.414213562373095*(3.0*fhat[7]+5.196152422706631*(fhat[6]+fhat[4])+9.0*fhat[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fhat[5]+3.0*(fhat[3]+fhat[1])+5.196152422706631*fhat[0])); 

  double fhatCtrl[4];  // fhatCtrl = anti-limited fhat evaluated at each control node on vx surface.
  fhatCtrl[0] = 0.06804138174397717*(1.732050807568877*fhat[5]-3.0*(fhat[3]+fhat[1])+5.196152422706631*fhat[0])*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = -0.06804138174397717*(1.732050807568877*fhat[5]+3.0*fhat[3]-1.0*(3.0*fhat[1]+5.196152422706631*fhat[0]))*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = -0.06804138174397717*(1.732050807568877*fhat[5]+3.0*(fhat[1]-1.0*fhat[3])-5.196152422706631*fhat[0])*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = 0.06804138174397717*(1.732050807568877*fhat[5]+3.0*(fhat[3]+fhat[1])+5.196152422706631*fhat[0])*limTheta(rCtrl[3],-1.0); 

  double fhatAL[4];  // fhatAL = mode coefficients of anti-limited f on surface.
  fhatAL[0] = 0.5*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.8660254037844386*(fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.8660254037844386*(fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 1.5*(fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 

  // Enforce limiters at surface quadrature nodes.
  double fhatALQuad[4]; 
  fhatALQuad[0] = std::max(0., std::min(0.5*((-0.5773502691896258*(1.732050807568877*fhatAL[1]-1.732050807568877*fhatAL[3]))-1.0*fhatAL[2]+fhatAL[0]), limQuad[0])); 
  fhatALQuad[1] = std::max(0., std::min(0.5*(0.5773502691896258*(1.732050807568877*fhatAL[1]-1.732050807568877*fhatAL[3])-1.0*fhatAL[2]+fhatAL[0]), limQuad[1])); 
  fhatALQuad[2] = std::max(0., std::min(0.5*((-1.0*(fhatAL[3]+fhatAL[1]))+fhatAL[2]+fhatAL[0]), limQuad[2])); 
  fhatALQuad[3] = std::max(0., std::min(0.5*(1.0*(fhatAL[3]+fhatAL[1])+fhatAL[2]+fhatAL[0]), limQuad[3])); 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // Begin surface update.
 
  double Gdiff[8]; 
  double Ghat[8]; 
  double incr2[8]; 

  if ( ((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]+0.05892556509887893*fr[5]+0.05892556509887893*fl[5]+0.1178511301977579*fr[4]-0.1178511301977579*fl[4]-0.1020620726159657*fr[3]-0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && (0.06804138174397717*fr[7]-0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]-0.05892556509887893*fr[5]-0.05892556509887893*fl[5]+0.1178511301977579*fr[4]-0.1178511301977579*fl[4]+0.1020620726159657*fr[3]+0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && (0.06804138174397717*fr[7]-0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]-0.05892556509887893*fr[5]-0.05892556509887893*fl[5]-0.1178511301977579*fr[4]+0.1178511301977579*fl[4]-0.1020620726159657*fr[3]-0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && ((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]+0.05892556509887893*fr[5]+0.05892556509887893*fl[5]-0.1178511301977579*fr[4]+0.1178511301977579*fl[4]+0.1020620726159657*fr[3]+0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) ) {

  incr2[2] = -0.1767766952966368*(2.0*nuVtSqSum[1]*fr[4]-2.0*nuVtSqSum[1]*fl[4]+2.0*nuVtSqSum[0]*fr[2]-2.0*nuVtSqSum[0]*fl[2]+((-1.732050807568877*fr[1])-1.732050807568877*fl[1])*nuVtSqSum[1]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*nuVtSqSum[0]); 
  incr2[4] = -0.1767766952966368*(2.0*nuVtSqSum[0]*fr[4]-2.0*nuVtSqSum[0]*fl[4]+2.0*nuVtSqSum[1]*fr[2]-2.0*nuVtSqSum[1]*fl[2]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0]*fr[1]-1.732050807568877*nuVtSqSum[0]*fl[1]); 
  incr2[6] = -0.1767766952966368*(2.0*nuVtSqSum[1]*fr[7]-2.0*nuVtSqSum[1]*fl[7]+2.0*nuVtSqSum[0]*fr[6]-2.0*nuVtSqSum[0]*fl[6]-1.732050807568877*nuVtSqSum[1]*fr[5]-1.732050807568877*nuVtSqSum[1]*fl[5]-1.732050807568877*nuVtSqSum[0]*fr[3]-1.732050807568877*nuVtSqSum[0]*fl[3]); 
  incr2[7] = -0.1767766952966368*(2.0*nuVtSqSum[0]*fr[7]-2.0*nuVtSqSum[0]*fl[7]+2.0*nuVtSqSum[1]*fr[6]-2.0*nuVtSqSum[1]*fl[6]-1.732050807568877*nuVtSqSum[0]*fr[5]-1.732050807568877*nuVtSqSum[0]*fl[5]-1.732050807568877*nuVtSqSum[1]*fr[3]-1.732050807568877*nuVtSqSum[1]*fl[3]); 


  Gdiff[0] = -0.0883883476483184*(8.660254037844386*nuVtSqSum[1]*fr[4]+8.660254037844386*nuVtSqSum[1]*fl[4]+8.660254037844386*nuVtSqSum[0]*fr[2]+8.660254037844386*nuVtSqSum[0]*fl[2]+(9.0*fl[1]-9.0*fr[1])*nuVtSqSum[1]+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[0]); 
  Gdiff[1] = -0.0883883476483184*(8.660254037844386*nuVtSqSum[0]*fr[4]+8.660254037844386*nuVtSqSum[0]*fl[4]+8.660254037844386*nuVtSqSum[1]*fr[2]+8.660254037844386*nuVtSqSum[1]*fl[2]+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[1]-9.0*nuVtSqSum[0]*fr[1]+9.0*nuVtSqSum[0]*fl[1]); 
  Gdiff[3] = -0.0883883476483184*(8.660254037844386*nuVtSqSum[1]*fr[7]+8.660254037844386*nuVtSqSum[1]*fl[7]+8.660254037844386*nuVtSqSum[0]*fr[6]+8.660254037844386*nuVtSqSum[0]*fl[6]-9.0*nuVtSqSum[1]*fr[5]+9.0*nuVtSqSum[1]*fl[5]-9.0*nuVtSqSum[0]*fr[3]+9.0*nuVtSqSum[0]*fl[3]); 
  Gdiff[5] = -0.0883883476483184*(8.660254037844386*nuVtSqSum[0]*fr[7]+8.660254037844386*nuVtSqSum[0]*fl[7]+8.660254037844386*nuVtSqSum[1]*fr[6]+8.660254037844386*nuVtSqSum[1]*fl[6]-9.0*nuVtSqSum[0]*fr[5]+9.0*nuVtSqSum[0]*fl[5]-9.0*nuVtSqSum[1]*fr[3]+9.0*nuVtSqSum[1]*fl[3]); 

  Ghat[0] = Gdiff[0]*rdv2L+0.7071067811865475*alphaDrSurf[1]*fhatAL[1]+0.7071067811865475*alphaDrSurf[0]*fhatAL[0]; 
  Ghat[1] = Gdiff[1]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fhatAL[1]+0.7071067811865475*fhatAL[0]*alphaDrSurf[1]; 
  Ghat[3] = Gdiff[3]*rdv2L+0.7071067811865475*alphaDrSurf[1]*fhatAL[3]+0.7071067811865475*alphaDrSurf[0]*fhatAL[2]; 
  Ghat[5] = Gdiff[5]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fhatAL[3]+0.7071067811865475*alphaDrSurf[1]*fhatAL[2]; 
  } else {

    double xBar[4];
    xBar[0] = ((-0.08930431353897002*fr[7])-0.08930431353897002*fl[7]+0.1546796083845572*fr[6]+0.1546796083845572*fl[6]+0.110485434560398*fr[5]-0.110485434560398*fl[5]+0.1546796083845572*fr[4]+0.1546796083845572*fl[4]-0.1913663861549357*fr[3]+0.1913663861549357*fl[3]-0.26791294061691*fr[2]-0.26791294061691*fl[2]-0.1913663861549357*fr[1]+0.1913663861549357*fl[1]+0.3314563036811939*fr[0]-0.3314563036811939*fl[0])/(0.5*(0.1020620726159657*fr[7]-0.1020620726159657*fl[7]-0.1767766952966368*fr[6]+0.1767766952966368*fl[6]-0.1767766952966368*fr[4]+0.1767766952966368*fl[4]+0.3061862178478971*fr[2]-0.3061862178478971*fl[2])+3.0*((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]+0.05892556509887893*fr[5]+0.05892556509887893*fl[5]+0.1178511301977579*fr[4]-0.1178511301977579*fl[4]-0.1020620726159657*fr[3]-0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0])); 
    xBar[1] = (0.08930431353897002*fr[7]+0.08930431353897002*fl[7]-0.1546796083845572*fr[6]-0.1546796083845572*fl[6]-0.110485434560398*fr[5]+0.110485434560398*fl[5]+0.1546796083845572*fr[4]+0.1546796083845572*fl[4]+0.1913663861549357*fr[3]-0.1913663861549357*fl[3]-0.26791294061691*fr[2]-0.26791294061691*fl[2]-0.1913663861549357*fr[1]+0.1913663861549357*fl[1]+0.3314563036811939*fr[0]-0.3314563036811939*fl[0])/(3.0*(0.06804138174397717*fr[7]-0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]-0.05892556509887893*fr[5]-0.05892556509887893*fl[5]+0.1178511301977579*fr[4]-0.1178511301977579*fl[4]+0.1020620726159657*fr[3]+0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0])+0.5*((-0.1020620726159657*fr[7])+0.1020620726159657*fl[7]+0.1767766952966368*fr[6]-0.1767766952966368*fl[6]-0.1767766952966368*fr[4]+0.1767766952966368*fl[4]+0.3061862178478971*fr[2]-0.3061862178478971*fl[2])); 
    xBar[2] = (0.08930431353897002*fr[7]+0.08930431353897002*fl[7]+0.1546796083845572*fr[6]+0.1546796083845572*fl[6]-0.110485434560398*fr[5]+0.110485434560398*fl[5]-0.1546796083845572*fr[4]-0.1546796083845572*fl[4]-0.1913663861549357*fr[3]+0.1913663861549357*fl[3]-0.26791294061691*fr[2]-0.26791294061691*fl[2]+0.1913663861549357*fr[1]-0.1913663861549357*fl[1]+0.3314563036811939*fr[0]-0.3314563036811939*fl[0])/(3.0*(0.06804138174397717*fr[7]-0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]-0.05892556509887893*fr[5]-0.05892556509887893*fl[5]-0.1178511301977579*fr[4]+0.1178511301977579*fl[4]-0.1020620726159657*fr[3]-0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0])+0.5*((-0.1020620726159657*fr[7])+0.1020620726159657*fl[7]-0.1767766952966368*fr[6]+0.1767766952966368*fl[6]+0.1767766952966368*fr[4]-0.1767766952966368*fl[4]+0.3061862178478971*fr[2]-0.3061862178478971*fl[2])); 
    xBar[3] = ((-0.08930431353897002*fr[7])-0.08930431353897002*fl[7]-0.1546796083845572*fr[6]-0.1546796083845572*fl[6]+0.110485434560398*fr[5]-0.110485434560398*fl[5]-0.1546796083845572*fr[4]-0.1546796083845572*fl[4]+0.1913663861549357*fr[3]-0.1913663861549357*fl[3]-0.26791294061691*fr[2]-0.26791294061691*fl[2]+0.1913663861549357*fr[1]-0.1913663861549357*fl[1]+0.3314563036811939*fr[0]-0.3314563036811939*fl[0])/(0.5*(0.1020620726159657*fr[7]-0.1020620726159657*fl[7]+0.1767766952966368*fr[6]-0.1767766952966368*fl[6]+0.1767766952966368*fr[4]-0.1767766952966368*fl[4]+0.3061862178478971*fr[2]-0.3061862178478971*fl[2])+3.0*((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]+0.05892556509887893*fr[5]+0.05892556509887893*fl[5]-0.1178511301977579*fr[4]+0.1178511301977579*fl[4]+0.1020620726159657*fr[3]+0.1020620726159657*fl[3]-0.2041241452319315*fr[2]+0.2041241452319315*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0])); 

    double xBarSq[4];
    xBarSq[0] = xBar[0]*xBar[0]; 
    xBarSq[1] = xBar[1]*xBar[1]; 
    xBarSq[2] = xBar[2]*xBar[2]; 
    xBarSq[3] = xBar[3]*xBar[3]; 

    double g1[4];
    g1[0] = (3.0*xBar[0])/(1.0-1.0*xBarSq[0])-(1.0*xBar[0]*xBarSq[0])/(1.0-1.0*xBarSq[0]); 
    g1[1] = (3.0*xBar[1])/(1.0-1.0*xBarSq[1])-(1.0*xBar[1]*xBarSq[1])/(1.0-1.0*xBarSq[1]); 
    g1[2] = (3.0*xBar[2])/(1.0-1.0*xBarSq[2])-(1.0*xBar[2]*xBarSq[2])/(1.0-1.0*xBarSq[2]); 
    g1[3] = (3.0*xBar[3])/(1.0-1.0*xBarSq[3])-(1.0*xBar[3]*xBarSq[3])/(1.0-1.0*xBarSq[3]); 

    double gBound[4];
    double gBoundP[4];

    if (std::abs(g1[0]) > 1.0e-15) {
      double g1Sq = g1[0]*g1[0];
      gBound[0] = (-(0.05103103630798288*g1[0]*fr[7])/std::sinh(g1[0]))+(0.05103103630798288*g1[0]*fl[7])/std::sinh(g1[0])+(0.08838834764831843*g1[0]*fr[6])/std::sinh(g1[0])-(0.08838834764831843*g1[0]*fl[6])/std::sinh(g1[0])+(0.05892556509887893*g1[0]*fr[5])/std::sinh(g1[0])+(0.05892556509887893*g1[0]*fl[5])/std::sinh(g1[0])+(0.08838834764831843*g1[0]*fr[4])/std::sinh(g1[0])-(0.08838834764831843*g1[0]*fl[4])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fr[3])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fl[3])/std::sinh(g1[0])-(0.1530931089239486*g1[0]*fr[2])/std::sinh(g1[0])+(0.1530931089239486*g1[0]*fl[2])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fr[1])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fl[1])/std::sinh(g1[0])+(0.1767766952966368*fr[0]*g1[0])/std::sinh(g1[0])+(0.1767766952966368*fl[0]*g1[0])/std::sinh(g1[0]); 
      gBoundP[0] = (-(0.05103103630798288*g1Sq*fr[7])/std::sinh(g1[0]))+(0.05103103630798288*g1Sq*fl[7])/std::sinh(g1[0])+(0.08838834764831843*g1Sq*fr[6])/std::sinh(g1[0])-(0.08838834764831843*g1Sq*fl[6])/std::sinh(g1[0])+(0.05892556509887893*g1Sq*fr[5])/std::sinh(g1[0])+(0.05892556509887893*g1Sq*fl[5])/std::sinh(g1[0])+(0.08838834764831843*g1Sq*fr[4])/std::sinh(g1[0])-(0.08838834764831843*g1Sq*fl[4])/std::sinh(g1[0])-(0.1020620726159657*g1Sq*fr[3])/std::sinh(g1[0])-(0.1020620726159657*g1Sq*fl[3])/std::sinh(g1[0])-(0.1530931089239486*g1Sq*fr[2])/std::sinh(g1[0])+(0.1530931089239486*g1Sq*fl[2])/std::sinh(g1[0])-(0.1020620726159657*g1Sq*fr[1])/std::sinh(g1[0])-(0.1020620726159657*g1Sq*fl[1])/std::sinh(g1[0])+(0.1767766952966368*fr[0]*g1Sq)/std::sinh(g1[0])+(0.1767766952966368*fl[0]*g1Sq)/std::sinh(g1[0]); 
    } else {
      gBound[0] = (-0.05103103630798288*fr[7])+0.05103103630798288*fl[7]+0.08838834764831843*fr[6]-0.08838834764831843*fl[6]+0.05892556509887893*fr[5]+0.05892556509887893*fl[5]+0.08838834764831843*fr[4]-0.08838834764831843*fl[4]-0.1020620726159657*fr[3]-0.1020620726159657*fl[3]-0.1530931089239486*fr[2]+0.1530931089239486*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]; 
    };

    if (std::abs(g1[1]) > 1.0e-15) {
      double g1Sq = g1[1]*g1[1];
      gBound[1] = (0.05103103630798288*g1[1]*fr[7])/std::sinh(g1[1])-(0.05103103630798288*g1[1]*fl[7])/std::sinh(g1[1])-(0.08838834764831843*g1[1]*fr[6])/std::sinh(g1[1])+(0.08838834764831843*g1[1]*fl[6])/std::sinh(g1[1])-(0.05892556509887893*g1[1]*fr[5])/std::sinh(g1[1])-(0.05892556509887893*g1[1]*fl[5])/std::sinh(g1[1])+(0.08838834764831843*g1[1]*fr[4])/std::sinh(g1[1])-(0.08838834764831843*g1[1]*fl[4])/std::sinh(g1[1])+(0.1020620726159657*g1[1]*fr[3])/std::sinh(g1[1])+(0.1020620726159657*g1[1]*fl[3])/std::sinh(g1[1])-(0.1530931089239486*g1[1]*fr[2])/std::sinh(g1[1])+(0.1530931089239486*g1[1]*fl[2])/std::sinh(g1[1])-(0.1020620726159657*fr[1]*g1[1])/std::sinh(g1[1])-(0.1020620726159657*fl[1]*g1[1])/std::sinh(g1[1])+(0.1767766952966368*fr[0]*g1[1])/std::sinh(g1[1])+(0.1767766952966368*fl[0]*g1[1])/std::sinh(g1[1]); 
      gBoundP[1] = (0.05103103630798288*g1Sq*fr[7])/std::sinh(g1[1])-(0.05103103630798288*g1Sq*fl[7])/std::sinh(g1[1])-(0.08838834764831843*g1Sq*fr[6])/std::sinh(g1[1])+(0.08838834764831843*g1Sq*fl[6])/std::sinh(g1[1])-(0.05892556509887893*g1Sq*fr[5])/std::sinh(g1[1])-(0.05892556509887893*g1Sq*fl[5])/std::sinh(g1[1])+(0.08838834764831843*g1Sq*fr[4])/std::sinh(g1[1])-(0.08838834764831843*g1Sq*fl[4])/std::sinh(g1[1])+(0.1020620726159657*g1Sq*fr[3])/std::sinh(g1[1])+(0.1020620726159657*g1Sq*fl[3])/std::sinh(g1[1])-(0.1530931089239486*g1Sq*fr[2])/std::sinh(g1[1])+(0.1530931089239486*g1Sq*fl[2])/std::sinh(g1[1])-(0.1020620726159657*fr[1]*g1Sq)/std::sinh(g1[1])-(0.1020620726159657*fl[1]*g1Sq)/std::sinh(g1[1])+(0.1767766952966368*fr[0]*g1Sq)/std::sinh(g1[1])+(0.1767766952966368*fl[0]*g1Sq)/std::sinh(g1[1]); 
    } else {
      gBound[1] = 0.05103103630798288*fr[7]-0.05103103630798288*fl[7]-0.08838834764831843*fr[6]+0.08838834764831843*fl[6]-0.05892556509887893*fr[5]-0.05892556509887893*fl[5]+0.08838834764831843*fr[4]-0.08838834764831843*fl[4]+0.1020620726159657*fr[3]+0.1020620726159657*fl[3]-0.1530931089239486*fr[2]+0.1530931089239486*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]; 
    };

    if (std::abs(g1[2]) > 1.0e-15) {
      double g1Sq = g1[2]*g1[2];
      gBound[2] = (0.05103103630798288*g1[2]*fr[7])/std::sinh(g1[2])-(0.05103103630798288*g1[2]*fl[7])/std::sinh(g1[2])+(0.08838834764831843*g1[2]*fr[6])/std::sinh(g1[2])-(0.08838834764831843*g1[2]*fl[6])/std::sinh(g1[2])-(0.05892556509887893*g1[2]*fr[5])/std::sinh(g1[2])-(0.05892556509887893*g1[2]*fl[5])/std::sinh(g1[2])-(0.08838834764831843*g1[2]*fr[4])/std::sinh(g1[2])+(0.08838834764831843*g1[2]*fl[4])/std::sinh(g1[2])-(0.1020620726159657*g1[2]*fr[3])/std::sinh(g1[2])-(0.1020620726159657*g1[2]*fl[3])/std::sinh(g1[2])-(0.1530931089239486*fr[2]*g1[2])/std::sinh(g1[2])+(0.1530931089239486*fl[2]*g1[2])/std::sinh(g1[2])+(0.1020620726159657*fr[1]*g1[2])/std::sinh(g1[2])+(0.1020620726159657*fl[1]*g1[2])/std::sinh(g1[2])+(0.1767766952966368*fr[0]*g1[2])/std::sinh(g1[2])+(0.1767766952966368*fl[0]*g1[2])/std::sinh(g1[2]); 
      gBoundP[2] = (0.05103103630798288*g1Sq*fr[7])/std::sinh(g1[2])-(0.05103103630798288*g1Sq*fl[7])/std::sinh(g1[2])+(0.08838834764831843*g1Sq*fr[6])/std::sinh(g1[2])-(0.08838834764831843*g1Sq*fl[6])/std::sinh(g1[2])-(0.05892556509887893*g1Sq*fr[5])/std::sinh(g1[2])-(0.05892556509887893*g1Sq*fl[5])/std::sinh(g1[2])-(0.08838834764831843*g1Sq*fr[4])/std::sinh(g1[2])+(0.08838834764831843*g1Sq*fl[4])/std::sinh(g1[2])-(0.1020620726159657*g1Sq*fr[3])/std::sinh(g1[2])-(0.1020620726159657*g1Sq*fl[3])/std::sinh(g1[2])-(0.1530931089239486*fr[2]*g1Sq)/std::sinh(g1[2])+(0.1530931089239486*fl[2]*g1Sq)/std::sinh(g1[2])+(0.1020620726159657*fr[1]*g1Sq)/std::sinh(g1[2])+(0.1020620726159657*fl[1]*g1Sq)/std::sinh(g1[2])+(0.1767766952966368*fr[0]*g1Sq)/std::sinh(g1[2])+(0.1767766952966368*fl[0]*g1Sq)/std::sinh(g1[2]); 
    } else {
      gBound[2] = 0.05103103630798288*fr[7]-0.05103103630798288*fl[7]+0.08838834764831843*fr[6]-0.08838834764831843*fl[6]-0.05892556509887893*fr[5]-0.05892556509887893*fl[5]-0.08838834764831843*fr[4]+0.08838834764831843*fl[4]-0.1020620726159657*fr[3]-0.1020620726159657*fl[3]-0.1530931089239486*fr[2]+0.1530931089239486*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]; 
    };

    if (std::abs(g1[3]) > 1.0e-15) {
      double g1Sq = g1[3]*g1[3];
      gBound[3] = (-(0.05103103630798288*g1[3]*fr[7])/std::sinh(g1[3]))+(0.05103103630798288*g1[3]*fl[7])/std::sinh(g1[3])-(0.08838834764831843*g1[3]*fr[6])/std::sinh(g1[3])+(0.08838834764831843*g1[3]*fl[6])/std::sinh(g1[3])+(0.05892556509887893*g1[3]*fr[5])/std::sinh(g1[3])+(0.05892556509887893*g1[3]*fl[5])/std::sinh(g1[3])-(0.08838834764831843*g1[3]*fr[4])/std::sinh(g1[3])+(0.08838834764831843*g1[3]*fl[4])/std::sinh(g1[3])+(0.1020620726159657*fr[3]*g1[3])/std::sinh(g1[3])+(0.1020620726159657*fl[3]*g1[3])/std::sinh(g1[3])-(0.1530931089239486*fr[2]*g1[3])/std::sinh(g1[3])+(0.1530931089239486*fl[2]*g1[3])/std::sinh(g1[3])+(0.1020620726159657*fr[1]*g1[3])/std::sinh(g1[3])+(0.1020620726159657*fl[1]*g1[3])/std::sinh(g1[3])+(0.1767766952966368*fr[0]*g1[3])/std::sinh(g1[3])+(0.1767766952966368*fl[0]*g1[3])/std::sinh(g1[3]); 
      gBoundP[3] = (-(0.05103103630798288*g1Sq*fr[7])/std::sinh(g1[3]))+(0.05103103630798288*g1Sq*fl[7])/std::sinh(g1[3])-(0.08838834764831843*g1Sq*fr[6])/std::sinh(g1[3])+(0.08838834764831843*g1Sq*fl[6])/std::sinh(g1[3])+(0.05892556509887893*g1Sq*fr[5])/std::sinh(g1[3])+(0.05892556509887893*g1Sq*fl[5])/std::sinh(g1[3])-(0.08838834764831843*g1Sq*fr[4])/std::sinh(g1[3])+(0.08838834764831843*g1Sq*fl[4])/std::sinh(g1[3])+(0.1020620726159657*fr[3]*g1Sq)/std::sinh(g1[3])+(0.1020620726159657*fl[3]*g1Sq)/std::sinh(g1[3])-(0.1530931089239486*fr[2]*g1Sq)/std::sinh(g1[3])+(0.1530931089239486*fl[2]*g1Sq)/std::sinh(g1[3])+(0.1020620726159657*fr[1]*g1Sq)/std::sinh(g1[3])+(0.1020620726159657*fl[1]*g1Sq)/std::sinh(g1[3])+(0.1767766952966368*fr[0]*g1Sq)/std::sinh(g1[3])+(0.1767766952966368*fl[0]*g1Sq)/std::sinh(g1[3]); 
    } else {
      gBound[3] = (-0.05103103630798288*fr[7])+0.05103103630798288*fl[7]-0.08838834764831843*fr[6]+0.08838834764831843*fl[6]+0.05892556509887893*fr[5]+0.05892556509887893*fl[5]-0.08838834764831843*fr[4]+0.08838834764831843*fl[4]+0.1020620726159657*fr[3]+0.1020620726159657*fl[3]-0.1530931089239486*fr[2]+0.1530931089239486*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]; 
    };


  incr2[2] = 0.25*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*gBound[3]+(3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*gBound[2]+((-3.0*gBound[1])-3.0*gBound[0])*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0]*gBound[1]+1.732050807568877*gBound[0]*nuVtSqSum[0]); 
  incr2[4] = 0.25*((1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*gBound[3]+(1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*gBound[2]+(1.732050807568877*gBound[1]+1.732050807568877*gBound[0])*nuVtSqSum[1]-3.0*nuVtSqSum[0]*gBound[1]-3.0*gBound[0]*nuVtSqSum[0]); 
  incr2[6] = 0.25*((5.196152422706631*nuVtSqSum[1]+3.0*nuVtSqSum[0])*gBound[3]+((-5.196152422706631*nuVtSqSum[1])-3.0*nuVtSqSum[0])*gBound[2]+(5.196152422706631*gBound[0]-5.196152422706631*gBound[1])*nuVtSqSum[1]+3.0*nuVtSqSum[0]*gBound[1]-3.0*gBound[0]*nuVtSqSum[0]); 
  incr2[7] = 0.25*((3.0*nuVtSqSum[1]+5.196152422706631*nuVtSqSum[0])*gBound[3]+((-3.0*nuVtSqSum[1])-5.196152422706631*nuVtSqSum[0])*gBound[2]+(3.0*gBound[1]-3.0*gBound[0])*nuVtSqSum[1]-5.196152422706631*nuVtSqSum[0]*gBound[1]+5.196152422706631*gBound[0]*nuVtSqSum[0]); 


  Gdiff[0] = 0.5*((1.732050807568877*nuVtSqSum[1]+nuVtSqSum[0])*gBoundP[3]+(1.732050807568877*nuVtSqSum[1]+nuVtSqSum[0])*gBoundP[2]+((-1.732050807568877*gBoundP[1])-1.732050807568877*gBoundP[0])*nuVtSqSum[1]+nuVtSqSum[0]*gBoundP[1]+gBoundP[0]*nuVtSqSum[0]); 
  Gdiff[1] = 0.5*((nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*gBoundP[3]+(nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*gBoundP[2]+(gBoundP[1]+gBoundP[0])*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0]*gBoundP[1]-1.732050807568877*gBoundP[0]*nuVtSqSum[0]); 
  Gdiff[3] = 0.5*((3.0*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0])*gBoundP[3]+((-3.0*nuVtSqSum[1])-1.732050807568877*nuVtSqSum[0])*gBoundP[2]+(3.0*gBoundP[0]-3.0*gBoundP[1])*nuVtSqSum[1]+1.732050807568877*nuVtSqSum[0]*gBoundP[1]-1.732050807568877*gBoundP[0]*nuVtSqSum[0]); 
  Gdiff[5] = 0.5*((1.732050807568877*nuVtSqSum[1]+3.0*nuVtSqSum[0])*gBoundP[3]+((-1.732050807568877*nuVtSqSum[1])-3.0*nuVtSqSum[0])*gBoundP[2]+(1.732050807568877*gBoundP[1]-1.732050807568877*gBoundP[0])*nuVtSqSum[1]-3.0*nuVtSqSum[0]*gBoundP[1]+3.0*gBoundP[0]*nuVtSqSum[0]); 

  Ghat[0] = Gdiff[0]*rdv2L+0.7071067811865475*alphaDrSurf[1]*fhatAL[1]+0.7071067811865475*alphaDrSurf[0]*fhatAL[0]; 
  Ghat[1] = Gdiff[1]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fhatAL[1]+0.7071067811865475*fhatAL[0]*alphaDrSurf[1]; 
  Ghat[3] = Gdiff[3]*rdv2L+0.7071067811865475*alphaDrSurf[1]*fhatAL[3]+0.7071067811865475*alphaDrSurf[0]*fhatAL[2]; 
  Ghat[5] = Gdiff[5]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fhatAL[3]+0.7071067811865475*alphaDrSurf[1]*fhatAL[2]; 
  };

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[1]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[3]; 
  incr1[7] = 0.8660254037844386*Ghat[5]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 

  return std::abs(wl[1]-(0.7071067811865475*nuUSum[0])/nuSum); 
} 
double GkLBOconstNuSurfPositivity1x2vTensor_Mu_P1(const double m_, const double cflL, const double cflR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[4]; 
  alphaDrSurf[0] = (4.0*wl[2]+2.0*dxvl[2])*nuSum; 

  double f0Quad[4]; 
  double f1Quad[4]; 
  double limQuad[4]; 
  double alphaQuad; 
  // Determine upwinding at each surface quadrature node.
  alphaQuad = -0.5*alphaDrSurf[0]; 
  if(alphaQuad > 0) {
     f0Quad[0] = 0.5*(fl[4]-1.0*(fl[2]+fl[1])+fl[0]); 
     f1Quad[0] = -0.5*(fl[7]-1.0*(fl[6]+fl[5])+fl[3]); 
     limQuad[0] = fl[0]/cflL*0.3535533905932737; 
  } else {
     f0Quad[0] = 0.5*(fr[4]-1.0*(fr[2]+fr[1])+fr[0]); 
     f1Quad[0] = 0.5*(fr[7]-1.0*(fr[6]+fr[5])+fr[3]); 
     limQuad[0] = fr[0]/cflR*0.3535533905932737; 
  } 
  alphaQuad = -0.5*alphaDrSurf[0]; 
  if(alphaQuad > 0) {
     f0Quad[1] = -0.5*(fl[4]+fl[2]-1.0*(fl[1]+fl[0])); 
     f1Quad[1] = 0.5*(fl[7]+fl[6]-1.0*(fl[5]+fl[3])); 
     limQuad[1] = fl[0]/cflL*0.3535533905932737; 
  } else {
     f0Quad[1] = -0.5*(fr[4]+fr[2]-1.0*(fr[1]+fr[0])); 
     f1Quad[1] = -0.5*(fr[7]+fr[6]-1.0*(fr[5]+fr[3])); 
     limQuad[1] = fr[0]/cflR*0.3535533905932737; 
  } 
  alphaQuad = -0.5*alphaDrSurf[0]; 
  if(alphaQuad > 0) {
     f0Quad[2] = -0.5*(fl[4]-1.0*fl[2]+fl[1]-1.0*fl[0]); 
     f1Quad[2] = 0.5*(fl[7]-1.0*fl[6]+fl[5]-1.0*fl[3]); 
     limQuad[2] = fl[0]/cflL*0.3535533905932737; 
  } else {
     f0Quad[2] = -0.5*(fr[4]-1.0*fr[2]+fr[1]-1.0*fr[0]); 
     f1Quad[2] = -0.5*(fr[7]-1.0*fr[6]+fr[5]-1.0*fr[3]); 
     limQuad[2] = fr[0]/cflR*0.3535533905932737; 
  } 
  alphaQuad = -0.5*alphaDrSurf[0]; 
  if(alphaQuad > 0) {
     f0Quad[3] = 0.5*(fl[4]+fl[2]+fl[1]+fl[0]); 
     f1Quad[3] = -0.5*(fl[7]+fl[6]+fl[5]+fl[3]); 
     limQuad[3] = fl[0]/cflL*0.3535533905932737; 
  } else {
     f0Quad[3] = 0.5*(fr[4]+fr[2]+fr[1]+fr[0]); 
     f1Quad[3] = 0.5*(fr[7]+fr[6]+fr[5]+fr[3]); 
     limQuad[3] = fr[0]/cflR*0.3535533905932737; 
  } 

  double fhat[8]; // (Volume) mode coefficients of fhat.
  fhat[0] = 0.5*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.5*(f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.5*(f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[3] = 0.5*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[4] = 0.5*(f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[5] = 0.5*(f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[6] = 0.5*(f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[7] = 0.5*(f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 

  double rCtrl[4];  // rCtrl=f1/f0 at each control node in dimensions other than vy.
  rCtrl[0] = (1.414213562373095*(3.0*fhat[7]-5.196152422706631*(fhat[6]+fhat[5])+9.0*fhat[3]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fhat[4]-3.0*(fhat[2]+fhat[1])+5.196152422706631*fhat[0])); 
  rCtrl[1] = -(1.414213562373095*(3.0*fhat[7]+5.196152422706631*fhat[6]-1.0*(5.196152422706631*fhat[5]+9.0*fhat[3])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fhat[4])+3.0*(fhat[1]-1.0*fhat[2])+5.196152422706631*fhat[0])); 
  rCtrl[2] = -(1.414213562373095*(3.0*fhat[7]+5.196152422706631*(fhat[5]-1.0*fhat[6])-9.0*fhat[3]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fhat[4])+3.0*(fhat[2]-1.0*fhat[1])+5.196152422706631*fhat[0])); 
  rCtrl[3] = (1.414213562373095*(3.0*fhat[7]+5.196152422706631*(fhat[6]+fhat[5])+9.0*fhat[3]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fhat[4]+3.0*(fhat[2]+fhat[1])+5.196152422706631*fhat[0])); 

  double fhatCtrl[4];  // fhatCtrl = anti-limited fhat evaluated at each control node on vy surface.
  fhatCtrl[0] = 0.06804138174397717*(1.732050807568877*fhat[4]-3.0*(fhat[2]+fhat[1])+5.196152422706631*fhat[0])*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = -0.06804138174397717*(1.732050807568877*fhat[4]+3.0*fhat[2]-1.0*(3.0*fhat[1]+5.196152422706631*fhat[0]))*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = -0.06804138174397717*(1.732050807568877*fhat[4]+3.0*(fhat[1]-1.0*fhat[2])-5.196152422706631*fhat[0])*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = 0.06804138174397717*(1.732050807568877*fhat[4]+3.0*(fhat[2]+fhat[1])+5.196152422706631*fhat[0])*limTheta(rCtrl[3],-1.0); 

  double fhatAL[4];  // fhatAL = mode coefficients of anti-limited f on surface.
  fhatAL[0] = 0.5*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.8660254037844386*(fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.8660254037844386*(fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 1.5*(fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 

  // Enforce limiters at surface quadrature nodes.
  double fhatALQuad[4]; 
  fhatALQuad[0] = std::max(0., std::min(0.5*((-0.5773502691896258*(1.732050807568877*fhatAL[1]-1.732050807568877*fhatAL[3]))-1.0*fhatAL[2]+fhatAL[0]), limQuad[0])); 
  fhatALQuad[1] = std::max(0., std::min(0.5*(0.5773502691896258*(1.732050807568877*fhatAL[1]-1.732050807568877*fhatAL[3])-1.0*fhatAL[2]+fhatAL[0]), limQuad[1])); 
  fhatALQuad[2] = std::max(0., std::min(0.5*((-1.0*(fhatAL[3]+fhatAL[1]))+fhatAL[2]+fhatAL[0]), limQuad[2])); 
  fhatALQuad[3] = std::max(0., std::min(0.5*(1.0*(fhatAL[3]+fhatAL[1])+fhatAL[2]+fhatAL[0]), limQuad[3])); 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // Begin surface update.
 
  double diffFac[2]; 
  diffFac[0] = nuVtSqSum[1]*(1.414213562373095*BmagInv[1]*wl[2]*m_+0.7071067811865475*BmagInv[1]*dxvl[2]*m_)+nuVtSqSum[0]*(1.414213562373095*BmagInv[0]*wl[2]*m_+0.7071067811865475*BmagInv[0]*dxvl[2]*m_); 
  diffFac[1] = nuVtSqSum[0]*(1.414213562373095*BmagInv[1]*wl[2]*m_+0.7071067811865475*BmagInv[1]*dxvl[2]*m_)+nuVtSqSum[1]*(1.414213562373095*BmagInv[0]*wl[2]*m_+0.7071067811865475*BmagInv[0]*dxvl[2]*m_); 

  double Gdiff[8]; 
  double Ghat[8]; 
  double incr2[8]; 

  if ( ((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]+0.1178511301977579*fr[5]-0.1178511301977579*fl[5]+0.05892556509887893*fr[4]+0.05892556509887893*fl[4]-0.2041241452319315*fr[3]+0.2041241452319315*fl[3]-0.1020620726159657*fr[2]-0.1020620726159657*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && (0.06804138174397717*fr[7]-0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]+0.1178511301977579*fr[5]-0.1178511301977579*fl[5]-0.05892556509887893*fr[4]-0.05892556509887893*fl[4]-0.2041241452319315*fr[3]+0.2041241452319315*fl[3]+0.1020620726159657*fr[2]+0.1020620726159657*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && (0.06804138174397717*fr[7]-0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]-0.1178511301977579*fr[5]+0.1178511301977579*fl[5]-0.05892556509887893*fr[4]-0.05892556509887893*fl[4]-0.2041241452319315*fr[3]+0.2041241452319315*fl[3]-0.1020620726159657*fr[2]-0.1020620726159657*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) && ((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]-0.1178511301977579*fr[5]+0.1178511301977579*fl[5]+0.05892556509887893*fr[4]+0.05892556509887893*fl[4]-0.2041241452319315*fr[3]+0.2041241452319315*fl[3]+0.1020620726159657*fr[2]+0.1020620726159657*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]>=0.0) ) {

  incr2[3] = -0.1767766952966368*(2.0*diffFac[1]*fr[5]-2.0*diffFac[1]*fl[5]+2.0*diffFac[0]*fr[3]-2.0*diffFac[0]*fl[3]-1.732050807568877*diffFac[1]*fr[1]-1.732050807568877*diffFac[1]*fl[1]-1.732050807568877*diffFac[0]*fr[0]-1.732050807568877*diffFac[0]*fl[0]); 
  incr2[5] = -0.1767766952966368*(2.0*diffFac[0]*fr[5]-2.0*diffFac[0]*fl[5]+2.0*diffFac[1]*fr[3]-2.0*diffFac[1]*fl[3]-1.732050807568877*diffFac[0]*fr[1]-1.732050807568877*diffFac[0]*fl[1]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*diffFac[1]); 
  incr2[6] = -0.1767766952966368*(2.0*diffFac[1]*fr[7]-2.0*diffFac[1]*fl[7]+2.0*diffFac[0]*fr[6]-2.0*diffFac[0]*fl[6]-1.732050807568877*diffFac[1]*fr[4]-1.732050807568877*diffFac[1]*fl[4]-1.732050807568877*diffFac[0]*fr[2]-1.732050807568877*diffFac[0]*fl[2]); 
  incr2[7] = -0.1767766952966368*(2.0*diffFac[0]*fr[7]-2.0*diffFac[0]*fl[7]+2.0*diffFac[1]*fr[6]-2.0*diffFac[1]*fl[6]-1.732050807568877*diffFac[0]*fr[4]-1.732050807568877*diffFac[0]*fl[4]-1.732050807568877*diffFac[1]*fr[2]-1.732050807568877*diffFac[1]*fl[2]); 


  Gdiff[0] = -0.0883883476483184*(8.660254037844386*diffFac[1]*fr[5]+8.660254037844386*diffFac[1]*fl[5]+8.660254037844386*diffFac[0]*fr[3]+8.660254037844386*diffFac[0]*fl[3]-9.0*diffFac[1]*fr[1]+9.0*diffFac[1]*fl[1]-9.0*diffFac[0]*fr[0]+9.0*diffFac[0]*fl[0]); 
  Gdiff[1] = -0.0883883476483184*(8.660254037844386*diffFac[0]*fr[5]+8.660254037844386*diffFac[0]*fl[5]+8.660254037844386*diffFac[1]*fr[3]+8.660254037844386*diffFac[1]*fl[3]-9.0*diffFac[0]*fr[1]+9.0*diffFac[0]*fl[1]+(9.0*fl[0]-9.0*fr[0])*diffFac[1]); 
  Gdiff[2] = -0.0883883476483184*(8.660254037844386*diffFac[1]*fr[7]+8.660254037844386*diffFac[1]*fl[7]+8.660254037844386*diffFac[0]*fr[6]+8.660254037844386*diffFac[0]*fl[6]-9.0*diffFac[1]*fr[4]+9.0*diffFac[1]*fl[4]-9.0*diffFac[0]*fr[2]+9.0*diffFac[0]*fl[2]); 
  Gdiff[4] = -0.0883883476483184*(8.660254037844386*diffFac[0]*fr[7]+8.660254037844386*diffFac[0]*fl[7]+8.660254037844386*diffFac[1]*fr[6]+8.660254037844386*diffFac[1]*fl[6]-9.0*diffFac[0]*fr[4]+9.0*diffFac[0]*fl[4]-9.0*diffFac[1]*fr[2]+9.0*diffFac[1]*fl[2]); 

  Ghat[0] = Gdiff[0]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fhatAL[0]; 
  Ghat[1] = Gdiff[1]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fhatAL[1]; 
  Ghat[2] = Gdiff[2]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fhatAL[2]; 
  Ghat[4] = Gdiff[4]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fhatAL[3]; 
  } else {

    double xBar[4];
    xBar[0] = ((-0.08930431353897002*fr[7])-0.08930431353897002*fl[7]+0.1546796083845572*fr[6]+0.1546796083845572*fl[6]+0.1546796083845572*fr[5]+0.1546796083845572*fl[5]+0.110485434560398*fr[4]-0.110485434560398*fl[4]-0.26791294061691*fr[3]-0.26791294061691*fl[3]-0.1913663861549357*fr[2]+0.1913663861549357*fl[2]-0.1913663861549357*fr[1]+0.1913663861549357*fl[1]+0.3314563036811939*fr[0]-0.3314563036811939*fl[0])/(0.5*(0.1020620726159657*fr[7]-0.1020620726159657*fl[7]-0.1767766952966368*fr[6]+0.1767766952966368*fl[6]-0.1767766952966368*fr[5]+0.1767766952966368*fl[5]+0.3061862178478971*fr[3]-0.3061862178478971*fl[3])+3.0*((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]+0.1178511301977579*fr[5]-0.1178511301977579*fl[5]+0.05892556509887893*fr[4]+0.05892556509887893*fl[4]-0.2041241452319315*fr[3]+0.2041241452319315*fl[3]-0.1020620726159657*fr[2]-0.1020620726159657*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0])); 
    xBar[1] = (0.08930431353897002*fr[7]+0.08930431353897002*fl[7]-0.1546796083845572*fr[6]-0.1546796083845572*fl[6]+0.1546796083845572*fr[5]+0.1546796083845572*fl[5]-0.110485434560398*fr[4]+0.110485434560398*fl[4]-0.26791294061691*fr[3]-0.26791294061691*fl[3]+0.1913663861549357*fr[2]-0.1913663861549357*fl[2]-0.1913663861549357*fr[1]+0.1913663861549357*fl[1]+0.3314563036811939*fr[0]-0.3314563036811939*fl[0])/(3.0*(0.06804138174397717*fr[7]-0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]+0.1178511301977579*fr[5]-0.1178511301977579*fl[5]-0.05892556509887893*fr[4]-0.05892556509887893*fl[4]-0.2041241452319315*fr[3]+0.2041241452319315*fl[3]+0.1020620726159657*fr[2]+0.1020620726159657*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0])+0.5*((-0.1020620726159657*fr[7])+0.1020620726159657*fl[7]+0.1767766952966368*fr[6]-0.1767766952966368*fl[6]-0.1767766952966368*fr[5]+0.1767766952966368*fl[5]+0.3061862178478971*fr[3]-0.3061862178478971*fl[3])); 
    xBar[2] = (0.08930431353897002*fr[7]+0.08930431353897002*fl[7]+0.1546796083845572*fr[6]+0.1546796083845572*fl[6]-0.1546796083845572*fr[5]-0.1546796083845572*fl[5]-0.110485434560398*fr[4]+0.110485434560398*fl[4]-0.26791294061691*fr[3]-0.26791294061691*fl[3]-0.1913663861549357*fr[2]+0.1913663861549357*fl[2]+0.1913663861549357*fr[1]-0.1913663861549357*fl[1]+0.3314563036811939*fr[0]-0.3314563036811939*fl[0])/(3.0*(0.06804138174397717*fr[7]-0.06804138174397717*fl[7]+0.1178511301977579*fr[6]-0.1178511301977579*fl[6]-0.1178511301977579*fr[5]+0.1178511301977579*fl[5]-0.05892556509887893*fr[4]-0.05892556509887893*fl[4]-0.2041241452319315*fr[3]+0.2041241452319315*fl[3]-0.1020620726159657*fr[2]-0.1020620726159657*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0])+0.5*((-0.1020620726159657*fr[7])+0.1020620726159657*fl[7]-0.1767766952966368*fr[6]+0.1767766952966368*fl[6]+0.1767766952966368*fr[5]-0.1767766952966368*fl[5]+0.3061862178478971*fr[3]-0.3061862178478971*fl[3])); 
    xBar[3] = ((-0.08930431353897002*fr[7])-0.08930431353897002*fl[7]-0.1546796083845572*fr[6]-0.1546796083845572*fl[6]-0.1546796083845572*fr[5]-0.1546796083845572*fl[5]+0.110485434560398*fr[4]-0.110485434560398*fl[4]-0.26791294061691*fr[3]-0.26791294061691*fl[3]+0.1913663861549357*fr[2]-0.1913663861549357*fl[2]+0.1913663861549357*fr[1]-0.1913663861549357*fl[1]+0.3314563036811939*fr[0]-0.3314563036811939*fl[0])/(0.5*(0.1020620726159657*fr[7]-0.1020620726159657*fl[7]+0.1767766952966368*fr[6]-0.1767766952966368*fl[6]+0.1767766952966368*fr[5]-0.1767766952966368*fl[5]+0.3061862178478971*fr[3]-0.3061862178478971*fl[3])+3.0*((-0.06804138174397717*fr[7])+0.06804138174397717*fl[7]-0.1178511301977579*fr[6]+0.1178511301977579*fl[6]-0.1178511301977579*fr[5]+0.1178511301977579*fl[5]+0.05892556509887893*fr[4]+0.05892556509887893*fl[4]-0.2041241452319315*fr[3]+0.2041241452319315*fl[3]+0.1020620726159657*fr[2]+0.1020620726159657*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0])); 

    double xBarSq[4];
    xBarSq[0] = xBar[0]*xBar[0]; 
    xBarSq[1] = xBar[1]*xBar[1]; 
    xBarSq[2] = xBar[2]*xBar[2]; 
    xBarSq[3] = xBar[3]*xBar[3]; 

    double g1[4];
    g1[0] = (3.0*xBar[0])/(1.0-1.0*xBarSq[0])-(1.0*xBar[0]*xBarSq[0])/(1.0-1.0*xBarSq[0]); 
    g1[1] = (3.0*xBar[1])/(1.0-1.0*xBarSq[1])-(1.0*xBar[1]*xBarSq[1])/(1.0-1.0*xBarSq[1]); 
    g1[2] = (3.0*xBar[2])/(1.0-1.0*xBarSq[2])-(1.0*xBar[2]*xBarSq[2])/(1.0-1.0*xBarSq[2]); 
    g1[3] = (3.0*xBar[3])/(1.0-1.0*xBarSq[3])-(1.0*xBar[3]*xBarSq[3])/(1.0-1.0*xBarSq[3]); 

    double gBound[4];
    double gBoundP[4];

    if (std::abs(g1[0]) > 1.0e-15) {
      double g1Sq = g1[0]*g1[0];
      gBound[0] = (-(0.05103103630798288*g1[0]*fr[7])/std::sinh(g1[0]))+(0.05103103630798288*g1[0]*fl[7])/std::sinh(g1[0])+(0.08838834764831843*g1[0]*fr[6])/std::sinh(g1[0])-(0.08838834764831843*g1[0]*fl[6])/std::sinh(g1[0])+(0.08838834764831843*g1[0]*fr[5])/std::sinh(g1[0])-(0.08838834764831843*g1[0]*fl[5])/std::sinh(g1[0])+(0.05892556509887893*g1[0]*fr[4])/std::sinh(g1[0])+(0.05892556509887893*g1[0]*fl[4])/std::sinh(g1[0])-(0.1530931089239486*g1[0]*fr[3])/std::sinh(g1[0])+(0.1530931089239486*g1[0]*fl[3])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fr[2])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fl[2])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fr[1])/std::sinh(g1[0])-(0.1020620726159657*g1[0]*fl[1])/std::sinh(g1[0])+(0.1767766952966368*fr[0]*g1[0])/std::sinh(g1[0])+(0.1767766952966368*fl[0]*g1[0])/std::sinh(g1[0]); 
      gBoundP[0] = (-(0.05103103630798288*g1Sq*fr[7])/std::sinh(g1[0]))+(0.05103103630798288*g1Sq*fl[7])/std::sinh(g1[0])+(0.08838834764831843*g1Sq*fr[6])/std::sinh(g1[0])-(0.08838834764831843*g1Sq*fl[6])/std::sinh(g1[0])+(0.08838834764831843*g1Sq*fr[5])/std::sinh(g1[0])-(0.08838834764831843*g1Sq*fl[5])/std::sinh(g1[0])+(0.05892556509887893*g1Sq*fr[4])/std::sinh(g1[0])+(0.05892556509887893*g1Sq*fl[4])/std::sinh(g1[0])-(0.1530931089239486*g1Sq*fr[3])/std::sinh(g1[0])+(0.1530931089239486*g1Sq*fl[3])/std::sinh(g1[0])-(0.1020620726159657*g1Sq*fr[2])/std::sinh(g1[0])-(0.1020620726159657*g1Sq*fl[2])/std::sinh(g1[0])-(0.1020620726159657*g1Sq*fr[1])/std::sinh(g1[0])-(0.1020620726159657*g1Sq*fl[1])/std::sinh(g1[0])+(0.1767766952966368*fr[0]*g1Sq)/std::sinh(g1[0])+(0.1767766952966368*fl[0]*g1Sq)/std::sinh(g1[0]); 
    } else {
      gBound[0] = (-0.05103103630798288*fr[7])+0.05103103630798288*fl[7]+0.08838834764831843*fr[6]-0.08838834764831843*fl[6]+0.08838834764831843*fr[5]-0.08838834764831843*fl[5]+0.05892556509887893*fr[4]+0.05892556509887893*fl[4]-0.1530931089239486*fr[3]+0.1530931089239486*fl[3]-0.1020620726159657*fr[2]-0.1020620726159657*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]; 
    };

    if (std::abs(g1[1]) > 1.0e-15) {
      double g1Sq = g1[1]*g1[1];
      gBound[1] = (0.05103103630798288*g1[1]*fr[7])/std::sinh(g1[1])-(0.05103103630798288*g1[1]*fl[7])/std::sinh(g1[1])-(0.08838834764831843*g1[1]*fr[6])/std::sinh(g1[1])+(0.08838834764831843*g1[1]*fl[6])/std::sinh(g1[1])+(0.08838834764831843*g1[1]*fr[5])/std::sinh(g1[1])-(0.08838834764831843*g1[1]*fl[5])/std::sinh(g1[1])-(0.05892556509887893*g1[1]*fr[4])/std::sinh(g1[1])-(0.05892556509887893*g1[1]*fl[4])/std::sinh(g1[1])-(0.1530931089239486*g1[1]*fr[3])/std::sinh(g1[1])+(0.1530931089239486*g1[1]*fl[3])/std::sinh(g1[1])+(0.1020620726159657*g1[1]*fr[2])/std::sinh(g1[1])+(0.1020620726159657*g1[1]*fl[2])/std::sinh(g1[1])-(0.1020620726159657*fr[1]*g1[1])/std::sinh(g1[1])-(0.1020620726159657*fl[1]*g1[1])/std::sinh(g1[1])+(0.1767766952966368*fr[0]*g1[1])/std::sinh(g1[1])+(0.1767766952966368*fl[0]*g1[1])/std::sinh(g1[1]); 
      gBoundP[1] = (0.05103103630798288*g1Sq*fr[7])/std::sinh(g1[1])-(0.05103103630798288*g1Sq*fl[7])/std::sinh(g1[1])-(0.08838834764831843*g1Sq*fr[6])/std::sinh(g1[1])+(0.08838834764831843*g1Sq*fl[6])/std::sinh(g1[1])+(0.08838834764831843*g1Sq*fr[5])/std::sinh(g1[1])-(0.08838834764831843*g1Sq*fl[5])/std::sinh(g1[1])-(0.05892556509887893*g1Sq*fr[4])/std::sinh(g1[1])-(0.05892556509887893*g1Sq*fl[4])/std::sinh(g1[1])-(0.1530931089239486*g1Sq*fr[3])/std::sinh(g1[1])+(0.1530931089239486*g1Sq*fl[3])/std::sinh(g1[1])+(0.1020620726159657*g1Sq*fr[2])/std::sinh(g1[1])+(0.1020620726159657*g1Sq*fl[2])/std::sinh(g1[1])-(0.1020620726159657*fr[1]*g1Sq)/std::sinh(g1[1])-(0.1020620726159657*fl[1]*g1Sq)/std::sinh(g1[1])+(0.1767766952966368*fr[0]*g1Sq)/std::sinh(g1[1])+(0.1767766952966368*fl[0]*g1Sq)/std::sinh(g1[1]); 
    } else {
      gBound[1] = 0.05103103630798288*fr[7]-0.05103103630798288*fl[7]-0.08838834764831843*fr[6]+0.08838834764831843*fl[6]+0.08838834764831843*fr[5]-0.08838834764831843*fl[5]-0.05892556509887893*fr[4]-0.05892556509887893*fl[4]-0.1530931089239486*fr[3]+0.1530931089239486*fl[3]+0.1020620726159657*fr[2]+0.1020620726159657*fl[2]-0.1020620726159657*fr[1]-0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]; 
    };

    if (std::abs(g1[2]) > 1.0e-15) {
      double g1Sq = g1[2]*g1[2];
      gBound[2] = (0.05103103630798288*g1[2]*fr[7])/std::sinh(g1[2])-(0.05103103630798288*g1[2]*fl[7])/std::sinh(g1[2])+(0.08838834764831843*g1[2]*fr[6])/std::sinh(g1[2])-(0.08838834764831843*g1[2]*fl[6])/std::sinh(g1[2])-(0.08838834764831843*g1[2]*fr[5])/std::sinh(g1[2])+(0.08838834764831843*g1[2]*fl[5])/std::sinh(g1[2])-(0.05892556509887893*g1[2]*fr[4])/std::sinh(g1[2])-(0.05892556509887893*g1[2]*fl[4])/std::sinh(g1[2])-(0.1530931089239486*g1[2]*fr[3])/std::sinh(g1[2])+(0.1530931089239486*g1[2]*fl[3])/std::sinh(g1[2])-(0.1020620726159657*fr[2]*g1[2])/std::sinh(g1[2])-(0.1020620726159657*fl[2]*g1[2])/std::sinh(g1[2])+(0.1020620726159657*fr[1]*g1[2])/std::sinh(g1[2])+(0.1020620726159657*fl[1]*g1[2])/std::sinh(g1[2])+(0.1767766952966368*fr[0]*g1[2])/std::sinh(g1[2])+(0.1767766952966368*fl[0]*g1[2])/std::sinh(g1[2]); 
      gBoundP[2] = (0.05103103630798288*g1Sq*fr[7])/std::sinh(g1[2])-(0.05103103630798288*g1Sq*fl[7])/std::sinh(g1[2])+(0.08838834764831843*g1Sq*fr[6])/std::sinh(g1[2])-(0.08838834764831843*g1Sq*fl[6])/std::sinh(g1[2])-(0.08838834764831843*g1Sq*fr[5])/std::sinh(g1[2])+(0.08838834764831843*g1Sq*fl[5])/std::sinh(g1[2])-(0.05892556509887893*g1Sq*fr[4])/std::sinh(g1[2])-(0.05892556509887893*g1Sq*fl[4])/std::sinh(g1[2])-(0.1530931089239486*g1Sq*fr[3])/std::sinh(g1[2])+(0.1530931089239486*g1Sq*fl[3])/std::sinh(g1[2])-(0.1020620726159657*fr[2]*g1Sq)/std::sinh(g1[2])-(0.1020620726159657*fl[2]*g1Sq)/std::sinh(g1[2])+(0.1020620726159657*fr[1]*g1Sq)/std::sinh(g1[2])+(0.1020620726159657*fl[1]*g1Sq)/std::sinh(g1[2])+(0.1767766952966368*fr[0]*g1Sq)/std::sinh(g1[2])+(0.1767766952966368*fl[0]*g1Sq)/std::sinh(g1[2]); 
    } else {
      gBound[2] = 0.05103103630798288*fr[7]-0.05103103630798288*fl[7]+0.08838834764831843*fr[6]-0.08838834764831843*fl[6]-0.08838834764831843*fr[5]+0.08838834764831843*fl[5]-0.05892556509887893*fr[4]-0.05892556509887893*fl[4]-0.1530931089239486*fr[3]+0.1530931089239486*fl[3]-0.1020620726159657*fr[2]-0.1020620726159657*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]; 
    };

    if (std::abs(g1[3]) > 1.0e-15) {
      double g1Sq = g1[3]*g1[3];
      gBound[3] = (-(0.05103103630798288*g1[3]*fr[7])/std::sinh(g1[3]))+(0.05103103630798288*g1[3]*fl[7])/std::sinh(g1[3])-(0.08838834764831843*g1[3]*fr[6])/std::sinh(g1[3])+(0.08838834764831843*g1[3]*fl[6])/std::sinh(g1[3])-(0.08838834764831843*g1[3]*fr[5])/std::sinh(g1[3])+(0.08838834764831843*g1[3]*fl[5])/std::sinh(g1[3])+(0.05892556509887893*g1[3]*fr[4])/std::sinh(g1[3])+(0.05892556509887893*g1[3]*fl[4])/std::sinh(g1[3])-(0.1530931089239486*fr[3]*g1[3])/std::sinh(g1[3])+(0.1530931089239486*fl[3]*g1[3])/std::sinh(g1[3])+(0.1020620726159657*fr[2]*g1[3])/std::sinh(g1[3])+(0.1020620726159657*fl[2]*g1[3])/std::sinh(g1[3])+(0.1020620726159657*fr[1]*g1[3])/std::sinh(g1[3])+(0.1020620726159657*fl[1]*g1[3])/std::sinh(g1[3])+(0.1767766952966368*fr[0]*g1[3])/std::sinh(g1[3])+(0.1767766952966368*fl[0]*g1[3])/std::sinh(g1[3]); 
      gBoundP[3] = (-(0.05103103630798288*g1Sq*fr[7])/std::sinh(g1[3]))+(0.05103103630798288*g1Sq*fl[7])/std::sinh(g1[3])-(0.08838834764831843*g1Sq*fr[6])/std::sinh(g1[3])+(0.08838834764831843*g1Sq*fl[6])/std::sinh(g1[3])-(0.08838834764831843*g1Sq*fr[5])/std::sinh(g1[3])+(0.08838834764831843*g1Sq*fl[5])/std::sinh(g1[3])+(0.05892556509887893*g1Sq*fr[4])/std::sinh(g1[3])+(0.05892556509887893*g1Sq*fl[4])/std::sinh(g1[3])-(0.1530931089239486*fr[3]*g1Sq)/std::sinh(g1[3])+(0.1530931089239486*fl[3]*g1Sq)/std::sinh(g1[3])+(0.1020620726159657*fr[2]*g1Sq)/std::sinh(g1[3])+(0.1020620726159657*fl[2]*g1Sq)/std::sinh(g1[3])+(0.1020620726159657*fr[1]*g1Sq)/std::sinh(g1[3])+(0.1020620726159657*fl[1]*g1Sq)/std::sinh(g1[3])+(0.1767766952966368*fr[0]*g1Sq)/std::sinh(g1[3])+(0.1767766952966368*fl[0]*g1Sq)/std::sinh(g1[3]); 
    } else {
      gBound[3] = (-0.05103103630798288*fr[7])+0.05103103630798288*fl[7]-0.08838834764831843*fr[6]+0.08838834764831843*fl[6]-0.08838834764831843*fr[5]+0.08838834764831843*fl[5]+0.05892556509887893*fr[4]+0.05892556509887893*fl[4]-0.1530931089239486*fr[3]+0.1530931089239486*fl[3]+0.1020620726159657*fr[2]+0.1020620726159657*fl[2]+0.1020620726159657*fr[1]+0.1020620726159657*fl[1]+0.1767766952966368*fr[0]+0.1767766952966368*fl[0]; 
    };


  incr2[3] = 0.25*((3.0*diffFac[1]+1.732050807568877*diffFac[0])*gBound[3]+(3.0*diffFac[1]+1.732050807568877*diffFac[0])*gBound[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*gBound[1]-3.0*gBound[0]*diffFac[1]+1.732050807568877*diffFac[0]*gBound[0]); 
  incr2[5] = 0.25*((1.732050807568877*diffFac[1]+3.0*diffFac[0])*gBound[3]+(1.732050807568877*diffFac[1]+3.0*diffFac[0])*gBound[2]+(1.732050807568877*diffFac[1]-3.0*diffFac[0])*gBound[1]+1.732050807568877*gBound[0]*diffFac[1]-3.0*diffFac[0]*gBound[0]); 
  incr2[6] = 0.25*((5.196152422706631*diffFac[1]+3.0*diffFac[0])*gBound[3]+((-5.196152422706631*diffFac[1])-3.0*diffFac[0])*gBound[2]+(3.0*diffFac[0]-5.196152422706631*diffFac[1])*gBound[1]+5.196152422706631*gBound[0]*diffFac[1]-3.0*diffFac[0]*gBound[0]); 
  incr2[7] = 0.25*((3.0*diffFac[1]+5.196152422706631*diffFac[0])*gBound[3]+((-3.0*diffFac[1])-5.196152422706631*diffFac[0])*gBound[2]+(3.0*diffFac[1]-5.196152422706631*diffFac[0])*gBound[1]-3.0*gBound[0]*diffFac[1]+5.196152422706631*diffFac[0]*gBound[0]); 


  Gdiff[0] = 0.5*((1.732050807568877*diffFac[1]+diffFac[0])*gBoundP[3]+(1.732050807568877*diffFac[1]+diffFac[0])*gBoundP[2]+(diffFac[0]-1.732050807568877*diffFac[1])*gBoundP[1]-1.732050807568877*gBoundP[0]*diffFac[1]+diffFac[0]*gBoundP[0]); 
  Gdiff[1] = 0.5*((diffFac[1]+1.732050807568877*diffFac[0])*gBoundP[3]+(diffFac[1]+1.732050807568877*diffFac[0])*gBoundP[2]+(diffFac[1]-1.732050807568877*diffFac[0])*gBoundP[1]+gBoundP[0]*diffFac[1]-1.732050807568877*diffFac[0]*gBoundP[0]); 
  Gdiff[2] = 0.5*((3.0*diffFac[1]+1.732050807568877*diffFac[0])*gBoundP[3]+((-3.0*diffFac[1])-1.732050807568877*diffFac[0])*gBoundP[2]+(1.732050807568877*diffFac[0]-3.0*diffFac[1])*gBoundP[1]+3.0*gBoundP[0]*diffFac[1]-1.732050807568877*diffFac[0]*gBoundP[0]); 
  Gdiff[4] = 0.5*((1.732050807568877*diffFac[1]+3.0*diffFac[0])*gBoundP[3]+((-1.732050807568877*diffFac[1])-3.0*diffFac[0])*gBoundP[2]+(1.732050807568877*diffFac[1]-3.0*diffFac[0])*gBoundP[1]-1.732050807568877*gBoundP[0]*diffFac[1]+3.0*diffFac[0]*gBoundP[0]); 

  Ghat[0] = Gdiff[0]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fhatAL[0]; 
  Ghat[1] = Gdiff[1]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fhatAL[1]; 
  Ghat[2] = Gdiff[2]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fhatAL[2]; 
  Ghat[4] = Gdiff[4]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fhatAL[3]; 
  };

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[1]; 
  incr1[6] = 0.8660254037844386*Ghat[2]; 
  incr1[7] = 0.8660254037844386*Ghat[4]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr1[5]*rdv2L-1.0*incr2[5]*rdvSq4L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 

  return std::abs(2.0*wl[2]); 
} 
