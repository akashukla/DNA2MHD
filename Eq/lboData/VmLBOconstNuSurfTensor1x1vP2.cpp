#include <VmLBOModDecl.h> 
double VmLBOconstNuSurf1x1vTensor_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:          Cell-center coordinates. 
  // dxv[2]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[3]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double favg[9]; 
  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 

  double fjump[9]; 
  fjump[0] = nuSum*vMuMidMax*(fl[0]-(1*fr[0])); 
  fjump[1] = nuSum*vMuMidMax*(fl[1]-(1*fr[1])); 
  fjump[2] = nuSum*vMuMidMax*(fl[2]-(-1*fr[2])); 
  fjump[3] = nuSum*vMuMidMax*(fl[3]-(-1*fr[3])); 
  fjump[4] = nuSum*vMuMidMax*(fl[4]-(1*fr[4])); 
  fjump[5] = nuSum*vMuMidMax*(fl[5]-(1*fr[5])); 
  fjump[6] = nuSum*vMuMidMax*(fl[6]-(-1*fr[6])); 
  fjump[7] = nuSum*vMuMidMax*(fl[7]-(1*fr[7])); 
  fjump[8] = nuSum*vMuMidMax*(fl[8]-(1*fr[8])); 

  double alphaDrag[3]; 
  alphaDrag[0] = 1.414213562373095*wl[1]*nuSum+0.7071067811865475*dxvl[1]*nuSum-1.0*sumNuUx[0]; 
  alphaDrag[1] = -1.0*sumNuUx[1]; 
  alphaDrag[2] = -1.0*sumNuUx[2]; 

  double Gdiff[9]; 
  double Ghat[9]; 
  double incr2[9]; 


  incr2[2] = 0.0078125*(38.34057902536163*nuVtSqSum[2]*fr[8]+38.34057902536163*nuVtSqSum[2]*fl[8]+38.34057902536163*nuVtSqSum[1]*fr[7]+38.34057902536163*nuVtSqSum[1]*fl[7]-55.15432893255071*nuVtSqSum[2]*fr[6]+55.15432893255071*nuVtSqSum[2]*fl[6]+38.34057902536163*nuVtSqSum[0]*fr[5]+38.34057902536163*nuVtSqSum[0]*fl[5]+39.19183588453087*nuVtSqSum[2]*fr[4]+39.19183588453087*nuVtSqSum[2]*fl[4]-55.15432893255071*nuVtSqSum[1]*fr[3]+55.15432893255071*nuVtSqSum[1]*fl[3]-55.15432893255071*nuVtSqSum[0]*fr[2]+55.15432893255071*nuVtSqSum[0]*fl[2]+(39.19183588453087*fr[1]+39.19183588453087*fl[1])*nuVtSqSum[1]+(39.19183588453087*fr[0]+39.19183588453087*fl[0])*nuVtSqSum[0]); 
  incr2[3] = 0.0015625*(171.4642819948225*nuVtSqSum[1]*fr[8]+171.4642819948225*nuVtSqSum[1]*fl[8]+(171.4642819948225*nuVtSqSum[2]+191.7028951268082*nuVtSqSum[0])*fr[7]+(171.4642819948225*nuVtSqSum[2]+191.7028951268082*nuVtSqSum[0])*fl[7]-246.6576574931337*nuVtSqSum[1]*fr[6]+246.6576574931337*nuVtSqSum[1]*fl[6]+191.7028951268082*nuVtSqSum[1]*fr[5]+191.7028951268082*nuVtSqSum[1]*fl[5]+175.2712184016533*nuVtSqSum[1]*fr[4]+175.2712184016533*nuVtSqSum[1]*fl[4]+((-246.6576574931337*nuVtSqSum[2])-275.7716446627535*nuVtSqSum[0])*fr[3]+(246.6576574931337*nuVtSqSum[2]+275.7716446627535*nuVtSqSum[0])*fl[3]+(175.2712184016533*fr[1]+175.2712184016533*fl[1])*nuVtSqSum[2]-275.7716446627535*nuVtSqSum[1]*fr[2]+275.7716446627535*nuVtSqSum[1]*fl[2]+(195.9591794226543*fr[0]+195.9591794226543*fl[0])*nuVtSqSum[1]+195.9591794226543*nuVtSqSum[0]*fr[1]+195.9591794226543*nuVtSqSum[0]*fl[1]); 
  incr2[5] = -0.0078125*(148.492424049175*nuVtSqSum[2]*fr[8]+148.492424049175*nuVtSqSum[2]*fl[8]+148.492424049175*nuVtSqSum[1]*fr[7]+148.492424049175*nuVtSqSum[1]*fl[7]-213.6117974270148*nuVtSqSum[2]*fr[6]+213.6117974270148*nuVtSqSum[2]*fl[6]+148.492424049175*nuVtSqSum[0]*fr[5]+148.492424049175*nuVtSqSum[0]*fl[5]+151.7893276880823*nuVtSqSum[2]*fr[4]+151.7893276880823*nuVtSqSum[2]*fl[4]-213.6117974270148*nuVtSqSum[1]*fr[3]+213.6117974270148*nuVtSqSum[1]*fl[3]-213.6117974270148*nuVtSqSum[0]*fr[2]+213.6117974270148*nuVtSqSum[0]*fl[2]+(151.7893276880823*fr[1]+151.7893276880823*fl[1])*nuVtSqSum[1]+(151.7893276880823*fr[0]+151.7893276880823*fl[0])*nuVtSqSum[0]); 
  incr2[6] = 2.232142857142857e-4*((857.3214099741125*nuVtSqSum[2]+1341.920265887657*nuVtSqSum[0])*fr[8]+(857.3214099741125*nuVtSqSum[2]+1341.920265887657*nuVtSqSum[0])*fl[8]+1200.249973963757*nuVtSqSum[1]*fr[7]+1200.249973963757*nuVtSqSum[1]*fl[7]+((-1233.288287465668*nuVtSqSum[2])-1930.401512639275*nuVtSqSum[0])*fr[6]+(1233.288287465668*nuVtSqSum[2]+1930.401512639275*nuVtSqSum[0])*fl[6]+1341.920265887657*nuVtSqSum[2]*fr[5]+1341.920265887657*nuVtSqSum[2]*fl[5]+(876.3560920082665*nuVtSqSum[2]+1371.71425595858*nuVtSqSum[0])*fr[4]+(876.3560920082665*nuVtSqSum[2]+1371.71425595858*nuVtSqSum[0])*fl[4]-1726.603602451936*nuVtSqSum[1]*fr[3]+1726.603602451936*nuVtSqSum[1]*fl[3]+((-1930.401512639275*fr[2])+1930.401512639275*fl[2]+1371.71425595858*fr[0]+1371.71425595858*fl[0])*nuVtSqSum[2]+(1226.898528811573*fr[1]+1226.898528811573*fl[1])*nuVtSqSum[1]); 
  incr2[7] = -0.0078125*(132.815661727072*nuVtSqSum[1]*fr[8]+132.815661727072*nuVtSqSum[1]*fl[8]+(132.815661727072*nuVtSqSum[2]+148.492424049175*nuVtSqSum[0])*fr[7]+(132.815661727072*nuVtSqSum[2]+148.492424049175*nuVtSqSum[0])*fl[7]-191.0601999370879*nuVtSqSum[1]*fr[6]+191.0601999370879*nuVtSqSum[1]*fl[6]+148.492424049175*nuVtSqSum[1]*fr[5]+148.492424049175*nuVtSqSum[1]*fl[5]+135.7645019878172*nuVtSqSum[1]*fr[4]+135.7645019878172*nuVtSqSum[1]*fl[4]+((-191.0601999370879*nuVtSqSum[2])-213.6117974270148*nuVtSqSum[0])*fr[3]+(191.0601999370879*nuVtSqSum[2]+213.6117974270148*nuVtSqSum[0])*fl[3]+(135.7645019878172*fr[1]+135.7645019878172*fl[1])*nuVtSqSum[2]-213.6117974270148*nuVtSqSum[1]*fr[2]+213.6117974270148*nuVtSqSum[1]*fl[2]+(151.7893276880823*fr[0]+151.7893276880823*fl[0])*nuVtSqSum[1]+151.7893276880823*nuVtSqSum[0]*fr[1]+151.7893276880823*nuVtSqSum[0]*fl[1]); 
  incr2[8] = -0.001116071428571429*((664.0783086353599*nuVtSqSum[2]+1039.446968344225*nuVtSqSum[0])*fr[8]+(664.0783086353599*nuVtSqSum[2]+1039.446968344225*nuVtSqSum[0])*fl[8]+929.7096320895039*nuVtSqSum[1]*fr[7]+929.7096320895039*nuVtSqSum[1]*fl[7]+((-955.3009996854395*nuVtSqSum[2])-1495.282581989103*nuVtSqSum[0])*fr[6]+(955.3009996854395*nuVtSqSum[2]+1495.282581989103*nuVtSqSum[0])*fl[6]+1039.446968344225*nuVtSqSum[2]*fr[5]+1039.446968344225*nuVtSqSum[2]*fl[5]+(678.8225099390861*nuVtSqSum[2]+1062.525293816576*nuVtSqSum[0])*fr[4]+(678.8225099390861*nuVtSqSum[2]+1062.525293816576*nuVtSqSum[0])*fl[4]-1337.421399559615*nuVtSqSum[1]*fr[3]+1337.421399559615*nuVtSqSum[1]*fl[3]+((-1495.282581989103*fr[2])+1495.282581989103*fl[2]+1062.525293816576*fr[0]+1062.525293816576*fl[0])*nuVtSqSum[2]+(950.3515139147205*fr[1]+950.3515139147205*fl[1])*nuVtSqSum[1]); 


  Gdiff[0] = 0.0125*(75.89466384404115*nuVtSqSum[2]*fr[8]-75.89466384404115*nuVtSqSum[2]*fl[8]+75.89466384404115*nuVtSqSum[1]*fr[7]-75.89466384404115*nuVtSqSum[1]*fl[7]-134.7219358530748*nuVtSqSum[2]*fr[6]-134.7219358530748*nuVtSqSum[2]*fl[6]+75.89466384404115*nuVtSqSum[0]*fr[5]-75.89466384404115*nuVtSqSum[0]*fl[5]+106.0660171779821*nuVtSqSum[2]*fr[4]-106.0660171779821*nuVtSqSum[2]*fl[4]-134.7219358530748*nuVtSqSum[1]*fr[3]-134.7219358530748*nuVtSqSum[1]*fl[3]-134.7219358530748*nuVtSqSum[0]*fr[2]-134.7219358530748*nuVtSqSum[0]*fl[2]+(106.0660171779821*fr[1]-106.0660171779821*fl[1])*nuVtSqSum[1]+(106.0660171779821*fr[0]-106.0660171779821*fl[0])*nuVtSqSum[0]); 
  Gdiff[1] = 0.0025*(339.411254969543*nuVtSqSum[1]*fr[8]-339.411254969543*nuVtSqSum[1]*fl[8]+(339.411254969543*nuVtSqSum[2]+379.4733192202057*nuVtSqSum[0])*fr[7]+((-339.411254969543*nuVtSqSum[2])-379.4733192202057*nuVtSqSum[0])*fl[7]-602.4948132556829*nuVtSqSum[1]*fr[6]-602.4948132556829*nuVtSqSum[1]*fl[6]+379.4733192202058*nuVtSqSum[1]*fr[5]-379.4733192202058*nuVtSqSum[1]*fl[5]+474.3416490252571*nuVtSqSum[1]*fr[4]-474.3416490252571*nuVtSqSum[1]*fl[4]+((-602.494813255683*nuVtSqSum[2])-673.609679265374*nuVtSqSum[0])*fr[3]+((-602.494813255683*nuVtSqSum[2])-673.609679265374*nuVtSqSum[0])*fl[3]+(474.3416490252571*fr[1]-474.3416490252571*fl[1])*nuVtSqSum[2]-673.609679265374*nuVtSqSum[1]*fr[2]-673.609679265374*nuVtSqSum[1]*fl[2]+(530.3300858899107*fr[0]-530.3300858899107*fl[0])*nuVtSqSum[1]+530.3300858899107*nuVtSqSum[0]*fr[1]-530.3300858899107*nuVtSqSum[0]*fl[1]); 
  Gdiff[4] = 3.571428571428572e-4*((1697.056274847715*nuVtSqSum[2]+2656.313234541441*nuVtSqSum[0])*fr[8]+((-1697.056274847715*nuVtSqSum[2])-2656.313234541441*nuVtSqSum[0])*fl[8]+2375.878784786801*nuVtSqSum[1]*fr[7]-2375.878784786801*nuVtSqSum[1]*fl[7]+((-3012.474066278414*nuVtSqSum[2])-4715.267754857619*nuVtSqSum[0])*fr[6]+((-3012.474066278414*nuVtSqSum[2])-4715.267754857619*nuVtSqSum[0])*fl[6]+2656.313234541441*nuVtSqSum[2]*fr[5]-2656.313234541441*nuVtSqSum[2]*fl[5]+(2371.708245126285*nuVtSqSum[2]+3712.310601229375*nuVtSqSum[0])*fr[4]+((-2371.708245126285*nuVtSqSum[2])-3712.310601229375*nuVtSqSum[0])*fl[4]-4217.463692789781*nuVtSqSum[1]*fr[3]-4217.463692789781*nuVtSqSum[1]*fl[3]+((-4715.267754857618*fr[2])-4715.267754857618*fl[2]+3712.310601229375*fr[0]-3712.310601229375*fl[0])*nuVtSqSum[2]+(3320.3915431768*fr[1]-3320.3915431768*fl[1])*nuVtSqSum[1]); 

  Ghat[0] = Gdiff[0]*rdv2L+alphaDrag[2]*(0.7905694150420947*favg[8]+0.6123724356957944*favg[6]+0.3535533905932737*favg[4])+alphaDrag[1]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])-1.118033988749895*fjump[5]+alphaDrag[0]*(0.7905694150420947*favg[5]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = Gdiff[1]*rdv2L+alphaDrag[1]*(0.7071067811865475*favg[8]+0.5477225575051661*favg[6]+0.7905694150420947*favg[5]+0.3162277660168379*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])-1.118033988749895*fjump[7]+alphaDrag[0]*(0.7905694150420948*favg[7]+0.6123724356957944*favg[3]+0.3535533905932737*favg[1])+alphaDrag[2]*(0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])-0.8660254037844386*fjump[3]-0.5*fjump[1]; 
  Ghat[4] = Gdiff[4]*rdv2L-1.118033988749895*fjump[8]+alphaDrag[0]*(0.7905694150420947*favg[8]+0.6123724356957944*favg[6]+0.3535533905932737*favg[4])+alphaDrag[2]*(0.5050762722761053*favg[8]+0.3912303982179757*favg[6]+0.7905694150420947*favg[5]+0.2258769757263128*favg[4]+0.6123724356957944*favg[2]+0.3535533905932737*favg[0])+alphaDrag[1]*(0.7071067811865475*favg[7]+0.5477225575051661*favg[3]+0.3162277660168379*favg[1])-0.8660254037844387*fjump[6]-0.5*fjump[4]; 

  double incr1[9]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -1.118033988749895*Ghat[0]; 
  incr1[6] = 0.8660254037844387*Ghat[4]; 
  incr1[7] = -1.118033988749895*Ghat[1]; 
  incr1[8] = -1.118033988749895*Ghat[4]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 
  outr[8] += incr2[8]*rdvSq4R+incr1[8]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr2[5]*rdvSq4L-1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr2[7]*rdvSq4L-1.0*incr1[7]*rdv2L; 
  outl[8] += incr2[8]*rdvSq4L-1.0*incr1[8]*rdv2L; 

  return std::abs((0.7905694150420947*sumNuUx[2])/nuSum-(0.7071067811865475*sumNuUx[0])/nuSum+wl[1]); 
} 