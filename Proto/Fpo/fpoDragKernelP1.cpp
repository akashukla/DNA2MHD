#include <math.h>
#include <fpoKernelsDecl.h>

void fpoDragKernelP1(const double dt, const double *dv, const double *f, const double *fL, const double *fR, const double *fT, const double *fB, const double *h, const double *hL, const double *hR, const double *hT, const double *hB, const int isTopEdge, const int isBotEdge, const int isLeftEdge, const int isRightEdge, double *fOut) {
  if (isLeftEdge) {
    fOut[0] = f[0] +-0.0015625*((840.0*fR[2]+840.0*f[2]-1208.370804016714*fR[1]+1208.370804016714*f[1]+858.6501033599193*fR[0]+858.6501033599193*f[0])*hR[2]+((-840.0*fR[2])-840.0*f[2]+1208.370804016714*fR[1]-1208.370804016714*f[1]-858.6501033599193*fR[0]-858.6501033599193*f[0])*h[2]+((-1491.098588289856*hR[1])-1491.098588289856*h[1]+1173.93568818739*hR[0]-1173.93568818739*h[0])*fR[2]+((-1491.098588289856*hR[1])-1491.098588289856*h[1]+1173.93568818739*hR[0]-1173.93568818739*h[0])*f[2]+(2145.0*fR[1]-2145.0*f[1]-1524.204710660612*fR[0]-1524.204710660612*f[0])*hR[1]+(2145.0*fR[1]-2145.0*f[1]-1524.204710660612*fR[0]-1524.204710660612*f[0])*h[1]+(1688.749537379655*h[0]-1688.749537379655*hR[0])*fR[1]+(1688.749537379655*hR[0]-1688.749537379655*h[0])*f[1]+(1200.0*fR[0]+1200.0*f[0])*hR[0]+((-1200.0*fR[0])-1200.0*f[0])*h[0])*gkyl_ipow(dv[0],-2)*dt;
    fOut[1] = f[1] +-0.0015625*((1454.922678357857*fR[2]+1454.922678357857*f[2]-2092.959626939803*fR[1]+2092.959626939803*f[1]+1487.225604943648*fR[0]+1487.225604943648*f[0])*hR[2]+((-1454.922678357857*fR[2])-1454.922678357857*f[2]+2092.959626939803*fR[1]-2092.959626939803*f[1]-1487.225604943648*fR[0]-1487.225604943648*f[0])*h[2]+((-2582.658514012257*hR[1])-2582.658514012257*h[1]+2033.316256758894*hR[0]-2033.316256758894*h[0])*fR[2]+((-2582.658514012257*hR[1])-2582.658514012257*h[1]+2033.316256758894*hR[0]-2033.316256758894*h[0])*f[2]+(3715.248982235241*fR[1]-3715.248982235241*f[1]-2640.0*fR[0]-2640.0*f[0])*hR[1]+(3715.248982235241*fR[1]-3715.248982235241*f[1]-2640.0*fR[0]-2640.0*f[0])*h[1]+(2925.0*h[0]-2925.0*hR[0])*fR[1]+(2925.0*hR[0]-2925.0*h[0])*f[1]+(2078.460969082652*fR[0]+2078.460969082652*f[0])*hR[0]+((-2078.460969082652*fR[0])-2078.460969082652*f[0])*h[0])*gkyl_ipow(dv[0],-2)*dt;
    fOut[2] = f[2] +0.0;
    fOut[3] = f[3] +0.0;
  } else if (isRightEdge) {
    fOut[0] = f[0] +-0.0015625*((840.0*fL[2]+840.0*f[2]+1208.370804016714*fL[1]-1208.370804016714*f[1]+858.6501033599193*fL[0]+858.6501033599193*f[0])*hL[2]+((-840.0*fL[2])-840.0*f[2]-1208.370804016714*fL[1]+1208.370804016714*f[1]-858.6501033599193*fL[0]-858.6501033599193*f[0])*h[2]+(1491.098588289856*hL[1]+1491.098588289856*h[1]+1173.93568818739*hL[0]-1173.93568818739*h[0])*fL[2]+(1491.098588289856*hL[1]+1491.098588289856*h[1]+1173.93568818739*hL[0]-1173.93568818739*h[0])*f[2]+(2145.0*fL[1]-2145.0*f[1]+1524.204710660612*fL[0]+1524.204710660612*f[0])*hL[1]+(2145.0*fL[1]-2145.0*f[1]+1524.204710660612*fL[0]+1524.204710660612*f[0])*h[1]+(1688.749537379655*hL[0]-1688.749537379655*h[0])*fL[1]+(1688.749537379655*h[0]-1688.749537379655*hL[0])*f[1]+(1200.0*fL[0]+1200.0*f[0])*hL[0]+((-1200.0*fL[0])-1200.0*f[0])*h[0])*gkyl_ipow(dv[0],-2)*dt;
    fOut[1] = f[1] +0.0015625*((1454.922678357857*fL[2]+1454.922678357857*f[2]+2092.959626939803*fL[1]-2092.959626939803*f[1]+1487.225604943648*fL[0]+1487.225604943648*f[0])*hL[2]+((-1454.922678357857*fL[2])-1454.922678357857*f[2]-2092.959626939803*fL[1]+2092.959626939803*f[1]-1487.225604943648*fL[0]-1487.225604943648*f[0])*h[2]+(2582.658514012257*hL[1]+2582.658514012257*h[1]+2033.316256758894*hL[0]-2033.316256758894*h[0])*fL[2]+(2582.658514012257*hL[1]+2582.658514012257*h[1]+2033.316256758894*hL[0]-2033.316256758894*h[0])*f[2]+(3715.248982235241*fL[1]-3715.248982235241*f[1]+2640.0*fL[0]+2640.0*f[0])*hL[1]+(3715.248982235241*fL[1]-3715.248982235241*f[1]+2640.0*fL[0]+2640.0*f[0])*h[1]+(2925.0*hL[0]-2925.0*h[0])*fL[1]+(2925.0*h[0]-2925.0*hL[0])*f[1]+(2078.460969082652*fL[0]+2078.460969082652*f[0])*hL[0]+((-2078.460969082652*fL[0])-2078.460969082652*f[0])*h[0])*gkyl_ipow(dv[0],-2)*dt;
    fOut[2] = f[2] +0.0;
    fOut[3] = f[3] +0.0;
  } else {
    fOut[0] = f[0] +-0.0015625*((840.0*fR[2]+840.0*f[2]-1208.370804016714*fR[1]+1208.370804016714*f[1]+858.6501033599193*fR[0]+858.6501033599193*f[0])*hR[2]+(840.0*fL[2]+840.0*f[2]+1208.370804016714*fL[1]-1208.370804016714*f[1]+858.6501033599193*fL[0]+858.6501033599193*f[0])*hL[2]+((-840.0*fR[2])-840.0*fL[2]-1680.0*f[2]+1208.370804016714*fR[1]-1208.370804016714*fL[1]-858.6501033599193*fR[0]-858.6501033599193*fL[0]-1717.300206719839*f[0])*h[2]+((-1491.098588289856*hR[1])-1491.098588289856*h[1]+1173.93568818739*hR[0]-1173.93568818739*h[0])*fR[2]+(1491.098588289856*hL[1]+1491.098588289856*h[1]+1173.93568818739*hL[0]-1173.93568818739*h[0])*fL[2]+((-1491.098588289856*hR[1])+1491.098588289856*hL[1]+1173.93568818739*hR[0]+1173.93568818739*hL[0]-2347.87137637478*h[0])*f[2]+(2145.0*fR[1]-2145.0*f[1]-1524.204710660612*fR[0]-1524.204710660612*f[0])*hR[1]+(2145.0*fL[1]-2145.0*f[1]+1524.204710660612*fL[0]+1524.204710660612*f[0])*hL[1]+(2145.0*fR[1]+2145.0*fL[1]-4290.0*f[1]-1524.204710660612*fR[0]+1524.204710660612*fL[0])*h[1]+(1688.749537379655*h[0]-1688.749537379655*hR[0])*fR[1]+(1688.749537379655*hL[0]-1688.749537379655*h[0])*fL[1]+(1688.749537379655*hR[0]-1688.749537379655*hL[0])*f[1]+(1200.0*fR[0]+1200.0*f[0])*hR[0]+(1200.0*fL[0]+1200.0*f[0])*hL[0]+((-1200.0*fR[0])-1200.0*fL[0]-2400.0*f[0])*h[0])*gkyl_ipow(dv[0],-2)*dt;
    fOut[1] = f[1] +-0.0015625*((1454.922678357857*fR[2]+1454.922678357857*f[2]-2092.959626939803*fR[1]+2092.959626939803*f[1]+1487.225604943648*fR[0]+1487.225604943648*f[0])*hR[2]+((-1454.922678357857*fL[2])-1454.922678357857*f[2]-2092.959626939803*fL[1]+2092.959626939803*f[1]-1487.225604943648*fL[0]-1487.225604943648*f[0])*hL[2]+((-1454.922678357857*fR[2])+1454.922678357857*fL[2]+2092.959626939803*fR[1]+2092.959626939803*fL[1]-4185.919253879607*f[1]-1487.225604943648*fR[0]+1487.225604943648*fL[0])*h[2]+((-2582.658514012257*hR[1])-2582.658514012257*h[1]+2033.316256758894*hR[0]-2033.316256758894*h[0])*fR[2]+((-2582.658514012257*hL[1])-2582.658514012257*h[1]-2033.316256758894*hL[0]+2033.316256758894*h[0])*fL[2]+((-2582.658514012257*hR[1])-2582.658514012257*hL[1]-5165.317028024515*h[1]+2033.316256758894*hR[0]-2033.316256758894*hL[0])*f[2]+(3715.248982235241*fR[1]-3715.248982235241*f[1]-2640.0*fR[0]-2640.0*f[0])*hR[1]+((-3715.248982235241*fL[1])+3715.248982235241*f[1]-2640.0*fL[0]-2640.0*f[0])*hL[1]+(3715.248982235241*fR[1]-3715.248982235241*fL[1]-2640.0*fR[0]-2640.0*fL[0]-5280.0*f[0])*h[1]+(2925.0*h[0]-2925.0*hR[0])*fR[1]+(2925.0*h[0]-2925.0*hL[0])*fL[1]+(2925.0*hR[0]+2925.0*hL[0]-5850.0*h[0])*f[1]+(2078.460969082652*fR[0]+2078.460969082652*f[0])*hR[0]+((-2078.460969082652*fL[0])-2078.460969082652*f[0])*hL[0]+(2078.460969082652*fL[0]-2078.460969082652*fR[0])*h[0])*gkyl_ipow(dv[0],-2)*dt;
    fOut[2] = f[2] +0.0;
    fOut[3] = f[3] +0.0;
  }

  fOut[0] += 0.0;
  fOut[1] += 0.03125*((66.40783086353598*f[1]+38.34057902536163*f[0])*hR[2]+(66.40783086353598*f[1]-38.34057902536163*f[0])*hL[2]+132.815661727072*f[1]*h[2]+((-95.53009996854395*f[1])-55.15432893255071*f[0])*hR[1]+(95.53009996854395*f[1]-55.15432893255071*f[0])*hL[1]+110.3086578651014*f[0]*h[1]+(67.8822509939086*hR[0]+67.8822509939086*hL[0]-135.7645019878172*h[0])*f[1]+39.19183588453087*f[0]*hR[0]-39.19183588453087*f[0]*hL[0])*gkyl_ipow(dv[0],-2)*dt;
  fOut[2] += 0.0;
  fOut[3] += 0.03125*((66.40783086353598*hR[2]+66.40783086353598*hL[2]+132.815661727072*h[2]-95.53009996854395*hR[1]+95.53009996854395*hL[1]+67.8822509939086*hR[0]+67.8822509939086*hL[0]-135.7645019878172*h[0])*f[3]+38.34057902536163*f[2]*hR[2]-38.34057902536163*f[2]*hL[2]+((-55.15432893255071*hR[1])-55.15432893255071*hL[1]+110.3086578651014*h[1]+39.19183588453087*hR[0]-39.19183588453087*hL[0])*f[2])*gkyl_ipow(dv[0],-2)*dt;

  if (isBotEdge) {
    fOut[0] += -0.0015625*((840.0*fT[2]+840.0*f[2]-1208.370804016714*fT[1]+1208.370804016714*f[1]+858.6501033599193*fT[0]+858.6501033599193*f[0])*hT[2]+((-840.0*fT[2])-840.0*f[2]+1208.370804016714*fT[1]-1208.370804016714*f[1]-858.6501033599193*fT[0]-858.6501033599193*f[0])*h[2]+((-1491.098588289856*hT[1])-1491.098588289856*h[1]+1173.93568818739*hT[0]-1173.93568818739*h[0])*fT[2]+((-1491.098588289856*hT[1])-1491.098588289856*h[1]+1173.93568818739*hT[0]-1173.93568818739*h[0])*f[2]+(2145.0*fT[1]-2145.0*f[1]-1524.204710660612*fT[0]-1524.204710660612*f[0])*hT[1]+(2145.0*fT[1]-2145.0*f[1]-1524.204710660612*fT[0]-1524.204710660612*f[0])*h[1]+(1688.749537379655*h[0]-1688.749537379655*hT[0])*fT[1]+(1688.749537379655*hT[0]-1688.749537379655*h[0])*f[1]+(1200.0*fT[0]+1200.0*f[0])*hT[0]+((-1200.0*fT[0])-1200.0*f[0])*h[0])*gkyl_ipow(dv[1],-2)*dt;
    fOut[1] += 0.0;
    fOut[2] += -0.0015625*((1454.922678357857*fT[2]+1454.922678357857*f[2]-2092.959626939803*fT[1]+2092.959626939803*f[1]+1487.225604943648*fT[0]+1487.225604943648*f[0])*hT[2]+((-1454.922678357857*fT[2])-1454.922678357857*f[2]+2092.959626939803*fT[1]-2092.959626939803*f[1]-1487.225604943648*fT[0]-1487.225604943648*f[0])*h[2]+((-2582.658514012257*hT[1])-2582.658514012257*h[1]+2033.316256758894*hT[0]-2033.316256758894*h[0])*fT[2]+((-2582.658514012257*hT[1])-2582.658514012257*h[1]+2033.316256758894*hT[0]-2033.316256758894*h[0])*f[2]+(3715.248982235241*fT[1]-3715.248982235241*f[1]-2640.0*fT[0]-2640.0*f[0])*hT[1]+(3715.248982235241*fT[1]-3715.248982235241*f[1]-2640.0*fT[0]-2640.0*f[0])*h[1]+(2925.0*h[0]-2925.0*hT[0])*fT[1]+(2925.0*hT[0]-2925.0*h[0])*f[1]+(2078.460969082652*fT[0]+2078.460969082652*f[0])*hT[0]+((-2078.460969082652*fT[0])-2078.460969082652*f[0])*h[0])*gkyl_ipow(dv[1],-2)*dt;
    fOut[3] += 0.0;
  } else if (isTopEdge) {
    fOut[0] += -0.0015625*((840.0*fB[2]+840.0*f[2]+1208.370804016714*fB[1]-1208.370804016714*f[1]+858.6501033599193*fB[0]+858.6501033599193*f[0])*hB[2]+((-840.0*fB[2])-840.0*f[2]-1208.370804016714*fB[1]+1208.370804016714*f[1]-858.6501033599193*fB[0]-858.6501033599193*f[0])*h[2]+(1491.098588289856*hB[1]+1491.098588289856*h[1]+1173.93568818739*hB[0]-1173.93568818739*h[0])*fB[2]+(1491.098588289856*hB[1]+1491.098588289856*h[1]+1173.93568818739*hB[0]-1173.93568818739*h[0])*f[2]+(2145.0*fB[1]-2145.0*f[1]+1524.204710660612*fB[0]+1524.204710660612*f[0])*hB[1]+(2145.0*fB[1]-2145.0*f[1]+1524.204710660612*fB[0]+1524.204710660612*f[0])*h[1]+(1688.749537379655*hB[0]-1688.749537379655*h[0])*fB[1]+(1688.749537379655*h[0]-1688.749537379655*hB[0])*f[1]+(1200.0*fB[0]+1200.0*f[0])*hB[0]+((-1200.0*fB[0])-1200.0*f[0])*h[0])*gkyl_ipow(dv[1],-2)*dt;
    fOut[1] += 0.0;
    fOut[2] += 0.0015625*((1454.922678357857*fB[2]+1454.922678357857*f[2]+2092.959626939803*fB[1]-2092.959626939803*f[1]+1487.225604943648*fB[0]+1487.225604943648*f[0])*hB[2]+((-1454.922678357857*fB[2])-1454.922678357857*f[2]-2092.959626939803*fB[1]+2092.959626939803*f[1]-1487.225604943648*fB[0]-1487.225604943648*f[0])*h[2]+(2582.658514012257*hB[1]+2582.658514012257*h[1]+2033.316256758894*hB[0]-2033.316256758894*h[0])*fB[2]+(2582.658514012257*hB[1]+2582.658514012257*h[1]+2033.316256758894*hB[0]-2033.316256758894*h[0])*f[2]+(3715.248982235241*fB[1]-3715.248982235241*f[1]+2640.0*fB[0]+2640.0*f[0])*hB[1]+(3715.248982235241*fB[1]-3715.248982235241*f[1]+2640.0*fB[0]+2640.0*f[0])*h[1]+(2925.0*hB[0]-2925.0*h[0])*fB[1]+(2925.0*h[0]-2925.0*hB[0])*f[1]+(2078.460969082652*fB[0]+2078.460969082652*f[0])*hB[0]+((-2078.460969082652*fB[0])-2078.460969082652*f[0])*h[0])*gkyl_ipow(dv[1],-2)*dt;
    fOut[3] += 0.0;
  } else {
    fOut[0] += -0.0015625*((840.0*fT[2]+840.0*f[2]-1208.370804016714*fT[1]+1208.370804016714*f[1]+858.6501033599193*fT[0]+858.6501033599193*f[0])*hT[2]+(840.0*fB[2]+840.0*f[2]+1208.370804016714*fB[1]-1208.370804016714*f[1]+858.6501033599193*fB[0]+858.6501033599193*f[0])*hB[2]+((-840.0*fT[2])-840.0*fB[2]-1680.0*f[2]+1208.370804016714*fT[1]-1208.370804016714*fB[1]-858.6501033599193*fT[0]-858.6501033599193*fB[0]-1717.300206719839*f[0])*h[2]+((-1491.098588289856*hT[1])-1491.098588289856*h[1]+1173.93568818739*hT[0]-1173.93568818739*h[0])*fT[2]+(1491.098588289856*hB[1]+1491.098588289856*h[1]+1173.93568818739*hB[0]-1173.93568818739*h[0])*fB[2]+((-1491.098588289856*hT[1])+1491.098588289856*hB[1]+1173.93568818739*hT[0]+1173.93568818739*hB[0]-2347.87137637478*h[0])*f[2]+(2145.0*fT[1]-2145.0*f[1]-1524.204710660612*fT[0]-1524.204710660612*f[0])*hT[1]+(2145.0*fB[1]-2145.0*f[1]+1524.204710660612*fB[0]+1524.204710660612*f[0])*hB[1]+(2145.0*fT[1]+2145.0*fB[1]-4290.0*f[1]-1524.204710660612*fT[0]+1524.204710660612*fB[0])*h[1]+(1688.749537379655*h[0]-1688.749537379655*hT[0])*fT[1]+(1688.749537379655*hB[0]-1688.749537379655*h[0])*fB[1]+(1688.749537379655*hT[0]-1688.749537379655*hB[0])*f[1]+(1200.0*fT[0]+1200.0*f[0])*hT[0]+(1200.0*fB[0]+1200.0*f[0])*hB[0]+((-1200.0*fT[0])-1200.0*fB[0]-2400.0*f[0])*h[0])*gkyl_ipow(dv[1],-2)*dt;
    fOut[1] += 0.0;
    fOut[2] += -0.0015625*((1454.922678357857*fT[2]+1454.922678357857*f[2]-2092.959626939803*fT[1]+2092.959626939803*f[1]+1487.225604943648*fT[0]+1487.225604943648*f[0])*hT[2]+((-1454.922678357857*fB[2])-1454.922678357857*f[2]-2092.959626939803*fB[1]+2092.959626939803*f[1]-1487.225604943648*fB[0]-1487.225604943648*f[0])*hB[2]+((-1454.922678357857*fT[2])+1454.922678357857*fB[2]+2092.959626939803*fT[1]+2092.959626939803*fB[1]-4185.919253879607*f[1]-1487.225604943648*fT[0]+1487.225604943648*fB[0])*h[2]+((-2582.658514012257*hT[1])-2582.658514012257*h[1]+2033.316256758894*hT[0]-2033.316256758894*h[0])*fT[2]+((-2582.658514012257*hB[1])-2582.658514012257*h[1]-2033.316256758894*hB[0]+2033.316256758894*h[0])*fB[2]+((-2582.658514012257*hT[1])-2582.658514012257*hB[1]-5165.317028024515*h[1]+2033.316256758894*hT[0]-2033.316256758894*hB[0])*f[2]+(3715.248982235241*fT[1]-3715.248982235241*f[1]-2640.0*fT[0]-2640.0*f[0])*hT[1]+((-3715.248982235241*fB[1])+3715.248982235241*f[1]-2640.0*fB[0]-2640.0*f[0])*hB[1]+(3715.248982235241*fT[1]-3715.248982235241*fB[1]-2640.0*fT[0]-2640.0*fB[0]-5280.0*f[0])*h[1]+(2925.0*h[0]-2925.0*hT[0])*fT[1]+(2925.0*h[0]-2925.0*hB[0])*fB[1]+(2925.0*hT[0]+2925.0*hB[0]-5850.0*h[0])*f[1]+(2078.460969082652*fT[0]+2078.460969082652*f[0])*hT[0]+((-2078.460969082652*fB[0])-2078.460969082652*f[0])*hB[0]+(2078.460969082652*fB[0]-2078.460969082652*fT[0])*h[0])*gkyl_ipow(dv[1],-2)*dt;
    fOut[3] += 0.0;
  }

  fOut[0] += 0.0;
  fOut[1] += 0.0;
  fOut[2] += 0.03125*((66.40783086353598*f[2]+38.34057902536163*f[0])*hT[2]+(66.40783086353598*f[2]-38.34057902536163*f[0])*hB[2]+132.815661727072*f[2]*h[2]+((-95.53009996854395*hT[1])+95.53009996854395*hB[1]+67.8822509939086*hT[0]+67.8822509939086*hB[0]-135.7645019878172*h[0])*f[2]-55.15432893255071*f[0]*hT[1]-55.15432893255071*f[0]*hB[1]+110.3086578651014*f[0]*h[1]+39.19183588453087*f[0]*hT[0]-39.19183588453087*f[0]*hB[0])*gkyl_ipow(dv[1],-2)*dt;
  fOut[3] += 0.03125*((66.40783086353598*hT[2]+66.40783086353598*hB[2]+132.815661727072*h[2]-95.53009996854395*hT[1]+95.53009996854395*hB[1]+67.8822509939086*hT[0]+67.8822509939086*hB[0]-135.7645019878172*h[0])*f[3]+38.34057902536163*f[1]*hT[2]-38.34057902536163*f[1]*hB[2]-55.15432893255071*f[1]*hT[1]-55.15432893255071*f[1]*hB[1]+110.3086578651014*f[1]*h[1]+(39.19183588453087*hT[0]-39.19183588453087*hB[0])*f[1])*gkyl_ipow(dv[1],-2)*dt;
}