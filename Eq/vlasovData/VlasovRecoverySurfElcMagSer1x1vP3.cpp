#include <VlasovModDecl.h> 
__host__ __device__ double VlasovRecoverySurfElcMag1x1vSer_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10l = 2/dxvl[1]; 
  double dv10r = 2/dxvr[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxvr[1], wv1 = wr[1]; 

  double Gleft[12]; 
  double Gright[12]; 
  double incr[12]; 

  Gleft[0] = 0.7071067811865475*(E0[3]*fl[8]+E0[2]*fl[4]+E0[1]*fl[1]+E0[0]*fl[0]); 
  Gleft[1] = 0.6210590034081186*(E0[2]*fl[8]+E0[3]*fl[4])+0.6324555320336759*(E0[1]*fl[4]+fl[1]*E0[2])+0.7071067811865475*(E0[0]*fl[1]+fl[0]*E0[1]); 
  Gleft[2] = 0.7071067811865475*(E0[3]*fl[10]+E0[2]*fl[6]+E0[1]*fl[3]+E0[0]*fl[2]); 
  Gleft[3] = 0.6210590034081186*(E0[2]*fl[10]+E0[3]*fl[6])+0.6324555320336759*(E0[1]*fl[6]+E0[2]*fl[3])+0.7071067811865475*(E0[0]*fl[3]+E0[1]*fl[2]); 
  Gleft[4] = (0.421637021355784*E0[3]+0.6210590034081186*E0[1])*fl[8]+(0.4517539514526256*E0[2]+0.7071067811865475*E0[0])*fl[4]+0.6210590034081186*fl[1]*E0[3]+0.7071067811865475*fl[0]*E0[2]+0.6324555320336759*E0[1]*fl[1]; 
  Gleft[5] = 0.7071067811865475*(E0[1]*fl[7]+E0[0]*fl[5]); 
  Gleft[6] = (0.421637021355784*E0[3]+0.6210590034081186*E0[1])*fl[10]+(0.4517539514526256*E0[2]+0.7071067811865475*E0[0])*fl[6]+(0.6210590034081186*E0[3]+0.6324555320336759*E0[1])*fl[3]+0.7071067811865475*E0[2]*fl[2]; 
  Gleft[7] = 0.6324555320336759*E0[2]*fl[7]+0.7071067811865475*(E0[0]*fl[7]+E0[1]*fl[5]); 
  Gleft[8] = (0.421637021355784*E0[2]+0.7071067811865475*E0[0])*fl[8]+(0.421637021355784*E0[3]+0.6210590034081186*E0[1])*fl[4]+0.7071067811865475*fl[0]*E0[3]+0.6210590034081186*fl[1]*E0[2]; 
  Gleft[9] = 0.7071067811865475*(E0[1]*fl[11]+E0[0]*fl[9]); 
  Gleft[10] = (0.421637021355784*E0[2]+0.7071067811865475*E0[0])*fl[10]+0.421637021355784*E0[3]*fl[6]+0.6210590034081186*(E0[1]*fl[6]+E0[2]*fl[3])+0.7071067811865475*fl[2]*E0[3]; 
  Gleft[11] = 0.6324555320336759*E0[2]*fl[11]+0.7071067811865475*(E0[0]*fl[11]+E0[1]*fl[9]); 

  Gright[0] = 0.7071067811865475*(E0[3]*fr[8]+E0[2]*fr[4]+E0[1]*fr[1]+E0[0]*fr[0]); 
  Gright[1] = 0.6210590034081186*(E0[2]*fr[8]+E0[3]*fr[4])+0.6324555320336759*(E0[1]*fr[4]+fr[1]*E0[2])+0.7071067811865475*(E0[0]*fr[1]+fr[0]*E0[1]); 
  Gright[2] = 0.7071067811865475*(E0[3]*fr[10]+E0[2]*fr[6]+E0[1]*fr[3]+E0[0]*fr[2]); 
  Gright[3] = 0.6210590034081186*(E0[2]*fr[10]+E0[3]*fr[6])+0.6324555320336759*(E0[1]*fr[6]+E0[2]*fr[3])+0.7071067811865475*(E0[0]*fr[3]+E0[1]*fr[2]); 
  Gright[4] = (0.421637021355784*E0[3]+0.6210590034081186*E0[1])*fr[8]+(0.4517539514526256*E0[2]+0.7071067811865475*E0[0])*fr[4]+0.6210590034081186*fr[1]*E0[3]+0.7071067811865475*fr[0]*E0[2]+0.6324555320336759*E0[1]*fr[1]; 
  Gright[5] = 0.7071067811865475*(E0[1]*fr[7]+E0[0]*fr[5]); 
  Gright[6] = (0.421637021355784*E0[3]+0.6210590034081186*E0[1])*fr[10]+(0.4517539514526256*E0[2]+0.7071067811865475*E0[0])*fr[6]+(0.6210590034081186*E0[3]+0.6324555320336759*E0[1])*fr[3]+0.7071067811865475*E0[2]*fr[2]; 
  Gright[7] = 0.6324555320336759*E0[2]*fr[7]+0.7071067811865475*(E0[0]*fr[7]+E0[1]*fr[5]); 
  Gright[8] = (0.421637021355784*E0[2]+0.7071067811865475*E0[0])*fr[8]+(0.421637021355784*E0[3]+0.6210590034081186*E0[1])*fr[4]+0.7071067811865475*fr[0]*E0[3]+0.6210590034081186*fr[1]*E0[2]; 
  Gright[9] = 0.7071067811865475*(E0[1]*fr[11]+E0[0]*fr[9]); 
  Gright[10] = (0.421637021355784*E0[2]+0.7071067811865475*E0[0])*fr[10]+0.421637021355784*E0[3]*fr[6]+0.6210590034081186*(E0[1]*fr[6]+E0[2]*fr[3])+0.7071067811865475*fr[2]*E0[3]; 
  Gright[11] = 0.6324555320336759*E0[2]*fr[11]+0.7071067811865475*(E0[0]*fr[11]+E0[1]*fr[9]); 

  incr[0] = (-0.1889822365046136*Gright[9])+0.1889822365046136*Gleft[9]+0.3493856214843422*(Gright[5]+Gleft[5])-0.3788861141556919*Gright[2]+0.3788861141556919*Gleft[2]+0.25*(Gright[0]+Gleft[0]); 
  incr[1] = (-0.1889822365046136*Gright[11])+0.1889822365046136*Gleft[11]+0.3493856214843422*(Gright[7]+Gleft[7])-0.3788861141556919*Gright[3]+0.3788861141556919*Gleft[3]+0.25*(Gright[1]+Gleft[1]); 
  incr[2] = 0.3273268353539885*Gright[9]-0.3273268353539885*Gleft[9]-0.605153647844909*(Gright[5]+Gleft[5])+0.65625*Gright[2]-0.65625*Gleft[2]-0.4330127018922193*(Gright[0]+Gleft[0]); 
  incr[3] = 0.3273268353539885*Gright[11]-0.3273268353539885*Gleft[11]-0.6051536478449089*(Gright[7]+Gleft[7])+0.65625*Gright[3]-0.65625*Gleft[3]-0.4330127018922193*(Gright[1]+Gleft[1]); 
  incr[4] = (-0.3788861141556919*Gright[6])+0.3788861141556919*Gleft[6]+0.25*(Gright[4]+Gleft[4]); 
  incr[5] = (-0.4225771273642583*Gright[9])+0.4225771273642583*Gleft[9]+0.78125*(Gright[5]+Gleft[5])-0.8472151069828725*Gright[2]+0.8472151069828725*Gleft[2]+0.5590169943749475*(Gright[0]+Gleft[0]); 
  incr[6] = 0.65625*Gright[6]-0.65625*Gleft[6]-0.4330127018922194*(Gright[4]+Gleft[4]); 
  incr[7] = (-0.4225771273642583*Gright[11])+0.4225771273642583*Gleft[11]+0.78125*(Gright[7]+Gleft[7])-0.8472151069828725*Gright[3]+0.8472151069828725*Gleft[3]+0.5590169943749476*(Gright[1]+Gleft[1]); 
  incr[8] = (-0.3788861141556919*Gright[10])+0.3788861141556919*Gleft[10]+0.25*(Gright[8]+Gleft[8]); 
  incr[9] = 0.5*Gright[9]-0.5*Gleft[9]-0.9243874661093152*(Gright[5]+Gleft[5])+1.00243843327159*Gright[2]-1.00243843327159*Gleft[2]-0.6614378277661477*(Gright[0]+Gleft[0]); 
  incr[10] = 0.65625*Gright[10]-0.65625*Gleft[10]-0.4330127018922193*(Gright[8]+Gleft[8]); 
  incr[11] = 0.5*Gright[11]-0.5*Gleft[11]-0.924387466109315*(Gright[7]+Gleft[7])+1.00243843327159*Gright[3]-1.00243843327159*Gleft[3]-0.6614378277661477*(Gright[1]+Gleft[1]); 

  outr[0] += incr[0]*dv10r; 
  outr[1] += incr[1]*dv10r; 
  outr[2] += incr[2]*dv10r; 
  outr[3] += incr[3]*dv10r; 
  outr[4] += incr[4]*dv10r; 
  outr[5] += incr[5]*dv10r; 
  outr[6] += incr[6]*dv10r; 
  outr[7] += incr[7]*dv10r; 
  outr[8] += incr[8]*dv10r; 
  outr[9] += incr[9]*dv10r; 
  outr[10] += incr[10]*dv10r; 
  outr[11] += incr[11]*dv10r; 
  outl[0] += -1.0*incr[0]*dv10l; 
  outl[1] += -1.0*incr[1]*dv10l; 
  outl[2] += incr[2]*dv10l; 
  outl[3] += incr[3]*dv10l; 
  outl[4] += -1.0*incr[4]*dv10l; 
  outl[5] += -1.0*incr[5]*dv10l; 
  outl[6] += incr[6]*dv10l; 
  outl[7] += -1.0*incr[7]*dv10l; 
  outl[8] += -1.0*incr[8]*dv10l; 
  outl[9] += incr[9]*dv10l; 
  outl[10] += incr[10]*dv10l; 
  outl[11] += incr[11]*dv10l; 
  const double amid = 0.7071067811865475*E0[0]-0.7905694150420947*E0[2]; 

  return std::abs(amid); 
} 
