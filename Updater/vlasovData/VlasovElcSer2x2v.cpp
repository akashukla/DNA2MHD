#include <VlasovModDecl.h> 
void VlasovVolElc2x2vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. E/f: Input electric-field/distribution function. out: Incremented output 
  double dv10 = 2/dxv[2]; 
  double dv11 = 2/dxv[3]; 
  out[3] += 0.8660254037844386*E[3]*f[5]*dv10+0.8660254037844386*E[2]*f[2]*dv10+0.8660254037844386*E[1]*f[1]*dv10+0.8660254037844386*E[0]*f[0]*dv10; 
  out[4] += 0.8660254037844386*E[3]*f[5]*dv11+0.8660254037844386*E[2]*f[2]*dv11+0.8660254037844386*E[1]*f[1]*dv11+0.8660254037844386*E[0]*f[0]*dv11; 
  out[6] += 0.8660254037844386*E[2]*f[5]*dv10+0.8660254037844386*f[2]*E[3]*dv10+0.8660254037844386*E[0]*f[1]*dv10+0.8660254037844386*f[0]*E[1]*dv10; 
  out[7] += 0.8660254037844386*E[1]*f[5]*dv10+0.8660254037844386*f[1]*E[3]*dv10+0.8660254037844386*E[0]*f[2]*dv10+0.8660254037844386*f[0]*E[2]*dv10; 
  out[8] += 0.8660254037844386*E[2]*f[5]*dv11+0.8660254037844386*f[2]*E[3]*dv11+0.8660254037844386*E[0]*f[1]*dv11+0.8660254037844386*f[0]*E[1]*dv11; 
  out[9] += 0.8660254037844386*E[1]*f[5]*dv11+0.8660254037844386*f[1]*E[3]*dv11+0.8660254037844386*E[0]*f[2]*dv11+0.8660254037844386*f[0]*E[2]*dv11; 
  out[10] += 0.8660254037844386*E[3]*f[11]*dv11+0.8660254037844386*E[2]*f[7]*dv11+0.8660254037844386*E[1]*f[6]*dv11+0.8660254037844386*E[0]*f[3]*dv11+0.8660254037844386*E[3]*f[12]*dv10+0.8660254037844386*E[2]*f[9]*dv10+0.8660254037844386*E[1]*f[8]*dv10+0.8660254037844386*E[0]*f[4]*dv10; 
  out[11] += 0.8660254037844386*E[0]*f[5]*dv10+0.8660254037844386*f[0]*E[3]*dv10+0.8660254037844386*E[1]*f[2]*dv10+0.8660254037844386*f[1]*E[2]*dv10; 
  out[12] += 0.8660254037844386*E[0]*f[5]*dv11+0.8660254037844386*f[0]*E[3]*dv11+0.8660254037844386*E[1]*f[2]*dv11+0.8660254037844386*f[1]*E[2]*dv11; 
  out[13] += 0.8660254037844386*E[2]*f[11]*dv11+0.8660254037844386*E[3]*f[7]*dv11+0.8660254037844386*E[0]*f[6]*dv11+0.8660254037844386*E[1]*f[3]*dv11+0.8660254037844386*E[2]*f[12]*dv10+0.8660254037844386*E[3]*f[9]*dv10+0.8660254037844386*E[0]*f[8]*dv10+0.8660254037844386*E[1]*f[4]*dv10; 
  out[14] += 0.8660254037844386*E[1]*f[11]*dv11+0.8660254037844386*E[0]*f[7]*dv11+0.8660254037844386*E[3]*f[6]*dv11+0.8660254037844386*E[2]*f[3]*dv11+0.8660254037844386*E[1]*f[12]*dv10+0.8660254037844386*E[0]*f[9]*dv10+0.8660254037844386*E[3]*f[8]*dv10+0.8660254037844386*E[2]*f[4]*dv10; 
  out[15] += 0.8660254037844386*E[0]*f[11]*dv11+0.8660254037844386*E[1]*f[7]*dv11+0.8660254037844386*E[2]*f[6]*dv11+0.8660254037844386*E[3]*f[3]*dv11+0.8660254037844386*E[0]*f[12]*dv10+0.8660254037844386*E[1]*f[9]*dv10+0.8660254037844386*E[2]*f[8]*dv10+0.8660254037844386*E[3]*f[4]*dv10; 
} 
void VlasovVolElc2x2vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. E/f: Input electric-field/distribution function. out: Incremented output 
  double dv10 = 2/dxv[2]; 
  double dv11 = 2/dxv[3]; 
  out[3] += 0.8660254037844386*E[7]*f[20]*dv10+0.8660254037844386*E[6]*f[19]*dv10+0.8660254037844386*E[5]*f[12]*dv10+0.8660254037844386*E[4]*f[11]*dv10+0.8660254037844386*E[3]*f[5]*dv10+0.8660254037844386*E[2]*f[2]*dv10+0.8660254037844386*E[1]*f[1]*dv10+0.8660254037844386*E[0]*f[0]*dv10; 
  out[4] += 0.8660254037844386*E[7]*f[20]*dv11+0.8660254037844386*E[6]*f[19]*dv11+0.8660254037844386*E[5]*f[12]*dv11+0.8660254037844386*E[4]*f[11]*dv11+0.8660254037844386*E[3]*f[5]*dv11+0.8660254037844386*E[2]*f[2]*dv11+0.8660254037844386*E[1]*f[1]*dv11+0.8660254037844386*E[0]*f[0]*dv11; 
  out[6] += 0.8660254037844386*E[5]*f[20]*dv10+0.7745966692414833*E[3]*f[19]*dv10+0.8660254037844386*E[7]*f[12]*dv10+0.7745966692414833*E[1]*f[11]*dv10+0.7745966692414833*f[5]*E[6]*dv10+0.8660254037844386*E[2]*f[5]*dv10+0.7745966692414833*f[1]*E[4]*dv10+0.8660254037844386*f[2]*E[3]*dv10+0.8660254037844386*E[0]*f[1]*dv10+0.8660254037844386*f[0]*E[1]*dv10; 
  out[7] += 0.7745966692414833*E[3]*f[20]*dv10+0.8660254037844386*E[4]*f[19]*dv10+0.7745966692414833*E[2]*f[12]*dv10+0.8660254037844386*E[6]*f[11]*dv10+0.7745966692414833*f[5]*E[7]*dv10+0.8660254037844386*E[1]*f[5]*dv10+0.7745966692414833*f[2]*E[5]*dv10+0.8660254037844386*f[1]*E[3]*dv10+0.8660254037844386*E[0]*f[2]*dv10+0.8660254037844386*f[0]*E[2]*dv10; 
  out[8] += 0.8660254037844386*E[5]*f[20]*dv11+0.7745966692414833*E[3]*f[19]*dv11+0.8660254037844386*E[7]*f[12]*dv11+0.7745966692414833*E[1]*f[11]*dv11+0.7745966692414833*f[5]*E[6]*dv11+0.8660254037844386*E[2]*f[5]*dv11+0.7745966692414833*f[1]*E[4]*dv11+0.8660254037844386*f[2]*E[3]*dv11+0.8660254037844386*E[0]*f[1]*dv11+0.8660254037844386*f[0]*E[1]*dv11; 
  out[9] += 0.7745966692414833*E[3]*f[20]*dv11+0.8660254037844386*E[4]*f[19]*dv11+0.7745966692414833*E[2]*f[12]*dv11+0.8660254037844386*E[6]*f[11]*dv11+0.7745966692414833*f[5]*E[7]*dv11+0.8660254037844386*E[1]*f[5]*dv11+0.7745966692414833*f[2]*E[5]*dv11+0.8660254037844386*f[1]*E[3]*dv11+0.8660254037844386*E[0]*f[2]*dv11+0.8660254037844386*f[0]*E[2]*dv11; 
  out[10] += 0.8660254037844386*E[7]*f[33]*dv11+0.8660254037844386*E[6]*f[32]*dv11+0.8660254037844386*E[5]*f[22]*dv11+0.8660254037844386*E[4]*f[21]*dv11+0.8660254037844386*E[3]*f[15]*dv11+0.8660254037844386*E[2]*f[7]*dv11+0.8660254037844386*E[1]*f[6]*dv11+0.8660254037844386*E[0]*f[3]*dv11+0.8660254037844386*E[7]*f[36]*dv10+0.8660254037844386*E[6]*f[35]*dv10+0.8660254037844386*E[5]*f[26]*dv10+0.8660254037844386*E[4]*f[25]*dv10+0.8660254037844386*E[3]*f[16]*dv10+0.8660254037844386*E[2]*f[9]*dv10+0.8660254037844386*E[1]*f[8]*dv10+0.8660254037844386*E[0]*f[4]*dv10; 
  out[13] += 1.936491673103709*E[7]*f[33]*dv10+1.936491673103709*E[6]*f[32]*dv10+1.936491673103709*E[5]*f[22]*dv10+1.936491673103709*E[4]*f[21]*dv10+1.936491673103709*E[3]*f[15]*dv10+1.936491673103709*E[2]*f[7]*dv10+1.936491673103709*E[1]*f[6]*dv10+1.936491673103709*E[0]*f[3]*dv10; 
  out[14] += 1.936491673103709*E[7]*f[36]*dv11+1.936491673103709*E[6]*f[35]*dv11+1.936491673103709*E[5]*f[26]*dv11+1.936491673103709*E[4]*f[25]*dv11+1.936491673103709*E[3]*f[16]*dv11+1.936491673103709*E[2]*f[9]*dv11+1.936491673103709*E[1]*f[8]*dv11+1.936491673103709*E[0]*f[4]*dv11; 
  out[15] += 0.6928203230275509*E[6]*f[20]*dv10+0.7745966692414833*E[2]*f[20]*dv10+0.6928203230275509*E[7]*f[19]*dv10+0.7745966692414833*E[1]*f[19]*dv10+0.7745966692414833*E[3]*f[12]*dv10+0.7745966692414833*E[3]*f[11]*dv10+0.7745966692414833*f[2]*E[7]*dv10+0.7745966692414833*f[1]*E[6]*dv10+0.7745966692414833*E[5]*f[5]*dv10+0.7745966692414833*E[4]*f[5]*dv10+0.8660254037844386*E[0]*f[5]*dv10+0.8660254037844386*f[0]*E[3]*dv10+0.8660254037844386*E[1]*f[2]*dv10+0.8660254037844386*f[1]*E[2]*dv10; 
  out[16] += 0.6928203230275509*E[6]*f[20]*dv11+0.7745966692414833*E[2]*f[20]*dv11+0.6928203230275509*E[7]*f[19]*dv11+0.7745966692414833*E[1]*f[19]*dv11+0.7745966692414833*E[3]*f[12]*dv11+0.7745966692414833*E[3]*f[11]*dv11+0.7745966692414833*f[2]*E[7]*dv11+0.7745966692414833*f[1]*E[6]*dv11+0.7745966692414833*E[5]*f[5]*dv11+0.7745966692414833*E[4]*f[5]*dv11+0.8660254037844386*E[0]*f[5]*dv11+0.8660254037844386*f[0]*E[3]*dv11+0.8660254037844386*E[1]*f[2]*dv11+0.8660254037844386*f[1]*E[2]*dv11; 
  out[17] += 0.8660254037844386*E[5]*f[33]*dv11+0.7745966692414833*E[3]*f[32]*dv11+0.8660254037844386*E[7]*f[22]*dv11+0.7745966692414833*E[1]*f[21]*dv11+0.7745966692414833*E[6]*f[15]*dv11+0.8660254037844386*E[2]*f[15]*dv11+0.8660254037844386*E[3]*f[7]*dv11+0.7745966692414833*E[4]*f[6]*dv11+0.8660254037844386*E[0]*f[6]*dv11+0.8660254037844386*E[1]*f[3]*dv11+0.8660254037844386*E[5]*f[36]*dv10+0.7745966692414833*E[3]*f[35]*dv10+0.8660254037844386*E[7]*f[26]*dv10+0.7745966692414833*E[1]*f[25]*dv10+0.7745966692414833*E[6]*f[16]*dv10+0.8660254037844386*E[2]*f[16]*dv10+0.8660254037844386*E[3]*f[9]*dv10+0.7745966692414833*E[4]*f[8]*dv10+0.8660254037844386*E[0]*f[8]*dv10+0.8660254037844386*E[1]*f[4]*dv10; 
  out[18] += 0.7745966692414833*E[3]*f[33]*dv11+0.8660254037844386*E[4]*f[32]*dv11+0.7745966692414833*E[2]*f[22]*dv11+0.8660254037844386*E[6]*f[21]*dv11+0.7745966692414833*E[7]*f[15]*dv11+0.8660254037844386*E[1]*f[15]*dv11+0.7745966692414833*E[5]*f[7]*dv11+0.8660254037844386*E[0]*f[7]*dv11+0.8660254037844386*E[3]*f[6]*dv11+0.8660254037844386*E[2]*f[3]*dv11+0.7745966692414833*E[3]*f[36]*dv10+0.8660254037844386*E[4]*f[35]*dv10+0.7745966692414833*E[2]*f[26]*dv10+0.8660254037844386*E[6]*f[25]*dv10+0.7745966692414833*E[7]*f[16]*dv10+0.8660254037844386*E[1]*f[16]*dv10+0.7745966692414833*E[5]*f[9]*dv10+0.8660254037844386*E[0]*f[9]*dv10+0.8660254037844386*E[3]*f[8]*dv10+0.8660254037844386*E[2]*f[4]*dv10; 
  out[21] += 0.7745966692414833*E[7]*f[20]*dv10+0.5532833351724881*E[6]*f[19]*dv10+0.8660254037844386*E[2]*f[19]*dv10+0.5532833351724881*E[4]*f[11]*dv10+0.8660254037844386*E[0]*f[11]*dv10+0.8660254037844386*f[2]*E[6]*dv10+0.7745966692414833*E[3]*f[5]*dv10+0.8660254037844386*f[0]*E[4]*dv10+0.7745966692414833*E[1]*f[1]*dv10; 
  out[22] += 0.5532833351724881*E[7]*f[20]*dv10+0.8660254037844386*E[1]*f[20]*dv10+0.7745966692414833*E[6]*f[19]*dv10+0.5532833351724881*E[5]*f[12]*dv10+0.8660254037844386*E[0]*f[12]*dv10+0.8660254037844386*f[1]*E[7]*dv10+0.7745966692414833*E[3]*f[5]*dv10+0.8660254037844386*f[0]*E[5]*dv10+0.7745966692414833*E[2]*f[2]*dv10; 
  out[23] += 1.936491673103709*E[5]*f[33]*dv10+1.732050807568877*E[3]*f[32]*dv10+1.936491673103709*E[7]*f[22]*dv10+1.732050807568877*E[1]*f[21]*dv10+1.732050807568877*E[6]*f[15]*dv10+1.936491673103709*E[2]*f[15]*dv10+1.936491673103709*E[3]*f[7]*dv10+1.732050807568877*E[4]*f[6]*dv10+1.936491673103709*E[0]*f[6]*dv10+1.936491673103709*E[1]*f[3]*dv10; 
  out[24] += 1.732050807568877*E[3]*f[33]*dv10+1.936491673103709*E[4]*f[32]*dv10+1.732050807568877*E[2]*f[22]*dv10+1.936491673103709*E[6]*f[21]*dv10+1.732050807568877*E[7]*f[15]*dv10+1.936491673103709*E[1]*f[15]*dv10+1.732050807568877*E[5]*f[7]*dv10+1.936491673103709*E[0]*f[7]*dv10+1.936491673103709*E[3]*f[6]*dv10+1.936491673103709*E[2]*f[3]*dv10; 
  out[25] += 0.7745966692414833*E[7]*f[20]*dv11+0.5532833351724881*E[6]*f[19]*dv11+0.8660254037844386*E[2]*f[19]*dv11+0.5532833351724881*E[4]*f[11]*dv11+0.8660254037844386*E[0]*f[11]*dv11+0.8660254037844386*f[2]*E[6]*dv11+0.7745966692414833*E[3]*f[5]*dv11+0.8660254037844386*f[0]*E[4]*dv11+0.7745966692414833*E[1]*f[1]*dv11; 
  out[26] += 0.5532833351724881*E[7]*f[20]*dv11+0.8660254037844386*E[1]*f[20]*dv11+0.7745966692414833*E[6]*f[19]*dv11+0.5532833351724881*E[5]*f[12]*dv11+0.8660254037844386*E[0]*f[12]*dv11+0.8660254037844386*f[1]*E[7]*dv11+0.7745966692414833*E[3]*f[5]*dv11+0.8660254037844386*f[0]*E[5]*dv11+0.7745966692414833*E[2]*f[2]*dv11; 
  out[27] += 0.8660254037844386*E[3]*f[34]*dv11+0.8660254037844386*E[2]*f[24]*dv11+0.8660254037844386*E[1]*f[23]*dv11+0.8660254037844386*E[0]*f[13]*dv11+1.936491673103709*E[7]*f[45]*dv10+1.936491673103709*E[6]*f[44]*dv10+1.936491673103709*E[5]*f[38]*dv10+1.936491673103709*E[4]*f[37]*dv10+1.936491673103709*E[3]*f[31]*dv10+1.936491673103709*E[2]*f[18]*dv10+1.936491673103709*E[1]*f[17]*dv10+1.936491673103709*E[0]*f[10]*dv10; 
  out[28] += 1.936491673103709*E[5]*f[36]*dv11+1.732050807568877*E[3]*f[35]*dv11+1.936491673103709*E[7]*f[26]*dv11+1.732050807568877*E[1]*f[25]*dv11+1.732050807568877*E[6]*f[16]*dv11+1.936491673103709*E[2]*f[16]*dv11+1.936491673103709*E[3]*f[9]*dv11+1.732050807568877*E[4]*f[8]*dv11+1.936491673103709*E[0]*f[8]*dv11+1.936491673103709*E[1]*f[4]*dv11; 
  out[29] += 1.732050807568877*E[3]*f[36]*dv11+1.936491673103709*E[4]*f[35]*dv11+1.732050807568877*E[2]*f[26]*dv11+1.936491673103709*E[6]*f[25]*dv11+1.732050807568877*E[7]*f[16]*dv11+1.936491673103709*E[1]*f[16]*dv11+1.732050807568877*E[5]*f[9]*dv11+1.936491673103709*E[0]*f[9]*dv11+1.936491673103709*E[3]*f[8]*dv11+1.936491673103709*E[2]*f[4]*dv11; 
  out[30] += 1.936491673103709*E[7]*f[45]*dv11+1.936491673103709*E[6]*f[44]*dv11+1.936491673103709*E[5]*f[38]*dv11+1.936491673103709*E[4]*f[37]*dv11+1.936491673103709*E[3]*f[31]*dv11+1.936491673103709*E[2]*f[18]*dv11+1.936491673103709*E[1]*f[17]*dv11+1.936491673103709*E[0]*f[10]*dv11+0.8660254037844386*E[3]*f[41]*dv10+0.8660254037844386*E[2]*f[29]*dv10+0.8660254037844386*E[1]*f[28]*dv10+0.8660254037844386*E[0]*f[14]*dv10; 
  out[31] += 0.6928203230275509*E[6]*f[33]*dv11+0.7745966692414833*E[2]*f[33]*dv11+0.6928203230275509*E[7]*f[32]*dv11+0.7745966692414833*E[1]*f[32]*dv11+0.7745966692414833*E[3]*f[22]*dv11+0.7745966692414833*E[3]*f[21]*dv11+0.7745966692414833*E[5]*f[15]*dv11+0.7745966692414833*E[4]*f[15]*dv11+0.8660254037844386*E[0]*f[15]*dv11+0.7745966692414833*E[7]*f[7]*dv11+0.8660254037844386*E[1]*f[7]*dv11+0.7745966692414833*E[6]*f[6]*dv11+0.8660254037844386*E[2]*f[6]*dv11+0.8660254037844386*E[3]*f[3]*dv11+0.6928203230275509*E[6]*f[36]*dv10+0.7745966692414833*E[2]*f[36]*dv10+0.6928203230275509*E[7]*f[35]*dv10+0.7745966692414833*E[1]*f[35]*dv10+0.7745966692414833*E[3]*f[26]*dv10+0.7745966692414833*E[3]*f[25]*dv10+0.7745966692414833*E[5]*f[16]*dv10+0.7745966692414833*E[4]*f[16]*dv10+0.8660254037844386*E[0]*f[16]*dv10+0.7745966692414833*E[7]*f[9]*dv10+0.8660254037844386*E[1]*f[9]*dv10+0.7745966692414833*E[6]*f[8]*dv10+0.8660254037844386*E[2]*f[8]*dv10+0.8660254037844386*E[3]*f[4]*dv10; 
  out[32] += 0.6928203230275509*E[3]*f[20]*dv10+0.7745966692414833*E[5]*f[19]*dv10+0.5532833351724881*E[4]*f[19]*dv10+0.8660254037844386*E[0]*f[19]*dv10+0.7745966692414833*E[6]*f[12]*dv10+0.5532833351724881*E[6]*f[11]*dv10+0.8660254037844386*E[2]*f[11]*dv10+0.6928203230275509*f[5]*E[7]*dv10+0.8660254037844386*f[0]*E[6]*dv10+0.7745966692414833*E[1]*f[5]*dv10+0.8660254037844386*f[2]*E[4]*dv10+0.7745966692414833*f[1]*E[3]*dv10; 
  out[33] += 0.5532833351724881*E[5]*f[20]*dv10+0.7745966692414833*E[4]*f[20]*dv10+0.8660254037844386*E[0]*f[20]*dv10+0.6928203230275509*E[3]*f[19]*dv10+0.5532833351724881*E[7]*f[12]*dv10+0.8660254037844386*E[1]*f[12]*dv10+0.7745966692414833*E[7]*f[11]*dv10+0.8660254037844386*f[0]*E[7]*dv10+0.6928203230275509*f[5]*E[6]*dv10+0.7745966692414833*E[2]*f[5]*dv10+0.8660254037844386*f[1]*E[5]*dv10+0.7745966692414833*f[2]*E[3]*dv10; 
  out[34] += 1.549193338482967*E[6]*f[33]*dv10+1.732050807568877*E[2]*f[33]*dv10+1.549193338482967*E[7]*f[32]*dv10+1.732050807568877*E[1]*f[32]*dv10+1.732050807568877*E[3]*f[22]*dv10+1.732050807568877*E[3]*f[21]*dv10+1.732050807568877*E[5]*f[15]*dv10+1.732050807568877*E[4]*f[15]*dv10+1.936491673103709*E[0]*f[15]*dv10+1.732050807568877*E[7]*f[7]*dv10+1.936491673103709*E[1]*f[7]*dv10+1.732050807568877*E[6]*f[6]*dv10+1.936491673103709*E[2]*f[6]*dv10+1.936491673103709*E[3]*f[3]*dv10; 
  out[35] += 0.6928203230275509*E[3]*f[20]*dv11+0.7745966692414833*E[5]*f[19]*dv11+0.5532833351724881*E[4]*f[19]*dv11+0.8660254037844386*E[0]*f[19]*dv11+0.7745966692414833*E[6]*f[12]*dv11+0.5532833351724881*E[6]*f[11]*dv11+0.8660254037844386*E[2]*f[11]*dv11+0.6928203230275509*f[5]*E[7]*dv11+0.8660254037844386*f[0]*E[6]*dv11+0.7745966692414833*E[1]*f[5]*dv11+0.8660254037844386*f[2]*E[4]*dv11+0.7745966692414833*f[1]*E[3]*dv11; 
  out[36] += 0.5532833351724881*E[5]*f[20]*dv11+0.7745966692414833*E[4]*f[20]*dv11+0.8660254037844386*E[0]*f[20]*dv11+0.6928203230275509*E[3]*f[19]*dv11+0.5532833351724881*E[7]*f[12]*dv11+0.8660254037844386*E[1]*f[12]*dv11+0.7745966692414833*E[7]*f[11]*dv11+0.8660254037844386*f[0]*E[7]*dv11+0.6928203230275509*f[5]*E[6]*dv11+0.7745966692414833*E[2]*f[5]*dv11+0.8660254037844386*f[1]*E[5]*dv11+0.7745966692414833*f[2]*E[3]*dv11; 
  out[37] += 0.7745966692414833*E[7]*f[33]*dv11+0.5532833351724881*E[6]*f[32]*dv11+0.8660254037844386*E[2]*f[32]*dv11+0.5532833351724881*E[4]*f[21]*dv11+0.8660254037844386*E[0]*f[21]*dv11+0.7745966692414833*E[3]*f[15]*dv11+0.8660254037844386*E[6]*f[7]*dv11+0.7745966692414833*E[1]*f[6]*dv11+0.8660254037844386*f[3]*E[4]*dv11+0.7745966692414833*E[7]*f[36]*dv10+0.5532833351724881*E[6]*f[35]*dv10+0.8660254037844386*E[2]*f[35]*dv10+0.5532833351724881*E[4]*f[25]*dv10+0.8660254037844386*E[0]*f[25]*dv10+0.7745966692414833*E[3]*f[16]*dv10+0.8660254037844386*E[6]*f[9]*dv10+0.7745966692414833*E[1]*f[8]*dv10+0.8660254037844386*E[4]*f[4]*dv10; 
  out[38] += 0.5532833351724881*E[7]*f[33]*dv11+0.8660254037844386*E[1]*f[33]*dv11+0.7745966692414833*E[6]*f[32]*dv11+0.5532833351724881*E[5]*f[22]*dv11+0.8660254037844386*E[0]*f[22]*dv11+0.7745966692414833*E[3]*f[15]*dv11+0.7745966692414833*E[2]*f[7]*dv11+0.8660254037844386*f[6]*E[7]*dv11+0.8660254037844386*f[3]*E[5]*dv11+0.5532833351724881*E[7]*f[36]*dv10+0.8660254037844386*E[1]*f[36]*dv10+0.7745966692414833*E[6]*f[35]*dv10+0.5532833351724881*E[5]*f[26]*dv10+0.8660254037844386*E[0]*f[26]*dv10+0.7745966692414833*E[3]*f[16]*dv10+0.7745966692414833*E[2]*f[9]*dv10+0.8660254037844386*E[7]*f[8]*dv10+0.8660254037844386*f[4]*E[5]*dv10; 
  out[39] += 0.7745966692414833*E[6]*f[34]*dv11+0.8660254037844386*E[2]*f[34]*dv11+0.8660254037844386*E[3]*f[24]*dv11+0.7745966692414833*E[4]*f[23]*dv11+0.8660254037844386*E[0]*f[23]*dv11+0.8660254037844386*E[1]*f[13]*dv11+1.936491673103709*E[5]*f[45]*dv10+1.732050807568877*E[3]*f[44]*dv10+1.936491673103709*E[7]*f[38]*dv10+1.732050807568877*E[1]*f[37]*dv10+1.732050807568877*E[6]*f[31]*dv10+1.936491673103709*E[2]*f[31]*dv10+1.936491673103709*E[3]*f[18]*dv10+1.732050807568877*E[4]*f[17]*dv10+1.936491673103709*E[0]*f[17]*dv10+1.936491673103709*E[1]*f[10]*dv10; 
  out[40] += 0.7745966692414833*E[7]*f[34]*dv11+0.8660254037844386*E[1]*f[34]*dv11+0.7745966692414833*E[5]*f[24]*dv11+0.8660254037844386*E[0]*f[24]*dv11+0.8660254037844386*E[3]*f[23]*dv11+0.8660254037844386*E[2]*f[13]*dv11+1.732050807568877*E[3]*f[45]*dv10+1.936491673103709*E[4]*f[44]*dv10+1.732050807568877*E[2]*f[38]*dv10+1.936491673103709*E[6]*f[37]*dv10+1.732050807568877*E[7]*f[31]*dv10+1.936491673103709*E[1]*f[31]*dv10+1.732050807568877*E[5]*f[18]*dv10+1.936491673103709*E[0]*f[18]*dv10+1.936491673103709*E[3]*f[17]*dv10+1.936491673103709*E[2]*f[10]*dv10; 
  out[41] += 1.549193338482967*E[6]*f[36]*dv11+1.732050807568877*E[2]*f[36]*dv11+1.549193338482967*E[7]*f[35]*dv11+1.732050807568877*E[1]*f[35]*dv11+1.732050807568877*E[3]*f[26]*dv11+1.732050807568877*E[3]*f[25]*dv11+1.732050807568877*E[5]*f[16]*dv11+1.732050807568877*E[4]*f[16]*dv11+1.936491673103709*E[0]*f[16]*dv11+1.732050807568877*E[7]*f[9]*dv11+1.936491673103709*E[1]*f[9]*dv11+1.732050807568877*E[6]*f[8]*dv11+1.936491673103709*E[2]*f[8]*dv11+1.936491673103709*E[3]*f[4]*dv11; 
  out[42] += 1.936491673103709*E[5]*f[45]*dv11+1.732050807568877*E[3]*f[44]*dv11+1.936491673103709*E[7]*f[38]*dv11+1.732050807568877*E[1]*f[37]*dv11+1.732050807568877*E[6]*f[31]*dv11+1.936491673103709*E[2]*f[31]*dv11+1.936491673103709*E[3]*f[18]*dv11+1.732050807568877*E[4]*f[17]*dv11+1.936491673103709*E[0]*f[17]*dv11+1.936491673103709*E[1]*f[10]*dv11+0.7745966692414833*E[6]*f[41]*dv10+0.8660254037844386*E[2]*f[41]*dv10+0.8660254037844386*E[3]*f[29]*dv10+0.7745966692414833*E[4]*f[28]*dv10+0.8660254037844386*E[0]*f[28]*dv10+0.8660254037844386*E[1]*f[14]*dv10; 
  out[43] += 1.732050807568877*E[3]*f[45]*dv11+1.936491673103709*E[4]*f[44]*dv11+1.732050807568877*E[2]*f[38]*dv11+1.936491673103709*E[6]*f[37]*dv11+1.732050807568877*E[7]*f[31]*dv11+1.936491673103709*E[1]*f[31]*dv11+1.732050807568877*E[5]*f[18]*dv11+1.936491673103709*E[0]*f[18]*dv11+1.936491673103709*E[3]*f[17]*dv11+1.936491673103709*E[2]*f[10]*dv11+0.7745966692414833*E[7]*f[41]*dv10+0.8660254037844386*E[1]*f[41]*dv10+0.7745966692414833*E[5]*f[29]*dv10+0.8660254037844386*E[0]*f[29]*dv10+0.8660254037844386*E[3]*f[28]*dv10+0.8660254037844386*E[2]*f[14]*dv10; 
  out[44] += 0.6928203230275509*E[3]*f[33]*dv11+0.7745966692414833*E[5]*f[32]*dv11+0.5532833351724881*E[4]*f[32]*dv11+0.8660254037844386*E[0]*f[32]*dv11+0.7745966692414833*E[6]*f[22]*dv11+0.5532833351724881*E[6]*f[21]*dv11+0.8660254037844386*E[2]*f[21]*dv11+0.6928203230275509*E[7]*f[15]*dv11+0.7745966692414833*E[1]*f[15]*dv11+0.8660254037844386*E[4]*f[7]*dv11+0.7745966692414833*E[3]*f[6]*dv11+0.8660254037844386*f[3]*E[6]*dv11+0.6928203230275509*E[3]*f[36]*dv10+0.7745966692414833*E[5]*f[35]*dv10+0.5532833351724881*E[4]*f[35]*dv10+0.8660254037844386*E[0]*f[35]*dv10+0.7745966692414833*E[6]*f[26]*dv10+0.5532833351724881*E[6]*f[25]*dv10+0.8660254037844386*E[2]*f[25]*dv10+0.6928203230275509*E[7]*f[16]*dv10+0.7745966692414833*E[1]*f[16]*dv10+0.8660254037844386*E[4]*f[9]*dv10+0.7745966692414833*E[3]*f[8]*dv10+0.8660254037844386*f[4]*E[6]*dv10; 
  out[45] += 0.5532833351724881*E[5]*f[33]*dv11+0.7745966692414833*E[4]*f[33]*dv11+0.8660254037844386*E[0]*f[33]*dv11+0.6928203230275509*E[3]*f[32]*dv11+0.5532833351724881*E[7]*f[22]*dv11+0.8660254037844386*E[1]*f[22]*dv11+0.7745966692414833*E[7]*f[21]*dv11+0.6928203230275509*E[6]*f[15]*dv11+0.7745966692414833*E[2]*f[15]*dv11+0.7745966692414833*E[3]*f[7]*dv11+0.8660254037844386*f[3]*E[7]*dv11+0.8660254037844386*E[5]*f[6]*dv11+0.5532833351724881*E[5]*f[36]*dv10+0.7745966692414833*E[4]*f[36]*dv10+0.8660254037844386*E[0]*f[36]*dv10+0.6928203230275509*E[3]*f[35]*dv10+0.5532833351724881*E[7]*f[26]*dv10+0.8660254037844386*E[1]*f[26]*dv10+0.7745966692414833*E[7]*f[25]*dv10+0.6928203230275509*E[6]*f[16]*dv10+0.7745966692414833*E[2]*f[16]*dv10+0.7745966692414833*E[3]*f[9]*dv10+0.8660254037844386*E[5]*f[8]*dv10+0.8660254037844386*f[4]*E[7]*dv10; 
  out[46] += 0.7745966692414833*E[5]*f[34]*dv11+0.7745966692414833*E[4]*f[34]*dv11+0.8660254037844386*E[0]*f[34]*dv11+0.7745966692414833*E[7]*f[24]*dv11+0.8660254037844386*E[1]*f[24]*dv11+0.7745966692414833*E[6]*f[23]*dv11+0.8660254037844386*E[2]*f[23]*dv11+0.8660254037844386*E[3]*f[13]*dv11+1.549193338482967*E[6]*f[45]*dv10+1.732050807568877*E[2]*f[45]*dv10+1.549193338482967*E[7]*f[44]*dv10+1.732050807568877*E[1]*f[44]*dv10+1.732050807568877*E[3]*f[38]*dv10+1.732050807568877*E[3]*f[37]*dv10+1.732050807568877*E[5]*f[31]*dv10+1.732050807568877*E[4]*f[31]*dv10+1.936491673103709*E[0]*f[31]*dv10+1.732050807568877*E[7]*f[18]*dv10+1.936491673103709*E[1]*f[18]*dv10+1.732050807568877*E[6]*f[17]*dv10+1.936491673103709*E[2]*f[17]*dv10+1.936491673103709*E[3]*f[10]*dv10; 
  out[47] += 1.549193338482967*E[6]*f[45]*dv11+1.732050807568877*E[2]*f[45]*dv11+1.549193338482967*E[7]*f[44]*dv11+1.732050807568877*E[1]*f[44]*dv11+1.732050807568877*E[3]*f[38]*dv11+1.732050807568877*E[3]*f[37]*dv11+1.732050807568877*E[5]*f[31]*dv11+1.732050807568877*E[4]*f[31]*dv11+1.936491673103709*E[0]*f[31]*dv11+1.732050807568877*E[7]*f[18]*dv11+1.936491673103709*E[1]*f[18]*dv11+1.732050807568877*E[6]*f[17]*dv11+1.936491673103709*E[2]*f[17]*dv11+1.936491673103709*E[3]*f[10]*dv11+0.7745966692414833*E[5]*f[41]*dv10+0.7745966692414833*E[4]*f[41]*dv10+0.8660254037844386*E[0]*f[41]*dv10+0.7745966692414833*E[7]*f[29]*dv10+0.8660254037844386*E[1]*f[29]*dv10+0.7745966692414833*E[6]*f[28]*dv10+0.8660254037844386*E[2]*f[28]*dv10+0.8660254037844386*E[3]*f[14]*dv10; 
} 