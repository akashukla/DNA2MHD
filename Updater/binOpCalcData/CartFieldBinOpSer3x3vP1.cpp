#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 

void CartFieldBinOpMultiply3x3vSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
}

void CartFieldBinOpConfPhaseMultiply3x3vSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[64]; 
  tmp[0] = 0.3535533905932738*A[7]*B[22]+0.3535533905932738*A[6]*B[9]+0.3535533905932738*A[5]*B[8]+0.3535533905932738*A[4]*B[7]+0.3535533905932738*A[3]*B[3]+0.3535533905932738*A[2]*B[2]+0.3535533905932738*A[1]*B[1]+0.3535533905932738*A[0]*B[0]; 
  tmp[1] = 0.3535533905932738*A[6]*B[22]+0.3535533905932738*A[7]*B[9]+0.3535533905932738*A[3]*B[8]+0.3535533905932738*A[2]*B[7]+0.3535533905932738*B[3]*A[5]+0.3535533905932738*B[2]*A[4]+0.3535533905932738*A[0]*B[1]+0.3535533905932738*B[0]*A[1]; 
  tmp[2] = 0.3535533905932738*A[5]*B[22]+0.3535533905932738*A[3]*B[9]+0.3535533905932738*A[7]*B[8]+0.3535533905932738*A[1]*B[7]+0.3535533905932738*B[3]*A[6]+0.3535533905932738*B[1]*A[4]+0.3535533905932738*A[0]*B[2]+0.3535533905932738*B[0]*A[2]; 
  tmp[3] = 0.3535533905932738*A[4]*B[22]+0.3535533905932738*A[2]*B[9]+0.3535533905932738*A[1]*B[8]+0.3535533905932738*A[7]*B[7]+0.3535533905932738*B[2]*A[6]+0.3535533905932738*B[1]*A[5]+0.3535533905932738*A[0]*B[3]+0.3535533905932738*B[0]*A[3]; 
  tmp[4] = 0.3535533905932738*A[7]*B[42]+0.3535533905932738*A[6]*B[25]+0.3535533905932738*A[5]*B[24]+0.3535533905932738*A[4]*B[23]+0.3535533905932738*A[3]*B[12]+0.3535533905932738*A[2]*B[11]+0.3535533905932738*A[1]*B[10]+0.3535533905932738*A[0]*B[4]; 
  tmp[5] = 0.3535533905932738*A[7]*B[43]+0.3535533905932738*A[6]*B[28]+0.3535533905932738*A[5]*B[27]+0.3535533905932738*A[4]*B[26]+0.3535533905932738*A[3]*B[15]+0.3535533905932738*A[2]*B[14]+0.3535533905932738*A[1]*B[13]+0.3535533905932738*A[0]*B[5]; 
  tmp[6] = 0.3535533905932738*A[7]*B[47]+0.3535533905932738*A[6]*B[34]+0.3535533905932738*A[5]*B[33]+0.3535533905932738*A[4]*B[32]+0.3535533905932738*A[3]*B[19]+0.3535533905932738*A[2]*B[18]+0.3535533905932738*A[1]*B[17]+0.3535533905932738*A[0]*B[6]; 
  tmp[7] = 0.3535533905932738*A[3]*B[22]+0.3535533905932738*A[5]*B[9]+0.3535533905932738*A[6]*B[8]+0.3535533905932738*A[0]*B[7]+0.3535533905932738*B[3]*A[7]+0.3535533905932738*B[0]*A[4]+0.3535533905932738*A[1]*B[2]+0.3535533905932738*B[1]*A[2]; 
  tmp[8] = 0.3535533905932738*A[2]*B[22]+0.3535533905932738*A[4]*B[9]+0.3535533905932738*A[0]*B[8]+0.3535533905932738*A[6]*B[7]+0.3535533905932738*B[2]*A[7]+0.3535533905932738*B[0]*A[5]+0.3535533905932738*A[1]*B[3]+0.3535533905932738*B[1]*A[3]; 
  tmp[9] = 0.3535533905932738*A[1]*B[22]+0.3535533905932738*A[0]*B[9]+0.3535533905932738*A[4]*B[8]+0.3535533905932738*A[5]*B[7]+0.3535533905932738*B[1]*A[7]+0.3535533905932738*B[0]*A[6]+0.3535533905932738*A[2]*B[3]+0.3535533905932738*B[2]*A[3]; 
  tmp[10] = 0.3535533905932738*A[6]*B[42]+0.3535533905932738*A[7]*B[25]+0.3535533905932738*A[3]*B[24]+0.3535533905932738*A[2]*B[23]+0.3535533905932738*A[5]*B[12]+0.3535533905932738*A[4]*B[11]+0.3535533905932738*A[0]*B[10]+0.3535533905932738*A[1]*B[4]; 
  tmp[11] = 0.3535533905932738*A[5]*B[42]+0.3535533905932738*A[3]*B[25]+0.3535533905932738*A[7]*B[24]+0.3535533905932738*A[1]*B[23]+0.3535533905932738*A[6]*B[12]+0.3535533905932738*A[0]*B[11]+0.3535533905932738*A[4]*B[10]+0.3535533905932738*A[2]*B[4]; 
  tmp[12] = 0.3535533905932738*A[4]*B[42]+0.3535533905932738*A[2]*B[25]+0.3535533905932738*A[1]*B[24]+0.3535533905932738*A[7]*B[23]+0.3535533905932738*A[0]*B[12]+0.3535533905932738*A[6]*B[11]+0.3535533905932738*A[5]*B[10]+0.3535533905932738*A[3]*B[4]; 
  tmp[13] = 0.3535533905932738*A[6]*B[43]+0.3535533905932738*A[7]*B[28]+0.3535533905932738*A[3]*B[27]+0.3535533905932738*A[2]*B[26]+0.3535533905932738*A[5]*B[15]+0.3535533905932738*A[4]*B[14]+0.3535533905932738*A[0]*B[13]+0.3535533905932738*A[1]*B[5]; 
  tmp[14] = 0.3535533905932738*A[5]*B[43]+0.3535533905932738*A[3]*B[28]+0.3535533905932738*A[7]*B[27]+0.3535533905932738*A[1]*B[26]+0.3535533905932738*A[6]*B[15]+0.3535533905932738*A[0]*B[14]+0.3535533905932738*A[4]*B[13]+0.3535533905932738*A[2]*B[5]; 
  tmp[15] = 0.3535533905932738*A[4]*B[43]+0.3535533905932738*A[2]*B[28]+0.3535533905932738*A[1]*B[27]+0.3535533905932738*A[7]*B[26]+0.3535533905932738*A[0]*B[15]+0.3535533905932738*A[6]*B[14]+0.3535533905932738*A[5]*B[13]+0.3535533905932738*A[3]*B[5]; 
  tmp[16] = 0.3535533905932738*A[7]*B[57]+0.3535533905932738*A[6]*B[46]+0.3535533905932738*A[5]*B[45]+0.3535533905932738*A[4]*B[44]+0.3535533905932738*A[3]*B[31]+0.3535533905932738*A[2]*B[30]+0.3535533905932738*A[1]*B[29]+0.3535533905932738*A[0]*B[16]; 
  tmp[17] = 0.3535533905932738*A[6]*B[47]+0.3535533905932738*A[7]*B[34]+0.3535533905932738*A[3]*B[33]+0.3535533905932738*A[2]*B[32]+0.3535533905932738*A[5]*B[19]+0.3535533905932738*A[4]*B[18]+0.3535533905932738*A[0]*B[17]+0.3535533905932738*A[1]*B[6]; 
  tmp[18] = 0.3535533905932738*A[5]*B[47]+0.3535533905932738*A[3]*B[34]+0.3535533905932738*A[7]*B[33]+0.3535533905932738*A[1]*B[32]+0.3535533905932738*A[6]*B[19]+0.3535533905932738*A[0]*B[18]+0.3535533905932738*A[4]*B[17]+0.3535533905932738*A[2]*B[6]; 
  tmp[19] = 0.3535533905932738*A[4]*B[47]+0.3535533905932738*A[2]*B[34]+0.3535533905932738*A[1]*B[33]+0.3535533905932738*A[7]*B[32]+0.3535533905932738*A[0]*B[19]+0.3535533905932738*A[6]*B[18]+0.3535533905932738*A[5]*B[17]+0.3535533905932738*A[3]*B[6]; 
  tmp[20] = 0.3535533905932738*A[7]*B[58]+0.3535533905932738*A[6]*B[50]+0.3535533905932738*A[5]*B[49]+0.3535533905932738*A[4]*B[48]+0.3535533905932738*A[3]*B[37]+0.3535533905932738*A[2]*B[36]+0.3535533905932738*A[1]*B[35]+0.3535533905932738*A[0]*B[20]; 
  tmp[21] = 0.3535533905932738*A[7]*B[59]+0.3535533905932738*A[6]*B[53]+0.3535533905932738*A[5]*B[52]+0.3535533905932738*A[4]*B[51]+0.3535533905932738*A[3]*B[40]+0.3535533905932738*A[2]*B[39]+0.3535533905932738*A[1]*B[38]+0.3535533905932738*A[0]*B[21]; 
  tmp[22] = 0.3535533905932738*A[0]*B[22]+0.3535533905932738*A[1]*B[9]+0.3535533905932738*A[2]*B[8]+0.3535533905932738*A[3]*B[7]+0.3535533905932738*B[0]*A[7]+0.3535533905932738*B[1]*A[6]+0.3535533905932738*B[2]*A[5]+0.3535533905932738*B[3]*A[4]; 
  tmp[23] = 0.3535533905932738*A[3]*B[42]+0.3535533905932738*A[5]*B[25]+0.3535533905932738*A[6]*B[24]+0.3535533905932738*A[0]*B[23]+0.3535533905932738*A[7]*B[12]+0.3535533905932738*A[1]*B[11]+0.3535533905932738*A[2]*B[10]+0.3535533905932738*A[4]*B[4]; 
  tmp[24] = 0.3535533905932738*A[2]*B[42]+0.3535533905932738*A[4]*B[25]+0.3535533905932738*A[0]*B[24]+0.3535533905932738*A[6]*B[23]+0.3535533905932738*A[1]*B[12]+0.3535533905932738*A[7]*B[11]+0.3535533905932738*A[3]*B[10]+0.3535533905932738*B[4]*A[5]; 
  tmp[25] = 0.3535533905932738*A[1]*B[42]+0.3535533905932738*A[0]*B[25]+0.3535533905932738*A[4]*B[24]+0.3535533905932738*A[5]*B[23]+0.3535533905932738*A[2]*B[12]+0.3535533905932738*A[3]*B[11]+0.3535533905932738*A[7]*B[10]+0.3535533905932738*B[4]*A[6]; 
  tmp[26] = 0.3535533905932738*A[3]*B[43]+0.3535533905932738*A[5]*B[28]+0.3535533905932738*A[6]*B[27]+0.3535533905932738*A[0]*B[26]+0.3535533905932738*A[7]*B[15]+0.3535533905932738*A[1]*B[14]+0.3535533905932738*A[2]*B[13]+0.3535533905932738*A[4]*B[5]; 
  tmp[27] = 0.3535533905932738*A[2]*B[43]+0.3535533905932738*A[4]*B[28]+0.3535533905932738*A[0]*B[27]+0.3535533905932738*A[6]*B[26]+0.3535533905932738*A[1]*B[15]+0.3535533905932738*A[7]*B[14]+0.3535533905932738*A[3]*B[13]+0.3535533905932738*A[5]*B[5]; 
  tmp[28] = 0.3535533905932738*A[1]*B[43]+0.3535533905932738*A[0]*B[28]+0.3535533905932738*A[4]*B[27]+0.3535533905932738*A[5]*B[26]+0.3535533905932738*A[2]*B[15]+0.3535533905932738*A[3]*B[14]+0.3535533905932738*A[7]*B[13]+0.3535533905932738*B[5]*A[6]; 
  tmp[29] = 0.3535533905932738*A[6]*B[57]+0.3535533905932738*A[7]*B[46]+0.3535533905932738*A[3]*B[45]+0.3535533905932738*A[2]*B[44]+0.3535533905932738*A[5]*B[31]+0.3535533905932738*A[4]*B[30]+0.3535533905932738*A[0]*B[29]+0.3535533905932738*A[1]*B[16]; 
  tmp[30] = 0.3535533905932738*A[5]*B[57]+0.3535533905932738*A[3]*B[46]+0.3535533905932738*A[7]*B[45]+0.3535533905932738*A[1]*B[44]+0.3535533905932738*A[6]*B[31]+0.3535533905932738*A[0]*B[30]+0.3535533905932738*A[4]*B[29]+0.3535533905932738*A[2]*B[16]; 
  tmp[31] = 0.3535533905932738*A[4]*B[57]+0.3535533905932738*A[2]*B[46]+0.3535533905932738*A[1]*B[45]+0.3535533905932738*A[7]*B[44]+0.3535533905932738*A[0]*B[31]+0.3535533905932738*A[6]*B[30]+0.3535533905932738*A[5]*B[29]+0.3535533905932738*A[3]*B[16]; 
  tmp[32] = 0.3535533905932738*A[3]*B[47]+0.3535533905932738*A[5]*B[34]+0.3535533905932738*A[6]*B[33]+0.3535533905932738*A[0]*B[32]+0.3535533905932738*A[7]*B[19]+0.3535533905932738*A[1]*B[18]+0.3535533905932738*A[2]*B[17]+0.3535533905932738*A[4]*B[6]; 
  tmp[33] = 0.3535533905932738*A[2]*B[47]+0.3535533905932738*A[4]*B[34]+0.3535533905932738*A[0]*B[33]+0.3535533905932738*A[6]*B[32]+0.3535533905932738*A[1]*B[19]+0.3535533905932738*A[7]*B[18]+0.3535533905932738*A[3]*B[17]+0.3535533905932738*A[5]*B[6]; 
  tmp[34] = 0.3535533905932738*A[1]*B[47]+0.3535533905932738*A[0]*B[34]+0.3535533905932738*A[4]*B[33]+0.3535533905932738*A[5]*B[32]+0.3535533905932738*A[2]*B[19]+0.3535533905932738*A[3]*B[18]+0.3535533905932738*A[7]*B[17]+0.3535533905932738*A[6]*B[6]; 
  tmp[35] = 0.3535533905932738*A[6]*B[58]+0.3535533905932738*A[7]*B[50]+0.3535533905932738*A[3]*B[49]+0.3535533905932738*A[2]*B[48]+0.3535533905932738*A[5]*B[37]+0.3535533905932738*A[4]*B[36]+0.3535533905932738*A[0]*B[35]+0.3535533905932738*A[1]*B[20]; 
  tmp[36] = 0.3535533905932738*A[5]*B[58]+0.3535533905932738*A[3]*B[50]+0.3535533905932738*A[7]*B[49]+0.3535533905932738*A[1]*B[48]+0.3535533905932738*A[6]*B[37]+0.3535533905932738*A[0]*B[36]+0.3535533905932738*A[4]*B[35]+0.3535533905932738*A[2]*B[20]; 
  tmp[37] = 0.3535533905932738*A[4]*B[58]+0.3535533905932738*A[2]*B[50]+0.3535533905932738*A[1]*B[49]+0.3535533905932738*A[7]*B[48]+0.3535533905932738*A[0]*B[37]+0.3535533905932738*A[6]*B[36]+0.3535533905932738*A[5]*B[35]+0.3535533905932738*A[3]*B[20]; 
  tmp[38] = 0.3535533905932738*A[6]*B[59]+0.3535533905932738*A[7]*B[53]+0.3535533905932738*A[3]*B[52]+0.3535533905932738*A[2]*B[51]+0.3535533905932738*A[5]*B[40]+0.3535533905932738*A[4]*B[39]+0.3535533905932738*A[0]*B[38]+0.3535533905932738*A[1]*B[21]; 
  tmp[39] = 0.3535533905932738*A[5]*B[59]+0.3535533905932738*A[3]*B[53]+0.3535533905932738*A[7]*B[52]+0.3535533905932738*A[1]*B[51]+0.3535533905932738*A[6]*B[40]+0.3535533905932738*A[0]*B[39]+0.3535533905932738*A[4]*B[38]+0.3535533905932738*A[2]*B[21]; 
  tmp[40] = 0.3535533905932738*A[4]*B[59]+0.3535533905932738*A[2]*B[53]+0.3535533905932738*A[1]*B[52]+0.3535533905932738*A[7]*B[51]+0.3535533905932738*A[0]*B[40]+0.3535533905932738*A[6]*B[39]+0.3535533905932738*A[5]*B[38]+0.3535533905932738*A[3]*B[21]; 
  tmp[41] = 0.3535533905932738*A[7]*B[63]+0.3535533905932738*A[6]*B[62]+0.3535533905932738*A[5]*B[61]+0.3535533905932738*A[4]*B[60]+0.3535533905932738*A[3]*B[56]+0.3535533905932738*A[2]*B[55]+0.3535533905932738*A[1]*B[54]+0.3535533905932738*A[0]*B[41]; 
  tmp[42] = 0.3535533905932738*A[0]*B[42]+0.3535533905932738*A[1]*B[25]+0.3535533905932738*A[2]*B[24]+0.3535533905932738*A[3]*B[23]+0.3535533905932738*A[4]*B[12]+0.3535533905932738*A[5]*B[11]+0.3535533905932738*A[6]*B[10]+0.3535533905932738*B[4]*A[7]; 
  tmp[43] = 0.3535533905932738*A[0]*B[43]+0.3535533905932738*A[1]*B[28]+0.3535533905932738*A[2]*B[27]+0.3535533905932738*A[3]*B[26]+0.3535533905932738*A[4]*B[15]+0.3535533905932738*A[5]*B[14]+0.3535533905932738*A[6]*B[13]+0.3535533905932738*B[5]*A[7]; 
  tmp[44] = 0.3535533905932738*A[3]*B[57]+0.3535533905932738*A[5]*B[46]+0.3535533905932738*A[6]*B[45]+0.3535533905932738*A[0]*B[44]+0.3535533905932738*A[7]*B[31]+0.3535533905932738*A[1]*B[30]+0.3535533905932738*A[2]*B[29]+0.3535533905932738*A[4]*B[16]; 
  tmp[45] = 0.3535533905932738*A[2]*B[57]+0.3535533905932738*A[4]*B[46]+0.3535533905932738*A[0]*B[45]+0.3535533905932738*A[6]*B[44]+0.3535533905932738*A[1]*B[31]+0.3535533905932738*A[7]*B[30]+0.3535533905932738*A[3]*B[29]+0.3535533905932738*A[5]*B[16]; 
  tmp[46] = 0.3535533905932738*A[1]*B[57]+0.3535533905932738*A[0]*B[46]+0.3535533905932738*A[4]*B[45]+0.3535533905932738*A[5]*B[44]+0.3535533905932738*A[2]*B[31]+0.3535533905932738*A[3]*B[30]+0.3535533905932738*A[7]*B[29]+0.3535533905932738*A[6]*B[16]; 
  tmp[47] = 0.3535533905932738*A[0]*B[47]+0.3535533905932738*A[1]*B[34]+0.3535533905932738*A[2]*B[33]+0.3535533905932738*A[3]*B[32]+0.3535533905932738*A[4]*B[19]+0.3535533905932738*A[5]*B[18]+0.3535533905932738*A[6]*B[17]+0.3535533905932738*B[6]*A[7]; 
  tmp[48] = 0.3535533905932738*A[3]*B[58]+0.3535533905932738*A[5]*B[50]+0.3535533905932738*A[6]*B[49]+0.3535533905932738*A[0]*B[48]+0.3535533905932738*A[7]*B[37]+0.3535533905932738*A[1]*B[36]+0.3535533905932738*A[2]*B[35]+0.3535533905932738*A[4]*B[20]; 
  tmp[49] = 0.3535533905932738*A[2]*B[58]+0.3535533905932738*A[4]*B[50]+0.3535533905932738*A[0]*B[49]+0.3535533905932738*A[6]*B[48]+0.3535533905932738*A[1]*B[37]+0.3535533905932738*A[7]*B[36]+0.3535533905932738*A[3]*B[35]+0.3535533905932738*A[5]*B[20]; 
  tmp[50] = 0.3535533905932738*A[1]*B[58]+0.3535533905932738*A[0]*B[50]+0.3535533905932738*A[4]*B[49]+0.3535533905932738*A[5]*B[48]+0.3535533905932738*A[2]*B[37]+0.3535533905932738*A[3]*B[36]+0.3535533905932738*A[7]*B[35]+0.3535533905932738*A[6]*B[20]; 
  tmp[51] = 0.3535533905932738*A[3]*B[59]+0.3535533905932738*A[5]*B[53]+0.3535533905932738*A[6]*B[52]+0.3535533905932738*A[0]*B[51]+0.3535533905932738*A[7]*B[40]+0.3535533905932738*A[1]*B[39]+0.3535533905932738*A[2]*B[38]+0.3535533905932738*A[4]*B[21]; 
  tmp[52] = 0.3535533905932738*A[2]*B[59]+0.3535533905932738*A[4]*B[53]+0.3535533905932738*A[0]*B[52]+0.3535533905932738*A[6]*B[51]+0.3535533905932738*A[1]*B[40]+0.3535533905932738*A[7]*B[39]+0.3535533905932738*A[3]*B[38]+0.3535533905932738*A[5]*B[21]; 
  tmp[53] = 0.3535533905932738*A[1]*B[59]+0.3535533905932738*A[0]*B[53]+0.3535533905932738*A[4]*B[52]+0.3535533905932738*A[5]*B[51]+0.3535533905932738*A[2]*B[40]+0.3535533905932738*A[3]*B[39]+0.3535533905932738*A[7]*B[38]+0.3535533905932738*A[6]*B[21]; 
  tmp[54] = 0.3535533905932738*A[6]*B[63]+0.3535533905932738*A[7]*B[62]+0.3535533905932738*A[3]*B[61]+0.3535533905932738*A[2]*B[60]+0.3535533905932738*A[5]*B[56]+0.3535533905932738*A[4]*B[55]+0.3535533905932738*A[0]*B[54]+0.3535533905932738*A[1]*B[41]; 
  tmp[55] = 0.3535533905932738*A[5]*B[63]+0.3535533905932738*A[3]*B[62]+0.3535533905932738*A[7]*B[61]+0.3535533905932738*A[1]*B[60]+0.3535533905932738*A[6]*B[56]+0.3535533905932738*A[0]*B[55]+0.3535533905932738*A[4]*B[54]+0.3535533905932738*A[2]*B[41]; 
  tmp[56] = 0.3535533905932738*A[4]*B[63]+0.3535533905932738*A[2]*B[62]+0.3535533905932738*A[1]*B[61]+0.3535533905932738*A[7]*B[60]+0.3535533905932738*A[0]*B[56]+0.3535533905932738*A[6]*B[55]+0.3535533905932738*A[5]*B[54]+0.3535533905932738*A[3]*B[41]; 
  tmp[57] = 0.3535533905932738*A[0]*B[57]+0.3535533905932738*A[1]*B[46]+0.3535533905932738*A[2]*B[45]+0.3535533905932738*A[3]*B[44]+0.3535533905932738*A[4]*B[31]+0.3535533905932738*A[5]*B[30]+0.3535533905932738*A[6]*B[29]+0.3535533905932738*A[7]*B[16]; 
  tmp[58] = 0.3535533905932738*A[0]*B[58]+0.3535533905932738*A[1]*B[50]+0.3535533905932738*A[2]*B[49]+0.3535533905932738*A[3]*B[48]+0.3535533905932738*A[4]*B[37]+0.3535533905932738*A[5]*B[36]+0.3535533905932738*A[6]*B[35]+0.3535533905932738*A[7]*B[20]; 
  tmp[59] = 0.3535533905932738*A[0]*B[59]+0.3535533905932738*A[1]*B[53]+0.3535533905932738*A[2]*B[52]+0.3535533905932738*A[3]*B[51]+0.3535533905932738*A[4]*B[40]+0.3535533905932738*A[5]*B[39]+0.3535533905932738*A[6]*B[38]+0.3535533905932738*A[7]*B[21]; 
  tmp[60] = 0.3535533905932738*A[3]*B[63]+0.3535533905932738*A[5]*B[62]+0.3535533905932738*A[6]*B[61]+0.3535533905932738*A[0]*B[60]+0.3535533905932738*A[7]*B[56]+0.3535533905932738*A[1]*B[55]+0.3535533905932738*A[2]*B[54]+0.3535533905932738*A[4]*B[41]; 
  tmp[61] = 0.3535533905932738*A[2]*B[63]+0.3535533905932738*A[4]*B[62]+0.3535533905932738*A[0]*B[61]+0.3535533905932738*A[6]*B[60]+0.3535533905932738*A[1]*B[56]+0.3535533905932738*A[7]*B[55]+0.3535533905932738*A[3]*B[54]+0.3535533905932738*A[5]*B[41]; 
  tmp[62] = 0.3535533905932738*A[1]*B[63]+0.3535533905932738*A[0]*B[62]+0.3535533905932738*A[4]*B[61]+0.3535533905932738*A[5]*B[60]+0.3535533905932738*A[2]*B[56]+0.3535533905932738*A[3]*B[55]+0.3535533905932738*A[7]*B[54]+0.3535533905932738*A[6]*B[41]; 
  tmp[63] = 0.3535533905932738*A[0]*B[63]+0.3535533905932738*A[1]*B[62]+0.3535533905932738*A[2]*B[61]+0.3535533905932738*A[3]*B[60]+0.3535533905932738*A[4]*B[56]+0.3535533905932738*A[5]*B[55]+0.3535533905932738*A[6]*B[54]+0.3535533905932738*A[7]*B[41]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<64; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseDivide3x3vSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if ((-1.837117307087383*A[7])+1.060660171779821*(A[6]+A[5]+A[4])-0.6123724356957944*(A[3]+A[2]+A[1])+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.837117307087383*A[7]+1.060660171779821*A[6]-1.060660171779821*(A[5]+A[4])-0.6123724356957944*(A[3]+A[2])+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.837117307087383*A[7]-1.060660171779821*A[6]+1.060660171779821*A[5]-1.060660171779821*A[4]-0.6123724356957944*A[3]+0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
  if ((-1.837117307087383*A[7])-1.060660171779821*(A[6]+A[5])+1.060660171779821*A[4]-0.6123724356957944*A[3]+0.6123724356957944*(A[2]+A[1])+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.837117307087383*A[7]-1.060660171779821*(A[6]+A[5])+1.060660171779821*A[4]+0.6123724356957944*A[3]-0.6123724356957944*(A[2]+A[1])+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
  if ((-1.837117307087383*A[7])-1.060660171779821*A[6]+1.060660171779821*A[5]-1.060660171779821*A[4]+0.6123724356957944*A[3]-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
  if ((-1.837117307087383*A[7])+1.060660171779821*A[6]-1.060660171779821*(A[5]+A[4])+0.6123724356957944*(A[3]+A[2])-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.837117307087383*A[7]+1.060660171779821*(A[6]+A[5]+A[4])+0.6123724356957944*(A[3]+A[2]+A[1])+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
 
  double As[8]; 
  double Bs[64]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    As[4] = 0.0; 
    As[5] = 0.0; 
    As[6] = 0.0; 
    As[7] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = 0.0; 
    Bs[3] = 0.0; 
    Bs[4] = B[4]; 
    Bs[5] = B[5]; 
    Bs[6] = B[6]; 
    Bs[7] = 0.0; 
    Bs[8] = 0.0; 
    Bs[9] = 0.0; 
    Bs[10] = 0.0; 
    Bs[11] = 0.0; 
    Bs[12] = 0.0; 
    Bs[13] = 0.0; 
    Bs[14] = 0.0; 
    Bs[15] = 0.0; 
    Bs[16] = B[16]; 
    Bs[17] = 0.0; 
    Bs[18] = 0.0; 
    Bs[19] = 0.0; 
    Bs[20] = B[20]; 
    Bs[21] = B[21]; 
    Bs[22] = 0.0; 
    Bs[23] = 0.0; 
    Bs[24] = 0.0; 
    Bs[25] = 0.0; 
    Bs[26] = 0.0; 
    Bs[27] = 0.0; 
    Bs[28] = 0.0; 
    Bs[29] = 0.0; 
    Bs[30] = 0.0; 
    Bs[31] = 0.0; 
    Bs[32] = 0.0; 
    Bs[33] = 0.0; 
    Bs[34] = 0.0; 
    Bs[35] = 0.0; 
    Bs[36] = 0.0; 
    Bs[37] = 0.0; 
    Bs[38] = 0.0; 
    Bs[39] = 0.0; 
    Bs[40] = 0.0; 
    Bs[41] = B[41]; 
    Bs[42] = 0.0; 
    Bs[43] = 0.0; 
    Bs[44] = 0.0; 
    Bs[45] = 0.0; 
    Bs[46] = 0.0; 
    Bs[47] = 0.0; 
    Bs[48] = 0.0; 
    Bs[49] = 0.0; 
    Bs[50] = 0.0; 
    Bs[51] = 0.0; 
    Bs[52] = 0.0; 
    Bs[53] = 0.0; 
    Bs[54] = 0.0; 
    Bs[55] = 0.0; 
    Bs[56] = 0.0; 
    Bs[57] = 0.0; 
    Bs[58] = 0.0; 
    Bs[59] = 0.0; 
    Bs[60] = 0.0; 
    Bs[61] = 0.0; 
    Bs[62] = 0.0; 
    Bs[63] = 0.0; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
    As[4] = A[4]; 
    As[5] = A[5]; 
    As[6] = A[6]; 
    As[7] = A[7]; 
    Bs[0] = B[0]; 
    Bs[1] = B[1]; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
    Bs[4] = B[4]; 
    Bs[5] = B[5]; 
    Bs[6] = B[6]; 
    Bs[7] = B[7]; 
    Bs[8] = B[8]; 
    Bs[9] = B[9]; 
    Bs[10] = B[10]; 
    Bs[11] = B[11]; 
    Bs[12] = B[12]; 
    Bs[13] = B[13]; 
    Bs[14] = B[14]; 
    Bs[15] = B[15]; 
    Bs[16] = B[16]; 
    Bs[17] = B[17]; 
    Bs[18] = B[18]; 
    Bs[19] = B[19]; 
    Bs[20] = B[20]; 
    Bs[21] = B[21]; 
    Bs[22] = B[22]; 
    Bs[23] = B[23]; 
    Bs[24] = B[24]; 
    Bs[25] = B[25]; 
    Bs[26] = B[26]; 
    Bs[27] = B[27]; 
    Bs[28] = B[28]; 
    Bs[29] = B[29]; 
    Bs[30] = B[30]; 
    Bs[31] = B[31]; 
    Bs[32] = B[32]; 
    Bs[33] = B[33]; 
    Bs[34] = B[34]; 
    Bs[35] = B[35]; 
    Bs[36] = B[36]; 
    Bs[37] = B[37]; 
    Bs[38] = B[38]; 
    Bs[39] = B[39]; 
    Bs[40] = B[40]; 
    Bs[41] = B[41]; 
    Bs[42] = B[42]; 
    Bs[43] = B[43]; 
    Bs[44] = B[44]; 
    Bs[45] = B[45]; 
    Bs[46] = B[46]; 
    Bs[47] = B[47]; 
    Bs[48] = B[48]; 
    Bs[49] = B[49]; 
    Bs[50] = B[50]; 
    Bs[51] = B[51]; 
    Bs[52] = B[52]; 
    Bs[53] = B[53]; 
    Bs[54] = B[54]; 
    Bs[55] = B[55]; 
    Bs[56] = B[56]; 
    Bs[57] = B[57]; 
    Bs[58] = B[58]; 
    Bs[59] = B[59]; 
    Bs[60] = B[60]; 
    Bs[61] = B[61]; 
    Bs[62] = B[62]; 
    Bs[63] = B[63]; 
  } 
 
  // Fill AEM matrix. 
  data->AEM_D = Eigen::MatrixXd::Zero(64,64); 
  data->AEM_D(0,0) = 0.3535533905932737*As[0]; 
  data->AEM_D(0,1) = 0.3535533905932737*As[1]; 
  data->AEM_D(0,2) = 0.3535533905932737*As[2]; 
  data->AEM_D(0,3) = 0.3535533905932737*As[3]; 
  data->AEM_D(0,7) = 0.3535533905932737*As[4]; 
  data->AEM_D(0,8) = 0.3535533905932737*As[1]; 
  data->AEM_D(0,9) = 0.3535533905932737*As[0]; 
  data->AEM_D(0,10) = 0.3535533905932737*As[4]; 
  data->AEM_D(0,11) = 0.3535533905932737*As[5]; 
  data->AEM_D(0,15) = 0.3535533905932737*As[2]; 
  data->AEM_D(0,16) = 0.3535533905932737*As[2]; 
  data->AEM_D(0,17) = 0.3535533905932737*As[4]; 
  data->AEM_D(0,18) = 0.3535533905932737*As[0]; 
  data->AEM_D(0,19) = 0.3535533905932737*As[6]; 
  data->AEM_D(0,23) = 0.3535533905932737*As[1]; 
  data->AEM_D(0,24) = 0.3535533905932737*As[3]; 
  data->AEM_D(0,25) = 0.3535533905932737*As[5]; 
  data->AEM_D(0,26) = 0.3535533905932737*As[6]; 
  data->AEM_D(0,27) = 0.3535533905932737*As[0]; 
  data->AEM_D(0,31) = 0.3535533905932737*As[7]; 
  data->AEM_D(0,36) = 0.3535533905932737*As[0]; 
  data->AEM_D(0,45) = 0.3535533905932737*As[0]; 
  data->AEM_D(0,54) = 0.3535533905932737*As[0]; 
  data->AEM_D(0,56) = 0.3535533905932737*As[4]; 
  data->AEM_D(0,57) = 0.3535533905932737*As[2]; 
  data->AEM_D(0,58) = 0.3535533905932737*As[1]; 
  data->AEM_D(0,59) = 0.3535533905932737*As[7]; 
  data->AEM_D(0,63) = 0.3535533905932737*As[0]; 
  data->AEM_D(1,0) = 0.3535533905932737*As[5]; 
  data->AEM_D(1,1) = 0.3535533905932737*As[3]; 
  data->AEM_D(1,2) = 0.3535533905932737*As[7]; 
  data->AEM_D(1,3) = 0.3535533905932737*As[1]; 
  data->AEM_D(1,7) = 0.3535533905932737*As[6]; 
  data->AEM_D(1,8) = 0.3535533905932737*As[6]; 
  data->AEM_D(1,9) = 0.3535533905932737*As[7]; 
  data->AEM_D(1,10) = 0.3535533905932737*As[3]; 
  data->AEM_D(1,11) = 0.3535533905932737*As[2]; 
  data->AEM_D(1,15) = 0.3535533905932737*As[5]; 
  data->AEM_D(1,20) = 0.3535533905932737*As[1]; 
  data->AEM_D(1,28) = 0.3535533905932737*As[2]; 
  data->AEM_D(1,36) = 0.3535533905932737*As[3]; 
  data->AEM_D(1,45) = 0.3535533905932737*As[1]; 
  data->AEM_D(1,53) = 0.3535533905932737*As[2]; 
  data->AEM_D(1,61) = 0.3535533905932737*As[3]; 
  data->AEM_D(2,14) = 0.3535533905932737*As[1]; 
  data->AEM_D(2,22) = 0.3535533905932737*As[2]; 
  data->AEM_D(2,30) = 0.3535533905932737*As[3]; 
  data->AEM_D(2,48) = 0.3535533905932737*As[7]; 
  data->AEM_D(2,49) = 0.3535533905932737*As[6]; 
  data->AEM_D(2,50) = 0.3535533905932737*As[5]; 
  data->AEM_D(2,51) = 0.3535533905932737*As[4]; 
  data->AEM_D(2,55) = 0.3535533905932737*As[3]; 
  data->AEM_D(2,60) = 0.3535533905932737*As[4]; 
  data->AEM_D(3,4) = 0.3535533905932737*As[5]; 
  data->AEM_D(3,12) = 0.3535533905932737*As[6]; 
  data->AEM_D(3,21) = 0.3535533905932737*As[4]; 
  data->AEM_D(3,29) = 0.3535533905932737*As[5]; 
  data->AEM_D(3,37) = 0.3535533905932737*As[6]; 
  data->AEM_D(4,6) = 0.3535533905932737*As[4]; 
  data->AEM_D(4,14) = 0.3535533905932737*As[5]; 
  data->AEM_D(4,22) = 0.3535533905932737*As[6]; 
  data->AEM_D(5,20) = 0.3535533905932737*As[7]; 
  data->AEM_D(5,29) = 0.3535533905932737*As[7]; 
  data->AEM_D(5,62) = 0.3535533905932737*As[7]; 
 
  // Fill BEV. 
  data->BEV_D << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6],Bs[7],Bs[8],Bs[9],Bs[10],Bs[11],Bs[12],Bs[13],Bs[14],Bs[15],Bs[16],Bs[17],Bs[18],Bs[19],Bs[20],Bs[21],Bs[22],Bs[23],Bs[24],Bs[25],Bs[26],Bs[27],Bs[28],Bs[29],Bs[30],Bs[31],Bs[32],Bs[33],Bs[34],Bs[35],Bs[36],Bs[37],Bs[38],Bs[39],Bs[40],Bs[41],Bs[42],Bs[43],Bs[44],Bs[45],Bs[46],Bs[47],Bs[48],Bs[49],Bs[50],Bs[51],Bs[52],Bs[53],Bs[54],Bs[55],Bs[56],Bs[57],Bs[58],Bs[59],Bs[60],Bs[61],Bs[62],Bs[63]; 
 
  // Solve the system of equations. 
  data->u_D = data->AEM_D.colPivHouseholderQr().solve(data->BEV_D); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,64,1) = data->u_D; 
 
} 
 