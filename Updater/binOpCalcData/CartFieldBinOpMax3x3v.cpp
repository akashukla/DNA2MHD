#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply3x3vMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[7]; 
  tmp[0] = 0.3535533905932737*A[3]*B[3]+0.3535533905932737*A[2]*B[2]+0.3535533905932737*A[1]*B[1]+0.3535533905932737*A[0]*B[0]; 
  tmp[1] = 0.3535533905932737*A[0]*B[1]+0.3535533905932737*B[0]*A[1]; 
  tmp[2] = 0.3535533905932737*A[0]*B[2]+0.3535533905932737*B[0]*A[2]; 
  tmp[3] = 0.3535533905932737*A[0]*B[3]+0.3535533905932737*B[0]*A[3]; 
  tmp[4] = 0.3535533905932737*A[0]*B[4]; 
  tmp[5] = 0.3535533905932737*A[0]*B[5]; 
  tmp[6] = 0.3535533905932737*A[0]*B[6]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<7; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply3x3vMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[28]; 
  tmp[0] = 0.3535533905932737*A[9]*B[24]+0.3535533905932737*A[8]*B[23]+0.3535533905932737*A[7]*B[22]+0.3535533905932737*A[6]*B[9]+0.3535533905932737*A[5]*B[8]+0.3535533905932737*A[4]*B[7]+0.3535533905932737*A[3]*B[3]+0.3535533905932737*A[2]*B[2]+0.3535533905932737*A[1]*B[1]+0.3535533905932737*A[0]*B[0]; 
  tmp[1] = 0.3162277660168379*A[1]*B[22]+0.3535533905932737*A[3]*B[8]+0.3535533905932737*A[2]*B[7]+0.3162277660168379*B[1]*A[7]+0.3535533905932737*B[3]*A[5]+0.3535533905932737*B[2]*A[4]+0.3535533905932737*A[0]*B[1]+0.3535533905932737*B[0]*A[1]; 
  tmp[2] = 0.3162277660168379*A[2]*B[23]+0.3535533905932737*A[3]*B[9]+0.3162277660168379*B[2]*A[8]+0.3535533905932737*A[1]*B[7]+0.3535533905932737*B[3]*A[6]+0.3535533905932737*B[1]*A[4]+0.3535533905932737*A[0]*B[2]+0.3535533905932737*B[0]*A[2]; 
  tmp[3] = 0.3162277660168379*A[3]*B[24]+0.3535533905932737*A[2]*B[9]+0.3162277660168379*B[3]*A[9]+0.3535533905932737*A[1]*B[8]+0.3535533905932737*B[2]*A[6]+0.3535533905932737*B[1]*A[5]+0.3535533905932737*A[0]*B[3]+0.3535533905932737*B[0]*A[3]; 
  tmp[4] = 0.3535533905932737*A[3]*B[12]+0.3535533905932737*A[2]*B[11]+0.3535533905932737*A[1]*B[10]+0.3535533905932737*A[0]*B[4]; 
  tmp[5] = 0.3535533905932737*A[3]*B[15]+0.3535533905932737*A[2]*B[14]+0.3535533905932737*A[1]*B[13]+0.3535533905932737*A[0]*B[5]; 
  tmp[6] = 0.3535533905932737*A[3]*B[19]+0.3535533905932737*A[2]*B[18]+0.3535533905932737*A[1]*B[17]+0.3535533905932737*A[0]*B[6]; 
  tmp[7] = 0.3162277660168379*A[4]*B[23]+0.3162277660168379*A[4]*B[22]+0.3535533905932737*A[5]*B[9]+0.3535533905932737*A[6]*B[8]+0.3162277660168379*B[7]*A[8]+0.3162277660168379*A[7]*B[7]+0.3535533905932737*A[0]*B[7]+0.3535533905932737*B[0]*A[4]+0.3535533905932737*A[1]*B[2]+0.3535533905932737*B[1]*A[2]; 
  tmp[8] = 0.3162277660168379*A[5]*B[24]+0.3162277660168379*A[5]*B[22]+0.3535533905932737*A[4]*B[9]+0.3162277660168379*B[8]*A[9]+0.3162277660168379*A[7]*B[8]+0.3535533905932737*A[0]*B[8]+0.3535533905932737*A[6]*B[7]+0.3535533905932737*B[0]*A[5]+0.3535533905932737*A[1]*B[3]+0.3535533905932737*B[1]*A[3]; 
  tmp[9] = 0.3162277660168379*A[6]*B[24]+0.3162277660168379*A[6]*B[23]+0.3162277660168379*A[9]*B[9]+0.3162277660168379*A[8]*B[9]+0.3535533905932737*A[0]*B[9]+0.3535533905932737*A[4]*B[8]+0.3535533905932737*A[5]*B[7]+0.3535533905932737*B[0]*A[6]+0.3535533905932737*A[2]*B[3]+0.3535533905932737*B[2]*A[3]; 
  tmp[10] = 0.3535533905932737*A[5]*B[12]+0.3535533905932737*A[4]*B[11]+0.3162277660168379*A[7]*B[10]+0.3535533905932737*A[0]*B[10]+0.3535533905932737*A[1]*B[4]; 
  tmp[11] = 0.3535533905932737*A[6]*B[12]+0.3162277660168379*A[8]*B[11]+0.3535533905932737*A[0]*B[11]+0.3535533905932737*A[4]*B[10]+0.3535533905932737*A[2]*B[4]; 
  tmp[12] = 0.3162277660168379*A[9]*B[12]+0.3535533905932737*A[0]*B[12]+0.3535533905932737*A[6]*B[11]+0.3535533905932737*A[5]*B[10]+0.3535533905932737*A[3]*B[4]; 
  tmp[13] = 0.3535533905932737*A[5]*B[15]+0.3535533905932737*A[4]*B[14]+0.3162277660168379*A[7]*B[13]+0.3535533905932737*A[0]*B[13]+0.3535533905932737*A[1]*B[5]; 
  tmp[14] = 0.3535533905932737*A[6]*B[15]+0.3162277660168379*A[8]*B[14]+0.3535533905932737*A[0]*B[14]+0.3535533905932737*A[4]*B[13]+0.3535533905932737*A[2]*B[5]; 
  tmp[15] = 0.3162277660168379*A[9]*B[15]+0.3535533905932737*A[0]*B[15]+0.3535533905932737*A[6]*B[14]+0.3535533905932737*A[5]*B[13]+0.3535533905932737*A[3]*B[5]; 
  tmp[16] = 0.3535533905932737*A[0]*B[16]; 
  tmp[17] = 0.3535533905932737*A[5]*B[19]+0.3535533905932737*A[4]*B[18]+0.3162277660168379*A[7]*B[17]+0.3535533905932737*A[0]*B[17]+0.3535533905932737*A[1]*B[6]; 
  tmp[18] = 0.3535533905932737*A[6]*B[19]+0.3162277660168379*A[8]*B[18]+0.3535533905932737*A[0]*B[18]+0.3535533905932737*A[4]*B[17]+0.3535533905932737*A[2]*B[6]; 
  tmp[19] = 0.3162277660168379*A[9]*B[19]+0.3535533905932737*A[0]*B[19]+0.3535533905932737*A[6]*B[18]+0.3535533905932737*A[5]*B[17]+0.3535533905932737*A[3]*B[6]; 
  tmp[20] = 0.3535533905932737*A[0]*B[20]; 
  tmp[21] = 0.3535533905932737*A[0]*B[21]; 
  tmp[22] = 0.2258769757263128*A[7]*B[22]+0.3535533905932737*A[0]*B[22]+0.3162277660168379*A[5]*B[8]+0.3162277660168379*A[4]*B[7]+0.3535533905932737*B[0]*A[7]+0.3162277660168379*A[1]*B[1]; 
  tmp[23] = 0.2258769757263128*A[8]*B[23]+0.3535533905932737*A[0]*B[23]+0.3162277660168379*A[6]*B[9]+0.3535533905932737*B[0]*A[8]+0.3162277660168379*A[4]*B[7]+0.3162277660168379*A[2]*B[2]; 
  tmp[24] = 0.2258769757263128*A[9]*B[24]+0.3535533905932737*A[0]*B[24]+0.3162277660168379*A[6]*B[9]+0.3535533905932737*B[0]*A[9]+0.3162277660168379*A[5]*B[8]+0.3162277660168379*A[3]*B[3]; 
  tmp[25] = 0.3535533905932737*A[0]*B[25]; 
  tmp[26] = 0.3535533905932737*A[0]*B[26]; 
  tmp[27] = 0.3535533905932737*A[0]*B[27]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<28; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpDivide3x3vMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
 
  double As[4]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(7,7); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(7);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(7);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.3535533905932737*As[0]; 
  AEM(0,1) = 0.3535533905932737*As[1]; 
  AEM(0,2) = 0.3535533905932737*As[2]; 
  AEM(0,3) = 0.3535533905932737*As[3]; 
  AEM(0,4) = 0.3535533905932737*As[1]; 
  AEM(0,5) = 0.3535533905932737*As[0]; 
  AEM(1,1) = 0.3535533905932737*As[2]; 
  AEM(1,3) = 0.3535533905932737*As[0]; 
  AEM(1,5) = 0.3535533905932737*As[3]; 
  AEM(2,1) = 0.3535533905932737*As[0]; 
 
  // Fill BEV. 
  BEV << B[0],B[1],B[2],B[3],B[4],B[5],B[6]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,7,1) = u; 
 
} 
 
void CartFieldBinOpDivide3x3vMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]+1.060660171779821*A[5]+1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]+1.060660171779821*A[5]+1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]+1.060660171779821*A[5]+1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]+1.060660171779821*A[5]+1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]-1.060660171779821*A[5]-1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]-1.060660171779821*A[5]-1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]-1.060660171779821*A[5]-1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]-1.060660171779821*A[5]-1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
 
  double As[10]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    As[4] = 0.0; 
    As[5] = 0.0; 
    As[6] = 0.0; 
    As[7] = 0.0; 
    As[8] = 0.0; 
    As[9] = 0.0; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
    As[4] = A[4]; 
    As[5] = A[5]; 
    As[6] = A[6]; 
    As[7] = A[7]; 
    As[8] = A[8]; 
    As[9] = A[9]; 
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(28,28); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(28);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(28);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.3535533905932737*As[0]; 
  AEM(0,1) = 0.3535533905932737*As[1]; 
  AEM(0,2) = 0.3535533905932737*As[2]; 
  AEM(0,3) = 0.3535533905932737*As[3]; 
  AEM(0,7) = 0.3535533905932737*As[4]; 
  AEM(0,8) = 0.3535533905932737*As[5]; 
  AEM(0,9) = 0.3535533905932737*As[6]; 
  AEM(0,10) = 0.3535533905932737*As[1]; 
  AEM(0,11) = 0.3162277660168379*As[7]+0.3535533905932737*As[0]; 
  AEM(0,12) = 0.3535533905932737*As[4]; 
  AEM(0,13) = 0.3535533905932737*As[5]; 
  AEM(0,17) = 0.3535533905932737*As[2]; 
  AEM(0,18) = 0.3535533905932737*As[3]; 
  AEM(0,20) = 0.3535533905932737*As[2]; 
  AEM(0,21) = 0.3535533905932737*As[4]; 
  AEM(0,22) = 0.3162277660168379*As[8]+0.3535533905932737*As[0]; 
  AEM(0,23) = 0.3535533905932737*As[6]; 
  AEM(0,27) = 0.3535533905932737*As[1]; 
  AEM(1,1) = 0.3535533905932737*As[3]; 
  AEM(1,2) = 0.3535533905932737*As[3]; 
  AEM(1,3) = 0.3535533905932737*As[5]; 
  AEM(1,4) = 0.3535533905932737*As[6]; 
  AEM(1,5) = 0.3162277660168379*As[9]+0.3535533905932737*As[0]; 
  AEM(1,10) = 0.3535533905932737*As[1]; 
  AEM(1,11) = 0.3535533905932737*As[2]; 
  AEM(1,16) = 0.3535533905932737*As[0]; 
  AEM(1,27) = 0.3535533905932737*As[0]; 
  AEM(2,10) = 0.3535533905932737*As[0]; 
  AEM(2,14) = 0.3535533905932737*As[4]; 
  AEM(2,15) = 0.3535533905932737*As[2]; 
  AEM(2,16) = 0.3535533905932737*As[1]; 
  AEM(2,21) = 0.3162277660168379*As[8]+0.3162277660168379*As[7]+0.3535533905932737*As[0]; 
  AEM(2,22) = 0.3535533905932737*As[6]; 
  AEM(2,23) = 0.3535533905932737*As[5]; 
  AEM(2,24) = 0.3535533905932737*As[5]; 
  AEM(2,25) = 0.3535533905932737*As[3]; 
  AEM(2,27) = 0.3535533905932737*As[1]; 
  AEM(3,3) = 0.3535533905932737*As[6]; 
  AEM(3,4) = 0.3162277660168379*As[9]+0.3162277660168379*As[7]+0.3535533905932737*As[0]; 
  AEM(3,5) = 0.3535533905932737*As[4]; 
  AEM(3,6) = 0.3535533905932737*As[6]; 
  AEM(3,8) = 0.3535533905932737*As[3]; 
  AEM(3,9) = 0.3535533905932737*As[2]; 
  AEM(3,13) = 0.3535533905932737*As[5]; 
  AEM(3,14) = 0.3535533905932737*As[4]; 
  AEM(3,15) = 0.3162277660168379*As[9]+0.3162277660168379*As[8]+0.3535533905932737*As[0]; 
  AEM(3,20) = 0.3535533905932737*As[1]; 
  AEM(4,2) = 0.3535533905932737*As[2]; 
  AEM(4,12) = 0.3535533905932737*As[3]; 
  AEM(4,23) = 0.3535533905932737*As[1]; 
  AEM(5,5) = 0.3535533905932737*As[2]; 
  AEM(5,15) = 0.3535533905932737*As[3]; 
  AEM(6,8) = 0.3535533905932737*As[1]; 
  AEM(6,18) = 0.3535533905932737*As[2]; 
  AEM(7,0) = 0.3535533905932737*As[3]; 
  AEM(7,24) = 0.3535533905932737*As[7]; 
  AEM(7,25) = 0.3162277660168379*As[1]; 
  AEM(8,3) = 0.3162277660168379*As[4]; 
  AEM(8,4) = 0.3162277660168379*As[5]; 
  AEM(8,6) = 0.3535533905932737*As[8]; 
  AEM(8,8) = 0.3162277660168379*As[2]; 
  AEM(8,13) = 0.3162277660168379*As[4]; 
  AEM(8,15) = 0.3162277660168379*As[6]; 
  AEM(8,16) = 0.3535533905932737*As[9]; 
  AEM(8,19) = 0.3162277660168379*As[3]; 
  AEM(8,24) = 0.3162277660168379*As[5]; 
  AEM(8,25) = 0.3162277660168379*As[6]; 
 
  // Fill BEV. 
  BEV << B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7],B[8],B[9],B[10],B[11],B[12],B[13],B[14],B[15],B[16],B[17],B[18],B[19],B[20],B[21],B[22],B[23],B[24],B[25],B[26],B[27]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,28,1) = u; 
 
} 
 