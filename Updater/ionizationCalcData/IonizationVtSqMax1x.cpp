#include <IonizationModDecl.h> 
#include <math.h> 
void IonizationTemp1xMax_P1(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[2]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[2]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.7071067811865476*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.4714045207910317*E*elemCharge)/m_; 
     vtSqIz[1] = 0.5*vtSq[1]; 
  }
 
  else { 
     vtSqIz[0] = 1.0e-10*vtSq[0]; 
     vtSqIz[1] = 1.0e-10*vtSq[1]; 
  }
 
} 
 
void IonizationTemp1xMax_P2(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[3]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[3]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.7071067811865476*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.4714045207910317*E*elemCharge)/m_; 
     vtSqIz[1] = 0.5*vtSq[1]; 
     vtSqIz[2] = 0.5*vtSq[2]; 
  }
 
  else { 
     vtSqIz[0] = 1.0e-10*vtSq[0]; 
     vtSqIz[1] = 1.0e-10*vtSq[1]; 
     vtSqIz[2] = 1.0e-10*vtSq[2]; 
  }
 
} 
 
void IonizationTemp1xMax_P3(const double elemCharge, const double m_, const double *vtSq, const double E, double *vtSqIz) 
{ 
  // E:   	   Voronov ionization energy [eV]. 
  // m_:           mass of electron. 
  // vtSq[4]:     electron squared thermal speed, sqrt(T/m) 
  // vtSqIz[4]:   ionization squared thermal speed. 
 
  double vtSq0 = 0.7071067811865476*vtSq[0]; 
  if (vtSq0 > 2.0/3.0*E*elemCharge/m_) { 
     vtSqIz[0] = 0.5*vtSq[0]-(0.4714045207910317*E*elemCharge)/m_; 
     vtSqIz[1] = 0.5*vtSq[1]; 
     vtSqIz[2] = 0.5*vtSq[2]; 
     vtSqIz[3] = 0.5*vtSq[3]; 
  }
 
  else { 
     vtSqIz[0] = 1.0e-10*vtSq[0]; 
     vtSqIz[1] = 1.0e-10*vtSq[1]; 
     vtSqIz[2] = 1.0e-10*vtSq[2]; 
     vtSqIz[3] = 1.0e-10*vtSq[3]; 
  }
 
} 
 
