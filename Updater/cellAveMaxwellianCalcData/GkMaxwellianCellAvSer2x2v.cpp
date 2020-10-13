#include <MaxwellianCellAvModDecl.h> 
#include <math.h> 
void GkMaxwellianCellAvSer2x2v_P1(const double m_, const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax) 
{ 
  // w[4]:      cell-center coordinates. 
  // m0[4]:     particle density. 
  // uPar[8]:   fluid velocity. 
  // vtSq[4]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.5*m0[0]; 
  double uParAv = 0.5*uPar[0]; 
  double vtSqAv = 0.5*vtSq[0]; 
  double bmagAv = 0.5*bmag[0]; 
  double vSqAv = 0.5*pow(uParAv,2)-1.0*w[2]*uParAv+(w[3]*bmagAv)/m_+0.5*pow(w[2],2); 
 
  fMax[0] = (0.253974543736964*bmagAv*m0Av)/(pow(vtSqAv,3/2)*exp(vSqAv/vtSqAv)); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
void GkMaxwellianCellAvSer2x2v_P2(const double m_, const double *w, const double *m0, const double *uPar, const double *vtSq, const double *bmag, double *fMax) 
{ 
  // w[4]:      cell-center coordinates. 
  // m0[8]:     particle density. 
  // uPar[16]:   fluid velocity. 
  // vtSq[8]:   squared thermal speed, sqrt(T/m). 
  // fMax: 	cell ave Maxwellian 
 
  double m0Av = 0.5*m0[0]; 
  double uParAv = 0.5*uPar[0]; 
  double vtSqAv = 0.5*vtSq[0]; 
  double bmagAv = 0.5*bmag[0]; 
  double vSqAv = 0.5*pow(uParAv,2)-1.0*w[2]*uParAv+(w[3]*bmagAv)/m_+0.5*pow(w[2],2); 
 
  fMax[0] = (0.253974543736964*bmagAv*m0Av)/(pow(vtSqAv,3/2)*exp(vSqAv/vtSqAv)); 
 
  if (m0Av <= 0 || vtSqAv <= 0 ) { 
    fMax[0] = 0.0;
  }
} 
