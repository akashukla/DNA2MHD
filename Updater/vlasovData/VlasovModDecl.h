#ifndef VLASOV_MOD_DELC_H 
#define VLASOV_MOD_DELC_H 
extern "C" { 
void VlasovVolStream1x1vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolStream1x1vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolStream1x2vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolStream1x2vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolStream1x3vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolStream1x3vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolStream2x2vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vMax_Y_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolStream2x2vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vMax_Y_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolStream2x3vMaxP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vMax_Y_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolStream2x3vMaxP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vMax_Y_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 



 
void VlasovVolStream1x1vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vSer_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolStream1x1vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x1vSer_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolStream1x2vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vSer_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolStream1x2vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x2vSer_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolStream1x3vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vSer_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolStream1x3vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream1x3vSer_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolStream2x2vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vSer_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vSer_Y_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolStream2x2vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x2vSer_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x2vSer_Y_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolStream2x3vSerP1(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vSer_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vSer_Y_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolStream2x3vSerP2(const double *w, const double *dxv, const double *f, double *out); 
void VlasovSurfStream2x3vSer_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfStream2x3vSer_Y_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr); 



 
void VlasovVolElc1x1vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x1vMax_X_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x1vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x1vMax_X_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc1x2vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x2vMax_X_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x2vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x2vMax_X_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc1x3vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x3vMax_X_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x3vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x3vMax_X_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc2x2vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x2vMax_X_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x2vMax_Y_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x2vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x2vMax_X_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x2vMax_Y_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc2x3vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x3vMax_X_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x3vMax_Y_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x3vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x3vMax_X_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x3vMax_Y_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 



 
void VlasovVolElc1x1vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x1vSer_X_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x1vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x1vSer_X_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc1x2vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x2vSer_X_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x2vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x2vSer_X_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc1x3vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x3vSer_X_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc1x3vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc1x3vSer_X_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc2x2vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x2vSer_X_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x2vSer_Y_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x2vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x2vSer_X_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x2vSer_Y_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 


void VlasovVolElc2x3vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x3vSer_X_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x3vSer_Y_P1(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 

void VlasovVolElc2x3vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out); 
void VlasovSurfElc2x3vSer_X_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 
void VlasovSurfElc2x3vSer_Y_P2(const double *w, const double *dxv, const double *E, const double *fl, const double *fr, double *outl, double *outr); 



 
} 
#endif 
