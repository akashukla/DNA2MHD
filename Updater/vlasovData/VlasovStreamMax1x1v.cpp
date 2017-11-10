void VlasovVolStream1x1vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  register double dv0dx0 = dxv[1]/dxv[0]; 
  register double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
} 
void VlasovVolStream1x1vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  register double dv0dx0 = dxv[1]/dxv[0]; 
  register double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[3] += 3.464101615137754*f[2]*w0dx0+0.8944271909999159*f[5]*dv0dx0+f[0]*dv0dx0; 
  out[4] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[3]*dv0dx0; 
} 
void VlasovVolStream1x1vMaxP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  register double dv0dx0 = dxv[1]/dxv[0]; 
  register double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[3] += 3.464101615137754*f[2]*w0dx0+0.8944271909999159*f[5]*dv0dx0+f[0]*dv0dx0; 
  out[4] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[3]*dv0dx0; 
  out[6] += 7.745966692414834*f[3]*w0dx0+2.0*f[7]*dv0dx0+2.23606797749979*f[1]*dv0dx0; 
  out[7] += 3.464101615137754*f[5]*w0dx0+0.8783100656536796*f[9]*dv0dx0+0.8944271909999159*f[2]*dv0dx0; 
  out[8] += 11.83215956619923*f[4]*w0dx0+5.291502622129181*f[0]*w0dx0+3.415650255319867*f[6]*dv0dx0+1.527525231651947*f[2]*dv0dx0; 
} 
void VlasovVolStream1x1vMaxP4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  register double dv0dx0 = dxv[1]/dxv[0]; 
  register double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[3] += 3.464101615137754*f[2]*w0dx0+0.8944271909999159*f[5]*dv0dx0+f[0]*dv0dx0; 
  out[4] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[3]*dv0dx0; 
  out[6] += 7.745966692414834*f[3]*w0dx0+2.0*f[7]*dv0dx0+2.23606797749979*f[1]*dv0dx0; 
  out[7] += 3.464101615137754*f[5]*w0dx0+0.8783100656536796*f[9]*dv0dx0+0.8944271909999159*f[2]*dv0dx0; 
  out[8] += 11.83215956619923*f[4]*w0dx0+5.291502622129181*f[0]*w0dx0+3.415650255319867*f[6]*dv0dx0+1.527525231651947*f[2]*dv0dx0; 
  out[10] += 7.745966692414834*f[7]*w0dx0+1.963961012123931*f[12]*dv0dx0+2.0*f[3]*dv0dx0; 
  out[11] += 11.83215956619923*f[6]*w0dx0+5.291502622129181*f[2]*w0dx0+3.055050463303894*f[10]*dv0dx0+1.366260102127947*f[5]*dv0dx0+3.415650255319867*f[4]*dv0dx0+1.527525231651947*f[0]*dv0dx0; 
  out[12] += 3.464101615137754*f[9]*w0dx0+0.8728715609439696*f[14]*dv0dx0+0.8783100656536796*f[5]*dv0dx0; 
  out[13] += 15.87450786638754*f[8]*w0dx0+10.39230484541326*f[1]*w0dx0+4.58257569495584*f[11]*dv0dx0+3.0*f[3]*dv0dx0; 
} 