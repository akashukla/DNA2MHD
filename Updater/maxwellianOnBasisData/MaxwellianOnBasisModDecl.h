#ifndef MAXWELLIANONBASIS_MOD_DECL_H 
#define MAXWELLIANONBASIS_MOD_DECL_H 
#include <cmath> 

extern "C" { 

void MaxwellianOnBasisGauss1x1vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x1vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x1vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x1vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x1vSer_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x1vSer_P3_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss1x1vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x1vSerUz_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x1vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss1x1vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x1vSerUz_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x1vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x2vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vSerUpar_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x2vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vSerUpar_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x2vSer_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vSerUpar_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vSer_P3_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss1x2vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x2vSerUz_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x2vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss1x2vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x2vSerUz_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x2vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x3vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x3vSerUpar_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x3vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x3vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x3vSerUpar_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x3vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x3vSer_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x3vSerUpar_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x3vSer_P3_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss2x2vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vSerUpar_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss2x2vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vSerUpar_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss2x2vSer_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vSerUpar_P3_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vSer_P3_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss2x2vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss2x2vSerUz_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss2x2vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss2x2vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss2x2vSerUz_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss2x2vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss2x3vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x3vSerUpar_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x3vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss2x3vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x3vSerUpar_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x3vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss3x3vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss3x3vSerUpar_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss3x3vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss3x3vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss3x3vSerUpar_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss3x3vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss3x2vSer_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss3x2vSerUz_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss3x2vSer_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss3x2vSer_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss3x2vSerUz_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss3x2vSer_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x1vMax_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x1vMax_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x1vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x1vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss1x1vMax_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x1vMaxUz_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x1vMax_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss1x1vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x1vMaxUz_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x1vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x2vMax_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vMaxUpar_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vMax_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x2vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vMaxUpar_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss1x2vMax_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x2vMaxUz_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x2vMax_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss1x2vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x2vMaxUz_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x2vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss2x2vMax_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vMaxUpar_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vMax_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss2x2vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vMaxUpar_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss2x2vMax_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss2x2vMaxUz_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss2x2vMax_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss2x2vMax_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss2x2vMaxUz_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss2x2vMax_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x1vTensor_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x1vTensor_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x1vTensor_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x1vTensor_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss1x1vTensor_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x1vTensorUz_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x1vTensor_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss1x1vTensor_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x1vTensorUz_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x1vTensor_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x2vTensor_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vTensorUpar_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vTensor_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss1x2vTensor_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vTensorUpar_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss1x2vTensor_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss1x2vTensor_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x2vTensorUz_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x2vTensor_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss1x2vTensor_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x2vTensorUz_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss1x2vTensor_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss2x2vTensor_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vTensorUpar_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vTensor_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void MaxwellianOnBasisGauss2x2vTensor_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vTensorUpar_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void MaxwellianOnBasisGauss2x2vTensor_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss2x2vTensor_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss2x2vTensorUz_P1_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss2x2vTensor_P1_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);
void GkMaxwellianOnBasisGauss2x2vTensor_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss2x2vTensorUz_P2_evAtConfOrd(const double *den, const double *flowU, const double *vtSq, const double *bmag, double *flowUOrd, double *vtSqOrd, double *fMFacOrd, double *bmagOrd);
void GkMaxwellianOnBasisGauss2x2vTensor_P2_phaseQuad(const double *flowUOrd, const double *vtSqOrd, const double *fMFacOrd, const double *bmagOrd, const double m_, const double *wc, const double *dxv, double *fMOut);

void MaxwellianInnerLoop(double *n, double *u, double *vtSq, double *bmag, double m_, double *fItr, double *weights, double *dz, double *zc, double *ordinates, double *basisAtOrdinates, double *phaseToConfOrdMap, int numPhaseBasis, int numConfOrds, int numPhaseOrds, int numConfDims, int numPhaseDims);
void GkMaxwellianInnerLoop(double *n, double *u, double *vtSq, double *bmag, double m_, double *fItr, double *weights, double *dz, double *zc, double *ordinates, double *basisAtOrdinates, double *phaseToConfOrdMap, int numPhaseBasis, int numConfOrds, int numPhaseOrds, int numConfDims, int numPhaseDims);
} 

#endif 
