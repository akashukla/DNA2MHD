-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in calculating charge exchange vRel and reaction rate.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[

void GkProdCXcellAvSer1x1v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkProdCXcellAvSer1x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer2x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer3x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer1x1v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkProdCXcellAvSer1x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer2x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer3x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer1x1v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkProdCXcellAvSer1x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer2x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvSer3x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax1x1v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkProdCXcellAvMax1x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax2x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax3x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax1x1v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkProdCXcellAvMax1x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax2x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax3x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax1x1v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 

void GkProdCXcellAvMax1x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax2x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 


void GkProdCXcellAvMax3x2v_P3(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX); 



]]
