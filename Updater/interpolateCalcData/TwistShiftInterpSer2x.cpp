#include <TwistShiftInterpModDecl.h> 
 
void TwistShiftInterp2xSer_P1(const double *xLimLo, const double *xLimUp, const double *fldSrc, double *fldDest) 
{ 
  // xLimUp:  1D DG expansion of the function xLimLo(y) giving the lower limit of the x integral.
  // xLimLo:  1D DG expansion of the function xLimUp(y) giving the upper limit of the x integral.
  // wSrc:    cell center of source grid cell.
  // wDest:   cell center of destination grid cell.
  // dxSrc:   cell length of source grid cell.
  // dxDest:  cell length of destination grid cell.
  // fldSrc:  source grid field.
  // fldDest: destination grid field.

  double eLo[2];
  double eUp[2];
  eLo[1] = std::max(-1.0,-1.0);
  eUp[1] = std::min( 1.0,1.0);

  double eLo1R2 = std::pow(eLo[1],2);
  double eLo1R3 = std::pow(eLo[1],3);
  double eLo1R4 = std::pow(eLo[1],4);
  double eLo1R5 = std::pow(eLo[1],5);
  double eLo1R6 = std::pow(eLo[1],6);
  double eUp1R2 = std::pow(eUp[1],2);
  double eUp1R3 = std::pow(eUp[1],3);
  double eUp1R4 = std::pow(eUp[1],4);
  double eUp1R5 = std::pow(eUp[1],5);
  double eUp1R6 = std::pow(eUp[1],6);
  double xLimLo0R2 = std::pow(xLimLo[0],2);
  double xLimLo1R2 = std::pow(xLimLo[1],2);
  double xLimLo0R3 = std::pow(xLimLo[0],3);
  double xLimLo1R3 = std::pow(xLimLo[1],3);
  double xLimUp0R2 = std::pow(xLimUp[0],2);
  double xLimUp1R2 = std::pow(xLimUp[1],2);
  double xLimUp0R3 = std::pow(xLimUp[0],3);
  double xLimUp1R3 = std::pow(xLimUp[1],3);

  fldDest[0] += 0.015625*(((9.0*eUp1R4-9.0*eLo1R4)*xLimUp1R2+(13.85640646055102*xLimUp[0]*eUp1R3-13.85640646055102*xLimUp[0]*eLo1R3)*xLimUp[1]+(9.0*eLo1R4-9.0*eUp1R4)*xLimLo1R2+(13.85640646055102*xLimLo[0]*eLo1R3-13.85640646055102*xLimLo[0]*eUp1R3)*xLimLo[1]+(6.0*xLimUp0R2-6.0*xLimLo0R2)*eUp1R2+(6.0*xLimLo0R2-6.0*xLimUp0R2)*eLo1R2)*fldSrc[3]+((11.31370849898477*eUp1R3-11.31370849898477*eLo1R3)*xLimUp[1]+(11.31370849898477*eLo1R3-11.31370849898477*eUp1R3)*xLimLo[1]+(9.797958971132715*xLimUp[0]-9.797958971132715*xLimLo[0])*eUp1R2+(9.797958971132715*xLimLo[0]-9.797958971132715*xLimUp[0])*eLo1R2)*fldSrc[2]+(6.928203230275509*eUp1R3-6.928203230275509*eLo1R3)*fldSrc[1]*xLimUp1R2+((12.0*xLimUp[0]*eUp1R2-12.0*xLimUp[0]*eLo1R2)*fldSrc[1]+9.797958971132715*fldSrc[0]*eUp1R2-9.797958971132715*fldSrc[0]*eLo1R2)*xLimUp[1]+(6.928203230275509*eLo1R3-6.928203230275509*eUp1R3)*fldSrc[1]*xLimLo1R2+((12.0*xLimLo[0]*eLo1R2-12.0*xLimLo[0]*eUp1R2)*fldSrc[1]-9.797958971132715*fldSrc[0]*eUp1R2+9.797958971132715*fldSrc[0]*eLo1R2)*xLimLo[1]+((6.928203230275509*xLimUp0R2-6.928203230275509*xLimLo0R2)*eUp[1]+(6.928203230275509*xLimLo0R2-6.928203230275509*xLimUp0R2)*eLo[1])*fldSrc[1]+(11.31370849898477*fldSrc[0]*xLimUp[0]-11.31370849898477*fldSrc[0]*xLimLo[0])*eUp[1]+(11.31370849898477*fldSrc[0]*xLimLo[0]-11.31370849898477*fldSrc[0]*xLimUp[0])*eLo[1]); 
  fldDest[1] += 0.003125*(((50.91168824543144*eUp1R5-50.91168824543144*eLo1R5)*xLimUp1R3+(110.227038425243*xLimUp[0]*eUp1R4-110.227038425243*xLimUp[0]*eLo1R4)*xLimUp1R2+(84.85281374238573*xLimUp0R2*eUp1R3-84.85281374238573*xLimUp0R2*eLo1R3)*xLimUp[1]+(50.91168824543144*eLo1R5-50.91168824543144*eUp1R5)*xLimLo1R3+(110.227038425243*xLimLo[0]*eLo1R4-110.227038425243*xLimLo[0]*eUp1R4)*xLimLo1R2+(84.85281374238573*xLimLo0R2*eLo1R3-84.85281374238573*xLimLo0R2*eUp1R3)*xLimLo[1]+(24.49489742783179*xLimUp0R3-24.49489742783179*xLimLo0R3)*eUp1R2+(24.49489742783179*xLimLo0R3-24.49489742783179*xLimUp0R3)*eLo1R2)*fldSrc[3]+((45.0*eUp1R4-45.0*eLo1R4)*xLimUp1R2+(69.28203230275508*xLimUp[0]*eUp1R3-69.28203230275508*xLimUp[0]*eLo1R3)*xLimUp[1]+(45.0*eLo1R4-45.0*eUp1R4)*xLimLo1R2+(69.28203230275508*xLimLo[0]*eLo1R3-69.28203230275508*xLimLo[0]*eUp1R3)*xLimLo[1]+(30.0*xLimUp0R2-30.0*xLimLo0R2)*eUp1R2+(30.0*xLimLo0R2-30.0*xLimUp0R2)*eLo1R2)*fldSrc[2]+(36.74234614174767*eUp1R4-36.74234614174767*eLo1R4)*fldSrc[1]*xLimUp1R3+((84.85281374238573*xLimUp[0]*eUp1R3-84.85281374238573*xLimUp[0]*eLo1R3)*fldSrc[1]+34.64101615137754*fldSrc[0]*eUp1R3-34.64101615137754*fldSrc[0]*eLo1R3)*xLimUp1R2+((73.48469228349535*xLimUp0R2*eUp1R2-73.48469228349535*xLimUp0R2*eLo1R2)*fldSrc[1]+60.0*fldSrc[0]*xLimUp[0]*eUp1R2-60.0*fldSrc[0]*xLimUp[0]*eLo1R2)*xLimUp[1]+(36.74234614174767*eLo1R4-36.74234614174767*eUp1R4)*fldSrc[1]*xLimLo1R3+((84.85281374238573*xLimLo[0]*eLo1R3-84.85281374238573*xLimLo[0]*eUp1R3)*fldSrc[1]-34.64101615137754*fldSrc[0]*eUp1R3+34.64101615137754*fldSrc[0]*eLo1R3)*xLimLo1R2+((73.48469228349535*xLimLo0R2*eLo1R2-73.48469228349535*xLimLo0R2*eUp1R2)*fldSrc[1]-60.0*fldSrc[0]*xLimLo[0]*eUp1R2+60.0*fldSrc[0]*xLimLo[0]*eLo1R2)*xLimLo[1]+((28.28427124746191*xLimUp0R3-28.28427124746191*xLimLo0R3)*eUp[1]+(28.28427124746191*xLimLo0R3-28.28427124746191*xLimUp0R3)*eLo[1])*fldSrc[1]+(34.64101615137754*fldSrc[0]*xLimUp0R2-34.64101615137754*fldSrc[0]*xLimLo0R2)*eUp[1]+(34.64101615137754*fldSrc[0]*xLimLo0R2-34.64101615137754*fldSrc[0]*xLimUp0R2)*eLo[1]); 
  fldDest[2] += 0.003125*(((62.35382907247956*eUp1R5-62.35382907247956*eLo1R5)*xLimUp1R2+(90.0*xLimUp[0]*eUp1R4-90.0*xLimUp[0]*eLo1R4)*xLimUp[1]+(62.35382907247956*eLo1R5-62.35382907247956*eUp1R5)*xLimLo1R2+(90.0*xLimLo[0]*eLo1R4-90.0*xLimLo[0]*eUp1R4)*xLimLo[1]+(34.64101615137754*xLimUp0R2-34.64101615137754*xLimLo0R2)*eUp1R3+(34.64101615137754*xLimLo0R2-34.64101615137754*xLimUp0R2)*eLo1R3)*fldSrc[3]+((73.48469228349535*eUp1R4-73.48469228349535*eLo1R4)*xLimUp[1]+(73.48469228349535*eLo1R4-73.48469228349535*eUp1R4)*xLimLo[1]+(56.56854249492383*xLimUp[0]-56.56854249492383*xLimLo[0])*eUp1R3+(56.56854249492383*xLimLo[0]-56.56854249492383*xLimUp[0])*eLo1R3)*fldSrc[2]+(45.0*eUp1R4-45.0*eLo1R4)*fldSrc[1]*xLimUp1R2+((69.28203230275508*xLimUp[0]*eUp1R3-69.28203230275508*xLimUp[0]*eLo1R3)*fldSrc[1]+56.56854249492383*fldSrc[0]*eUp1R3-56.56854249492383*fldSrc[0]*eLo1R3)*xLimUp[1]+(45.0*eLo1R4-45.0*eUp1R4)*fldSrc[1]*xLimLo1R2+((69.28203230275508*xLimLo[0]*eLo1R3-69.28203230275508*xLimLo[0]*eUp1R3)*fldSrc[1]-56.56854249492383*fldSrc[0]*eUp1R3+56.56854249492383*fldSrc[0]*eLo1R3)*xLimLo[1]+((30.0*xLimUp0R2-30.0*xLimLo0R2)*eUp1R2+(30.0*xLimLo0R2-30.0*xLimUp0R2)*eLo1R2)*fldSrc[1]+(48.98979485566358*fldSrc[0]*xLimUp[0]-48.98979485566358*fldSrc[0]*xLimLo[0])*eUp1R2+(48.98979485566358*fldSrc[0]*xLimLo[0]-48.98979485566358*fldSrc[0]*xLimUp[0])*eLo1R2); 
  fldDest[3] += 0.003125*(((73.48469228349535*eUp1R6-73.48469228349535*eLo1R6)*xLimUp1R3+(152.7350647362943*xLimUp[0]*eUp1R5-152.7350647362943*xLimUp[0]*eLo1R5)*xLimUp1R2+(110.227038425243*xLimUp0R2*eUp1R4-110.227038425243*xLimUp0R2*eLo1R4)*xLimUp[1]+(73.48469228349535*eLo1R6-73.48469228349535*eUp1R6)*xLimLo1R3+(152.7350647362943*xLimLo[0]*eLo1R5-152.7350647362943*xLimLo[0]*eUp1R5)*xLimLo1R2+(110.227038425243*xLimLo0R2*eLo1R4-110.227038425243*xLimLo0R2*eUp1R4)*xLimLo[1]+(28.28427124746191*xLimUp0R3-28.28427124746191*xLimLo0R3)*eUp1R3+(28.28427124746191*xLimLo0R3-28.28427124746191*xLimUp0R3)*eLo1R3)*fldSrc[3]+((62.35382907247956*eUp1R5-62.35382907247956*eLo1R5)*xLimUp1R2+(90.0*xLimUp[0]*eUp1R4-90.0*xLimUp[0]*eLo1R4)*xLimUp[1]+(62.35382907247956*eLo1R5-62.35382907247956*eUp1R5)*xLimLo1R2+(90.0*xLimLo[0]*eLo1R4-90.0*xLimLo[0]*eUp1R4)*xLimLo[1]+(34.64101615137754*xLimUp0R2-34.64101615137754*xLimLo0R2)*eUp1R3+(34.64101615137754*xLimLo0R2-34.64101615137754*xLimUp0R2)*eLo1R3)*fldSrc[2]+(50.91168824543144*eUp1R5-50.91168824543144*eLo1R5)*fldSrc[1]*xLimUp1R3+((110.227038425243*xLimUp[0]*eUp1R4-110.227038425243*xLimUp[0]*eLo1R4)*fldSrc[1]+45.0*fldSrc[0]*eUp1R4-45.0*fldSrc[0]*eLo1R4)*xLimUp1R2+((84.85281374238573*xLimUp0R2*eUp1R3-84.85281374238573*xLimUp0R2*eLo1R3)*fldSrc[1]+69.28203230275508*fldSrc[0]*xLimUp[0]*eUp1R3-69.28203230275508*fldSrc[0]*xLimUp[0]*eLo1R3)*xLimUp[1]+(50.91168824543144*eLo1R5-50.91168824543144*eUp1R5)*fldSrc[1]*xLimLo1R3+((110.227038425243*xLimLo[0]*eLo1R4-110.227038425243*xLimLo[0]*eUp1R4)*fldSrc[1]-45.0*fldSrc[0]*eUp1R4+45.0*fldSrc[0]*eLo1R4)*xLimLo1R2+((84.85281374238573*xLimLo0R2*eLo1R3-84.85281374238573*xLimLo0R2*eUp1R3)*fldSrc[1]-69.28203230275508*fldSrc[0]*xLimLo[0]*eUp1R3+69.28203230275508*fldSrc[0]*xLimLo[0]*eLo1R3)*xLimLo[1]+((24.49489742783179*xLimUp0R3-24.49489742783179*xLimLo0R3)*eUp1R2+(24.49489742783179*xLimLo0R3-24.49489742783179*xLimUp0R3)*eLo1R2)*fldSrc[1]+(30.0*fldSrc[0]*xLimUp0R2-30.0*fldSrc[0]*xLimLo0R2)*eUp1R2+(30.0*fldSrc[0]*xLimLo0R2-30.0*fldSrc[0]*xLimUp0R2)*eLo1R2); 

}

void TwistShiftInterp_limXvarYfixed2xSer_P1(const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double *fldSrc, double *fldDest) 
{ 
  // xLimUp:  1D DG expansion of the function xLimLo(y) giving the lower limit of the x integral.
  // xLimLo:  1D DG expansion of the function xLimUp(y) giving the upper limit of the x integral.
  // yLimUp:  lower limit of the y integral in [-1,1] normalized space.
  // yLimLo:  upper limit of the y integral in [-1,1] normalized space.
  // fldSrc:  source grid field.
  // fldDest: destination grid field.

  double yLimLoR2 = std::pow(yLimLo,2);
  double yLimLoR3 = std::pow(yLimLo,3);
  double yLimLoR4 = std::pow(yLimLo,4);
  double yLimLoR5 = std::pow(yLimLo,5);
  double yLimLoR6 = std::pow(yLimLo,6);
  double yLimUpR2 = std::pow(yLimUp,2);
  double yLimUpR3 = std::pow(yLimUp,3);
  double yLimUpR4 = std::pow(yLimUp,4);
  double yLimUpR5 = std::pow(yLimUp,5);
  double yLimUpR6 = std::pow(yLimUp,6);
  double xLimLo0R2 = std::pow(xLimLo[0],2);
  double xLimLo1R2 = std::pow(xLimLo[1],2);
  double xLimLo0R3 = std::pow(xLimLo[0],3);
  double xLimLo1R3 = std::pow(xLimLo[1],3);
  double xLimUp0R2 = std::pow(xLimUp[0],2);
  double xLimUp1R2 = std::pow(xLimUp[1],2);
  double xLimUp0R3 = std::pow(xLimUp[0],3);
  double xLimUp1R3 = std::pow(xLimUp[1],3);

  fldDest[0] += 0.015625*((9.0*xLimUp1R2-9.0*xLimLo1R2)*fldSrc[3]*yLimUpR4+((13.85640646055102*xLimUp[0]*xLimUp[1]-13.85640646055102*xLimLo[0]*xLimLo[1])*fldSrc[3]+(11.31370849898477*xLimUp[1]-11.31370849898477*xLimLo[1])*fldSrc[2]+6.928203230275509*fldSrc[1]*xLimUp1R2-6.928203230275509*fldSrc[1]*xLimLo1R2)*yLimUpR3+((6.0*xLimUp0R2-6.0*xLimLo0R2)*fldSrc[3]+(9.797958971132715*xLimUp[0]-9.797958971132715*xLimLo[0])*fldSrc[2]+(12.0*xLimUp[0]*fldSrc[1]+9.797958971132715*fldSrc[0])*xLimUp[1]+((-12.0*xLimLo[0]*fldSrc[1])-9.797958971132715*fldSrc[0])*xLimLo[1])*yLimUpR2+((6.928203230275509*xLimUp0R2-6.928203230275509*xLimLo0R2)*fldSrc[1]+11.31370849898477*fldSrc[0]*xLimUp[0]-11.31370849898477*fldSrc[0]*xLimLo[0])*yLimUp+(9.0*xLimLo1R2-9.0*xLimUp1R2)*fldSrc[3]*yLimLoR4+((13.85640646055102*xLimLo[0]*xLimLo[1]-13.85640646055102*xLimUp[0]*xLimUp[1])*fldSrc[3]+(11.31370849898477*xLimLo[1]-11.31370849898477*xLimUp[1])*fldSrc[2]-6.928203230275509*fldSrc[1]*xLimUp1R2+6.928203230275509*fldSrc[1]*xLimLo1R2)*yLimLoR3+((6.0*xLimLo0R2-6.0*xLimUp0R2)*fldSrc[3]+(9.797958971132715*xLimLo[0]-9.797958971132715*xLimUp[0])*fldSrc[2]+((-12.0*xLimUp[0]*fldSrc[1])-9.797958971132715*fldSrc[0])*xLimUp[1]+(12.0*xLimLo[0]*fldSrc[1]+9.797958971132715*fldSrc[0])*xLimLo[1])*yLimLoR2+((6.928203230275509*xLimLo0R2-6.928203230275509*xLimUp0R2)*fldSrc[1]-11.31370849898477*fldSrc[0]*xLimUp[0]+11.31370849898477*fldSrc[0]*xLimLo[0])*yLimLo); 
  fldDest[1] += 0.003125*((50.91168824543144*xLimUp1R3-50.91168824543144*xLimLo1R3)*fldSrc[3]*yLimUpR5+((110.227038425243*xLimUp[0]*xLimUp1R2-110.227038425243*xLimLo[0]*xLimLo1R2)*fldSrc[3]+(45.0*xLimUp1R2-45.0*xLimLo1R2)*fldSrc[2]+36.74234614174767*fldSrc[1]*xLimUp1R3-36.74234614174767*fldSrc[1]*xLimLo1R3)*yLimUpR4+((84.85281374238573*xLimUp0R2*xLimUp[1]-84.85281374238573*xLimLo0R2*xLimLo[1])*fldSrc[3]+(69.28203230275508*xLimUp[0]*xLimUp[1]-69.28203230275508*xLimLo[0]*xLimLo[1])*fldSrc[2]+(84.85281374238573*xLimUp[0]*fldSrc[1]+34.64101615137754*fldSrc[0])*xLimUp1R2+((-84.85281374238573*xLimLo[0]*fldSrc[1])-34.64101615137754*fldSrc[0])*xLimLo1R2)*yLimUpR3+((24.49489742783179*xLimUp0R3-24.49489742783179*xLimLo0R3)*fldSrc[3]+(30.0*xLimUp0R2-30.0*xLimLo0R2)*fldSrc[2]+(73.48469228349535*xLimUp0R2*fldSrc[1]+60.0*fldSrc[0]*xLimUp[0])*xLimUp[1]+((-73.48469228349535*xLimLo0R2*fldSrc[1])-60.0*fldSrc[0]*xLimLo[0])*xLimLo[1])*yLimUpR2+((28.28427124746191*xLimUp0R3-28.28427124746191*xLimLo0R3)*fldSrc[1]+34.64101615137754*fldSrc[0]*xLimUp0R2-34.64101615137754*fldSrc[0]*xLimLo0R2)*yLimUp+(50.91168824543144*xLimLo1R3-50.91168824543144*xLimUp1R3)*fldSrc[3]*yLimLoR5+((110.227038425243*xLimLo[0]*xLimLo1R2-110.227038425243*xLimUp[0]*xLimUp1R2)*fldSrc[3]+(45.0*xLimLo1R2-45.0*xLimUp1R2)*fldSrc[2]-36.74234614174767*fldSrc[1]*xLimUp1R3+36.74234614174767*fldSrc[1]*xLimLo1R3)*yLimLoR4+((84.85281374238573*xLimLo0R2*xLimLo[1]-84.85281374238573*xLimUp0R2*xLimUp[1])*fldSrc[3]+(69.28203230275508*xLimLo[0]*xLimLo[1]-69.28203230275508*xLimUp[0]*xLimUp[1])*fldSrc[2]+((-84.85281374238573*xLimUp[0]*fldSrc[1])-34.64101615137754*fldSrc[0])*xLimUp1R2+(84.85281374238573*xLimLo[0]*fldSrc[1]+34.64101615137754*fldSrc[0])*xLimLo1R2)*yLimLoR3+((24.49489742783179*xLimLo0R3-24.49489742783179*xLimUp0R3)*fldSrc[3]+(30.0*xLimLo0R2-30.0*xLimUp0R2)*fldSrc[2]+((-73.48469228349535*xLimUp0R2*fldSrc[1])-60.0*fldSrc[0]*xLimUp[0])*xLimUp[1]+(73.48469228349535*xLimLo0R2*fldSrc[1]+60.0*fldSrc[0]*xLimLo[0])*xLimLo[1])*yLimLoR2+((28.28427124746191*xLimLo0R3-28.28427124746191*xLimUp0R3)*fldSrc[1]-34.64101615137754*fldSrc[0]*xLimUp0R2+34.64101615137754*fldSrc[0]*xLimLo0R2)*yLimLo); 
  fldDest[2] += 0.003125*((62.35382907247956*xLimUp1R2-62.35382907247956*xLimLo1R2)*fldSrc[3]*yLimUpR5+((90.0*xLimUp[0]*xLimUp[1]-90.0*xLimLo[0]*xLimLo[1])*fldSrc[3]+(73.48469228349535*xLimUp[1]-73.48469228349535*xLimLo[1])*fldSrc[2]+45.0*fldSrc[1]*xLimUp1R2-45.0*fldSrc[1]*xLimLo1R2)*yLimUpR4+((34.64101615137754*xLimUp0R2-34.64101615137754*xLimLo0R2)*fldSrc[3]+(56.56854249492383*xLimUp[0]-56.56854249492383*xLimLo[0])*fldSrc[2]+(69.28203230275508*xLimUp[0]*fldSrc[1]+56.56854249492383*fldSrc[0])*xLimUp[1]+((-69.28203230275508*xLimLo[0]*fldSrc[1])-56.56854249492383*fldSrc[0])*xLimLo[1])*yLimUpR3+((30.0*xLimUp0R2-30.0*xLimLo0R2)*fldSrc[1]+48.98979485566358*fldSrc[0]*xLimUp[0]-48.98979485566358*fldSrc[0]*xLimLo[0])*yLimUpR2+(62.35382907247956*xLimLo1R2-62.35382907247956*xLimUp1R2)*fldSrc[3]*yLimLoR5+((90.0*xLimLo[0]*xLimLo[1]-90.0*xLimUp[0]*xLimUp[1])*fldSrc[3]+(73.48469228349535*xLimLo[1]-73.48469228349535*xLimUp[1])*fldSrc[2]-45.0*fldSrc[1]*xLimUp1R2+45.0*fldSrc[1]*xLimLo1R2)*yLimLoR4+((34.64101615137754*xLimLo0R2-34.64101615137754*xLimUp0R2)*fldSrc[3]+(56.56854249492383*xLimLo[0]-56.56854249492383*xLimUp[0])*fldSrc[2]+((-69.28203230275508*xLimUp[0]*fldSrc[1])-56.56854249492383*fldSrc[0])*xLimUp[1]+(69.28203230275508*xLimLo[0]*fldSrc[1]+56.56854249492383*fldSrc[0])*xLimLo[1])*yLimLoR3+((30.0*xLimLo0R2-30.0*xLimUp0R2)*fldSrc[1]-48.98979485566358*fldSrc[0]*xLimUp[0]+48.98979485566358*fldSrc[0]*xLimLo[0])*yLimLoR2); 
  fldDest[3] += 0.003125*((73.48469228349535*xLimUp1R3-73.48469228349535*xLimLo1R3)*fldSrc[3]*yLimUpR6+((152.7350647362943*xLimUp[0]*xLimUp1R2-152.7350647362943*xLimLo[0]*xLimLo1R2)*fldSrc[3]+(62.35382907247956*xLimUp1R2-62.35382907247956*xLimLo1R2)*fldSrc[2]+50.91168824543144*fldSrc[1]*xLimUp1R3-50.91168824543144*fldSrc[1]*xLimLo1R3)*yLimUpR5+((110.227038425243*xLimUp0R2*xLimUp[1]-110.227038425243*xLimLo0R2*xLimLo[1])*fldSrc[3]+(90.0*xLimUp[0]*xLimUp[1]-90.0*xLimLo[0]*xLimLo[1])*fldSrc[2]+(110.227038425243*xLimUp[0]*fldSrc[1]+45.0*fldSrc[0])*xLimUp1R2+((-110.227038425243*xLimLo[0]*fldSrc[1])-45.0*fldSrc[0])*xLimLo1R2)*yLimUpR4+((28.28427124746191*xLimUp0R3-28.28427124746191*xLimLo0R3)*fldSrc[3]+(34.64101615137754*xLimUp0R2-34.64101615137754*xLimLo0R2)*fldSrc[2]+(84.85281374238573*xLimUp0R2*fldSrc[1]+69.28203230275508*fldSrc[0]*xLimUp[0])*xLimUp[1]+((-84.85281374238573*xLimLo0R2*fldSrc[1])-69.28203230275508*fldSrc[0]*xLimLo[0])*xLimLo[1])*yLimUpR3+((24.49489742783179*xLimUp0R3-24.49489742783179*xLimLo0R3)*fldSrc[1]+30.0*fldSrc[0]*xLimUp0R2-30.0*fldSrc[0]*xLimLo0R2)*yLimUpR2+(73.48469228349535*xLimLo1R3-73.48469228349535*xLimUp1R3)*fldSrc[3]*yLimLoR6+((152.7350647362943*xLimLo[0]*xLimLo1R2-152.7350647362943*xLimUp[0]*xLimUp1R2)*fldSrc[3]+(62.35382907247956*xLimLo1R2-62.35382907247956*xLimUp1R2)*fldSrc[2]-50.91168824543144*fldSrc[1]*xLimUp1R3+50.91168824543144*fldSrc[1]*xLimLo1R3)*yLimLoR5+((110.227038425243*xLimLo0R2*xLimLo[1]-110.227038425243*xLimUp0R2*xLimUp[1])*fldSrc[3]+(90.0*xLimLo[0]*xLimLo[1]-90.0*xLimUp[0]*xLimUp[1])*fldSrc[2]+((-110.227038425243*xLimUp[0]*fldSrc[1])-45.0*fldSrc[0])*xLimUp1R2+(110.227038425243*xLimLo[0]*fldSrc[1]+45.0*fldSrc[0])*xLimLo1R2)*yLimLoR4+((28.28427124746191*xLimLo0R3-28.28427124746191*xLimUp0R3)*fldSrc[3]+(34.64101615137754*xLimLo0R2-34.64101615137754*xLimUp0R2)*fldSrc[2]+((-84.85281374238573*xLimUp0R2*fldSrc[1])-69.28203230275508*fldSrc[0]*xLimUp[0])*xLimUp[1]+(84.85281374238573*xLimLo0R2*fldSrc[1]+69.28203230275508*fldSrc[0]*xLimLo[0])*xLimLo[1])*yLimLoR3+((24.49489742783179*xLimLo0R3-24.49489742783179*xLimUp0R3)*fldSrc[1]-30.0*fldSrc[0]*xLimUp0R2+30.0*fldSrc[0]*xLimLo0R2)*yLimLoR2); 

}

void TwistShiftInterp_xLimDG2xSer_P1(const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyLim, const double ycLim, const double *fldSrc, double *fldDest) 
{ 
  // xLimUp:  1D DG expansion of the function xLimLo(y) giving the lower limit of the x integral.
  // xLimLo:  1D DG expansion of the function xLimUp(y) giving the upper limit of the x integral.
  // yLimUp:  lower limit of the y integral in [-1,1] normalized space.
  // yLimLo:  upper limit of the y integral in [-1,1] normalized space.
  // dyLimSpace:  length of the subregion in which the DG expansion of the yLimLo/yLimUp functions are defined.
  // ycLimSpace:  logical center of the subregion in which the DG expansion of the yLimLo/yLimUp functions are defined.
  // fldDest: destination grid field.

  double yLimLoR2 = std::pow(yLimLo,2);
  double yLimLoR3 = std::pow(yLimLo,3);
  double yLimLoR4 = std::pow(yLimLo,4);
  double yLimLoR5 = std::pow(yLimLo,5);
  double yLimLoR6 = std::pow(yLimLo,6);
  double yLimUpR2 = std::pow(yLimUp,2);
  double yLimUpR3 = std::pow(yLimUp,3);
  double yLimUpR4 = std::pow(yLimUp,4);
  double yLimUpR5 = std::pow(yLimUp,5);
  double yLimUpR6 = std::pow(yLimUp,6);
  double dyLimR2 = std::pow(dyLim,2);
  double dyLimR3 = std::pow(dyLim,3);
  double ycLimR2 = std::pow(ycLim,2);
  double ycLimR3 = std::pow(ycLim,3);
  double xLimLo0R2 = std::pow(xLimLo[0],2);
  double xLimLo1R2 = std::pow(xLimLo[1],2);
  double xLimLo0R3 = std::pow(xLimLo[0],3);
  double xLimLo1R3 = std::pow(xLimLo[1],3);
  double xLimUp0R2 = std::pow(xLimUp[0],2);
  double xLimUp1R2 = std::pow(xLimUp[1],2);
  double xLimUp0R3 = std::pow(xLimUp[0],3);
  double xLimUp1R3 = std::pow(xLimUp[1],3);

  fldDest[0] += (0.03125*(((36.0*xLimUp1R2-36.0*xLimLo1R2)*fldSrc[3]*yLimUpR2+(41.56921938165305*fldSrc[1]*xLimUp1R2-41.56921938165305*fldSrc[1]*xLimLo1R2)*yLimUp+(36.0*xLimLo1R2-36.0*xLimUp1R2)*fldSrc[3]*yLimLoR2+(41.56921938165305*fldSrc[1]*xLimLo1R2-41.56921938165305*fldSrc[1]*xLimUp1R2)*yLimLo)*ycLimR2+((48.0*xLimLo1R2-48.0*xLimUp1R2)*fldSrc[3]*yLimUpR3+(((20.78460969082652*xLimLo[0]*xLimLo[1]-20.78460969082652*xLimUp[0]*xLimUp[1])*fldSrc[3]+(16.97056274847715*xLimLo[1]-16.97056274847715*xLimUp[1])*fldSrc[2])*dyLim-41.56921938165305*fldSrc[1]*xLimUp1R2+41.56921938165305*fldSrc[1]*xLimLo1R2)*yLimUpR2+(((-24.0*xLimUp[0]*fldSrc[1])-19.59591794226543*fldSrc[0])*xLimUp[1]+(24.0*xLimLo[0]*fldSrc[1]+19.59591794226543*fldSrc[0])*xLimLo[1])*dyLim*yLimUp+(48.0*xLimUp1R2-48.0*xLimLo1R2)*fldSrc[3]*yLimLoR3+(((20.78460969082652*xLimUp[0]*xLimUp[1]-20.78460969082652*xLimLo[0]*xLimLo[1])*fldSrc[3]+(16.97056274847715*xLimUp[1]-16.97056274847715*xLimLo[1])*fldSrc[2])*dyLim+41.56921938165305*fldSrc[1]*xLimUp1R2-41.56921938165305*fldSrc[1]*xLimLo1R2)*yLimLoR2+((24.0*xLimUp[0]*fldSrc[1]+19.59591794226543*fldSrc[0])*xLimUp[1]+((-24.0*xLimLo[0]*fldSrc[1])-19.59591794226543*fldSrc[0])*xLimLo[1])*dyLim*yLimLo)*ycLim+(18.0*xLimUp1R2-18.0*xLimLo1R2)*fldSrc[3]*yLimUpR4+(((13.85640646055102*xLimUp[0]*xLimUp[1]-13.85640646055102*xLimLo[0]*xLimLo[1])*fldSrc[3]+(11.31370849898477*xLimUp[1]-11.31370849898477*xLimLo[1])*fldSrc[2])*dyLim+13.85640646055102*fldSrc[1]*xLimUp1R2-13.85640646055102*fldSrc[1]*xLimLo1R2)*yLimUpR3+(((3.0*xLimUp0R2-3.0*xLimLo0R2)*fldSrc[3]+(4.898979485566357*xLimUp[0]-4.898979485566357*xLimLo[0])*fldSrc[2])*dyLimR2+((12.0*xLimUp[0]*fldSrc[1]+9.797958971132715*fldSrc[0])*xLimUp[1]+((-12.0*xLimLo[0]*fldSrc[1])-9.797958971132715*fldSrc[0])*xLimLo[1])*dyLim)*yLimUpR2+((3.464101615137754*xLimUp0R2-3.464101615137754*xLimLo0R2)*fldSrc[1]+5.656854249492382*fldSrc[0]*xLimUp[0]-5.656854249492382*fldSrc[0]*xLimLo[0])*dyLimR2*yLimUp+(18.0*xLimLo1R2-18.0*xLimUp1R2)*fldSrc[3]*yLimLoR4+(((13.85640646055102*xLimLo[0]*xLimLo[1]-13.85640646055102*xLimUp[0]*xLimUp[1])*fldSrc[3]+(11.31370849898477*xLimLo[1]-11.31370849898477*xLimUp[1])*fldSrc[2])*dyLim-13.85640646055102*fldSrc[1]*xLimUp1R2+13.85640646055102*fldSrc[1]*xLimLo1R2)*yLimLoR3+(((3.0*xLimLo0R2-3.0*xLimUp0R2)*fldSrc[3]+(4.898979485566357*xLimLo[0]-4.898979485566357*xLimUp[0])*fldSrc[2])*dyLimR2+(((-12.0*xLimUp[0]*fldSrc[1])-9.797958971132715*fldSrc[0])*xLimUp[1]+(12.0*xLimLo[0]*fldSrc[1]+9.797958971132715*fldSrc[0])*xLimLo[1])*dyLim)*yLimLoR2+((3.464101615137754*xLimLo0R2-3.464101615137754*xLimUp0R2)*fldSrc[1]-5.656854249492382*fldSrc[0]*xLimUp[0]+5.656854249492382*fldSrc[0]*xLimLo[0])*dyLimR2*yLimLo))/dyLimR2; 
  fldDest[1] += -(0.00625*(((509.1168824543145*xLimUp1R3-509.1168824543145*xLimLo1R3)*fldSrc[3]*yLimUpR2+(587.877538267963*fldSrc[1]*xLimUp1R3-587.877538267963*fldSrc[1]*xLimLo1R3)*yLimUp+(509.1168824543145*xLimLo1R3-509.1168824543145*xLimUp1R3)*fldSrc[3]*yLimLoR2+(587.877538267963*fldSrc[1]*xLimLo1R3-587.877538267963*fldSrc[1]*xLimUp1R3)*yLimLo)*ycLimR3+((1018.233764908629*xLimLo1R3-1018.233764908629*xLimUp1R3)*fldSrc[3]*yLimUpR3+(((440.9081537009721*xLimLo[0]*xLimLo1R2-440.9081537009721*xLimUp[0]*xLimUp1R2)*fldSrc[3]+(180.0*xLimLo1R2-180.0*xLimUp1R2)*fldSrc[2])*dyLim-881.8163074019443*fldSrc[1]*xLimUp1R3+881.8163074019443*fldSrc[1]*xLimLo1R3)*yLimUpR2+(((-509.1168824543145*xLimUp[0]*fldSrc[1])-207.8460969082653*fldSrc[0])*xLimUp1R2+(509.1168824543145*xLimLo[0]*fldSrc[1]+207.8460969082653*fldSrc[0])*xLimLo1R2)*dyLim*yLimUp+(1018.233764908629*xLimUp1R3-1018.233764908629*xLimLo1R3)*fldSrc[3]*yLimLoR3+(((440.9081537009721*xLimUp[0]*xLimUp1R2-440.9081537009721*xLimLo[0]*xLimLo1R2)*fldSrc[3]+(180.0*xLimUp1R2-180.0*xLimLo1R2)*fldSrc[2])*dyLim+881.8163074019443*fldSrc[1]*xLimUp1R3-881.8163074019443*fldSrc[1]*xLimLo1R3)*yLimLoR2+((509.1168824543145*xLimUp[0]*fldSrc[1]+207.8460969082653*fldSrc[0])*xLimUp1R2+((-509.1168824543145*xLimLo[0]*fldSrc[1])-207.8460969082653*fldSrc[0])*xLimLo1R2)*dyLim*yLimLo)*ycLimR2+((763.6753236814716*xLimUp1R3-763.6753236814716*xLimLo1R3)*fldSrc[3]*yLimUpR4+(((587.877538267963*xLimUp[0]*xLimUp1R2-587.877538267963*xLimLo[0]*xLimLo1R2)*fldSrc[3]+(240.0*xLimUp1R2-240.0*xLimLo1R2)*fldSrc[2])*dyLim+587.877538267963*fldSrc[1]*xLimUp1R3-587.877538267963*fldSrc[1]*xLimLo1R3)*yLimUpR3+(((127.2792206135786*xLimUp0R2*xLimUp[1]-127.2792206135786*xLimLo0R2*xLimLo[1])*fldSrc[3]+(103.9230484541326*xLimUp[0]*xLimUp[1]-103.9230484541326*xLimLo[0]*xLimLo[1])*fldSrc[2])*dyLimR2+((509.1168824543145*xLimUp[0]*fldSrc[1]+207.8460969082653*fldSrc[0])*xLimUp1R2+((-509.1168824543145*xLimLo[0]*fldSrc[1])-207.8460969082653*fldSrc[0])*xLimLo1R2)*dyLim)*yLimUpR2+((146.9693845669907*xLimUp0R2*fldSrc[1]+120.0*fldSrc[0]*xLimUp[0])*xLimUp[1]+((-146.9693845669907*xLimLo0R2*fldSrc[1])-120.0*fldSrc[0]*xLimLo[0])*xLimLo[1])*dyLimR2*yLimUp+(763.6753236814716*xLimLo1R3-763.6753236814716*xLimUp1R3)*fldSrc[3]*yLimLoR4+(((587.877538267963*xLimLo[0]*xLimLo1R2-587.877538267963*xLimUp[0]*xLimUp1R2)*fldSrc[3]+(240.0*xLimLo1R2-240.0*xLimUp1R2)*fldSrc[2])*dyLim-587.877538267963*fldSrc[1]*xLimUp1R3+587.877538267963*fldSrc[1]*xLimLo1R3)*yLimLoR3+(((127.2792206135786*xLimLo0R2*xLimLo[1]-127.2792206135786*xLimUp0R2*xLimUp[1])*fldSrc[3]+(103.9230484541326*xLimLo[0]*xLimLo[1]-103.9230484541326*xLimUp[0]*xLimUp[1])*fldSrc[2])*dyLimR2+(((-509.1168824543145*xLimUp[0]*fldSrc[1])-207.8460969082653*fldSrc[0])*xLimUp1R2+(509.1168824543145*xLimLo[0]*fldSrc[1]+207.8460969082653*fldSrc[0])*xLimLo1R2)*dyLim)*yLimLoR2+(((-146.9693845669907*xLimUp0R2*fldSrc[1])-120.0*fldSrc[0]*xLimUp[0])*xLimUp[1]+(146.9693845669907*xLimLo0R2*fldSrc[1]+120.0*fldSrc[0]*xLimLo[0])*xLimLo[1])*dyLimR2*yLimLo)*ycLim+(203.6467529817258*xLimLo1R3-203.6467529817258*xLimUp1R3)*fldSrc[3]*yLimUpR5+(((220.454076850486*xLimLo[0]*xLimLo1R2-220.454076850486*xLimUp[0]*xLimUp1R2)*fldSrc[3]+(90.0*xLimLo1R2-90.0*xLimUp1R2)*fldSrc[2])*dyLim-146.9693845669907*fldSrc[1]*xLimUp1R3+146.9693845669907*fldSrc[1]*xLimLo1R3)*yLimUpR4+(((84.85281374238573*xLimLo0R2*xLimLo[1]-84.85281374238573*xLimUp0R2*xLimUp[1])*fldSrc[3]+(69.28203230275508*xLimLo[0]*xLimLo[1]-69.28203230275508*xLimUp[0]*xLimUp[1])*fldSrc[2])*dyLimR2+(((-169.7056274847715*xLimUp[0]*fldSrc[1])-69.28203230275508*fldSrc[0])*xLimUp1R2+(169.7056274847715*xLimLo[0]*fldSrc[1]+69.28203230275508*fldSrc[0])*xLimLo1R2)*dyLim)*yLimUpR3+(((12.24744871391589*xLimLo0R3-12.24744871391589*xLimUp0R3)*fldSrc[3]+(15.0*xLimLo0R2-15.0*xLimUp0R2)*fldSrc[2])*dyLimR3+(((-73.48469228349535*xLimUp0R2*fldSrc[1])-60.0*fldSrc[0]*xLimUp[0])*xLimUp[1]+(73.48469228349535*xLimLo0R2*fldSrc[1]+60.0*fldSrc[0]*xLimLo[0])*xLimLo[1])*dyLimR2)*yLimUpR2+((14.14213562373095*xLimLo0R3-14.14213562373095*xLimUp0R3)*fldSrc[1]-17.32050807568877*fldSrc[0]*xLimUp0R2+17.32050807568877*fldSrc[0]*xLimLo0R2)*dyLimR3*yLimUp+(203.6467529817258*xLimUp1R3-203.6467529817258*xLimLo1R3)*fldSrc[3]*yLimLoR5+(((220.454076850486*xLimUp[0]*xLimUp1R2-220.454076850486*xLimLo[0]*xLimLo1R2)*fldSrc[3]+(90.0*xLimUp1R2-90.0*xLimLo1R2)*fldSrc[2])*dyLim+146.9693845669907*fldSrc[1]*xLimUp1R3-146.9693845669907*fldSrc[1]*xLimLo1R3)*yLimLoR4+(((84.85281374238573*xLimUp0R2*xLimUp[1]-84.85281374238573*xLimLo0R2*xLimLo[1])*fldSrc[3]+(69.28203230275508*xLimUp[0]*xLimUp[1]-69.28203230275508*xLimLo[0]*xLimLo[1])*fldSrc[2])*dyLimR2+((169.7056274847715*xLimUp[0]*fldSrc[1]+69.28203230275508*fldSrc[0])*xLimUp1R2+((-169.7056274847715*xLimLo[0]*fldSrc[1])-69.28203230275508*fldSrc[0])*xLimLo1R2)*dyLim)*yLimLoR3+(((12.24744871391589*xLimUp0R3-12.24744871391589*xLimLo0R3)*fldSrc[3]+(15.0*xLimUp0R2-15.0*xLimLo0R2)*fldSrc[2])*dyLimR3+((73.48469228349535*xLimUp0R2*fldSrc[1]+60.0*fldSrc[0]*xLimUp[0])*xLimUp[1]+((-73.48469228349535*xLimLo0R2*fldSrc[1])-60.0*fldSrc[0]*xLimLo[0])*xLimLo[1])*dyLimR2)*yLimLoR2+((14.14213562373095*xLimUp0R3-14.14213562373095*xLimLo0R3)*fldSrc[1]+17.32050807568877*fldSrc[0]*xLimUp0R2-17.32050807568877*fldSrc[0]*xLimLo0R2)*dyLimR3*yLimLo))/dyLimR3; 
  fldDest[2] += (0.00625*(((207.8460969082653*xLimUp1R2-207.8460969082653*xLimLo1R2)*fldSrc[3]*yLimUpR3+(180.0*fldSrc[1]*xLimUp1R2-180.0*fldSrc[1]*xLimLo1R2)*yLimUpR2+(207.8460969082653*xLimLo1R2-207.8460969082653*xLimUp1R2)*fldSrc[3]*yLimLoR3+(180.0*fldSrc[1]*xLimLo1R2-180.0*fldSrc[1]*xLimUp1R2)*yLimLoR2)*ycLimR2+((311.7691453623978*xLimLo1R2-311.7691453623978*xLimUp1R2)*fldSrc[3]*yLimUpR4+(((120.0*xLimLo[0]*xLimLo[1]-120.0*xLimUp[0]*xLimUp[1])*fldSrc[3]+(97.97958971132716*xLimLo[1]-97.97958971132716*xLimUp[1])*fldSrc[2])*dyLim-240.0*fldSrc[1]*xLimUp1R2+240.0*fldSrc[1]*xLimLo1R2)*yLimUpR3+(((-103.9230484541326*xLimUp[0]*fldSrc[1])-84.85281374238573*fldSrc[0])*xLimUp[1]+(103.9230484541326*xLimLo[0]*fldSrc[1]+84.85281374238573*fldSrc[0])*xLimLo[1])*dyLim*yLimUpR2+(311.7691453623978*xLimUp1R2-311.7691453623978*xLimLo1R2)*fldSrc[3]*yLimLoR4+(((120.0*xLimUp[0]*xLimUp[1]-120.0*xLimLo[0]*xLimLo[1])*fldSrc[3]+(97.97958971132716*xLimUp[1]-97.97958971132716*xLimLo[1])*fldSrc[2])*dyLim+240.0*fldSrc[1]*xLimUp1R2-240.0*fldSrc[1]*xLimLo1R2)*yLimLoR3+((103.9230484541326*xLimUp[0]*fldSrc[1]+84.85281374238573*fldSrc[0])*xLimUp[1]+((-103.9230484541326*xLimLo[0]*fldSrc[1])-84.85281374238573*fldSrc[0])*xLimLo[1])*dyLim*yLimLoR2)*ycLim+(124.7076581449591*xLimUp1R2-124.7076581449591*xLimLo1R2)*fldSrc[3]*yLimUpR5+(((90.0*xLimUp[0]*xLimUp[1]-90.0*xLimLo[0]*xLimLo[1])*fldSrc[3]+(73.48469228349535*xLimUp[1]-73.48469228349535*xLimLo[1])*fldSrc[2])*dyLim+90.0*fldSrc[1]*xLimUp1R2-90.0*fldSrc[1]*xLimLo1R2)*yLimUpR4+(((17.32050807568877*xLimUp0R2-17.32050807568877*xLimLo0R2)*fldSrc[3]+(28.28427124746191*xLimUp[0]-28.28427124746191*xLimLo[0])*fldSrc[2])*dyLimR2+((69.28203230275508*xLimUp[0]*fldSrc[1]+56.56854249492383*fldSrc[0])*xLimUp[1]+((-69.28203230275508*xLimLo[0]*fldSrc[1])-56.56854249492383*fldSrc[0])*xLimLo[1])*dyLim)*yLimUpR3+((15.0*xLimUp0R2-15.0*xLimLo0R2)*fldSrc[1]+24.49489742783179*fldSrc[0]*xLimUp[0]-24.49489742783179*fldSrc[0]*xLimLo[0])*dyLimR2*yLimUpR2+(124.7076581449591*xLimLo1R2-124.7076581449591*xLimUp1R2)*fldSrc[3]*yLimLoR5+(((90.0*xLimLo[0]*xLimLo[1]-90.0*xLimUp[0]*xLimUp[1])*fldSrc[3]+(73.48469228349535*xLimLo[1]-73.48469228349535*xLimUp[1])*fldSrc[2])*dyLim-90.0*fldSrc[1]*xLimUp1R2+90.0*fldSrc[1]*xLimLo1R2)*yLimLoR4+(((17.32050807568877*xLimLo0R2-17.32050807568877*xLimUp0R2)*fldSrc[3]+(28.28427124746191*xLimLo[0]-28.28427124746191*xLimUp[0])*fldSrc[2])*dyLimR2+(((-69.28203230275508*xLimUp[0]*fldSrc[1])-56.56854249492383*fldSrc[0])*xLimUp[1]+(69.28203230275508*xLimLo[0]*fldSrc[1]+56.56854249492383*fldSrc[0])*xLimLo[1])*dyLim)*yLimLoR3+((15.0*xLimLo0R2-15.0*xLimUp0R2)*fldSrc[1]-24.49489742783179*fldSrc[0]*xLimUp[0]+24.49489742783179*fldSrc[0]*xLimLo[0])*dyLimR2*yLimLoR2))/dyLimR2; 
  fldDest[3] += -(0.00625*(((587.877538267963*xLimUp1R3-587.877538267963*xLimLo1R3)*fldSrc[3]*yLimUpR3+(509.1168824543145*fldSrc[1]*xLimUp1R3-509.1168824543145*fldSrc[1]*xLimLo1R3)*yLimUpR2+(587.877538267963*xLimLo1R3-587.877538267963*xLimUp1R3)*fldSrc[3]*yLimLoR3+(509.1168824543145*fldSrc[1]*xLimLo1R3-509.1168824543145*fldSrc[1]*xLimUp1R3)*yLimLoR2)*ycLimR3+((1322.724461102916*xLimLo1R3-1322.724461102916*xLimUp1R3)*fldSrc[3]*yLimUpR4+(((509.1168824543145*xLimLo[0]*xLimLo1R2-509.1168824543145*xLimUp[0]*xLimUp1R2)*fldSrc[3]+(207.8460969082653*xLimLo1R2-207.8460969082653*xLimUp1R2)*fldSrc[2])*dyLim-1018.233764908629*fldSrc[1]*xLimUp1R3+1018.233764908629*fldSrc[1]*xLimLo1R3)*yLimUpR3+(((-440.9081537009721*xLimUp[0]*fldSrc[1])-180.0*fldSrc[0])*xLimUp1R2+(440.9081537009721*xLimLo[0]*fldSrc[1]+180.0*fldSrc[0])*xLimLo1R2)*dyLim*yLimUpR2+(1322.724461102916*xLimUp1R3-1322.724461102916*xLimLo1R3)*fldSrc[3]*yLimLoR4+(((509.1168824543145*xLimUp[0]*xLimUp1R2-509.1168824543145*xLimLo[0]*xLimLo1R2)*fldSrc[3]+(207.8460969082653*xLimUp1R2-207.8460969082653*xLimLo1R2)*fldSrc[2])*dyLim+1018.233764908629*fldSrc[1]*xLimUp1R3-1018.233764908629*fldSrc[1]*xLimLo1R3)*yLimLoR3+((440.9081537009721*xLimUp[0]*fldSrc[1]+180.0*fldSrc[0])*xLimUp1R2+((-440.9081537009721*xLimLo[0]*fldSrc[1])-180.0*fldSrc[0])*xLimLo1R2)*dyLim*yLimLoR2)*ycLimR2+((1058.179568882333*xLimUp1R3-1058.179568882333*xLimLo1R3)*fldSrc[3]*yLimUpR5+(((763.6753236814716*xLimUp[0]*xLimUp1R2-763.6753236814716*xLimLo[0]*xLimLo1R2)*fldSrc[3]+(311.7691453623978*xLimUp1R2-311.7691453623978*xLimLo1R2)*fldSrc[2])*dyLim+763.6753236814716*fldSrc[1]*xLimUp1R3-763.6753236814716*fldSrc[1]*xLimLo1R3)*yLimUpR4+(((146.9693845669907*xLimUp0R2*xLimUp[1]-146.9693845669907*xLimLo0R2*xLimLo[1])*fldSrc[3]+(120.0*xLimUp[0]*xLimUp[1]-120.0*xLimLo[0]*xLimLo[1])*fldSrc[2])*dyLimR2+((587.877538267963*xLimUp[0]*fldSrc[1]+240.0*fldSrc[0])*xLimUp1R2+((-587.877538267963*xLimLo[0]*fldSrc[1])-240.0*fldSrc[0])*xLimLo1R2)*dyLim)*yLimUpR3+((127.2792206135786*xLimUp0R2*fldSrc[1]+103.9230484541326*fldSrc[0]*xLimUp[0])*xLimUp[1]+((-127.2792206135786*xLimLo0R2*fldSrc[1])-103.9230484541326*fldSrc[0]*xLimLo[0])*xLimLo[1])*dyLimR2*yLimUpR2+(1058.179568882333*xLimLo1R3-1058.179568882333*xLimUp1R3)*fldSrc[3]*yLimLoR5+(((763.6753236814716*xLimLo[0]*xLimLo1R2-763.6753236814716*xLimUp[0]*xLimUp1R2)*fldSrc[3]+(311.7691453623978*xLimLo1R2-311.7691453623978*xLimUp1R2)*fldSrc[2])*dyLim-763.6753236814716*fldSrc[1]*xLimUp1R3+763.6753236814716*fldSrc[1]*xLimLo1R3)*yLimLoR4+(((146.9693845669907*xLimLo0R2*xLimLo[1]-146.9693845669907*xLimUp0R2*xLimUp[1])*fldSrc[3]+(120.0*xLimLo[0]*xLimLo[1]-120.0*xLimUp[0]*xLimUp[1])*fldSrc[2])*dyLimR2+(((-587.877538267963*xLimUp[0]*fldSrc[1])-240.0*fldSrc[0])*xLimUp1R2+(587.877538267963*xLimLo[0]*fldSrc[1]+240.0*fldSrc[0])*xLimLo1R2)*dyLim)*yLimLoR3+(((-127.2792206135786*xLimUp0R2*fldSrc[1])-103.9230484541326*fldSrc[0]*xLimUp[0])*xLimUp[1]+(127.2792206135786*xLimLo0R2*fldSrc[1]+103.9230484541326*fldSrc[0]*xLimLo[0])*xLimLo[1])*dyLimR2*yLimLoR2)*ycLim+(293.9387691339815*xLimLo1R3-293.9387691339815*xLimUp1R3)*fldSrc[3]*yLimUpR6+(((305.4701294725887*xLimLo[0]*xLimLo1R2-305.4701294725887*xLimUp[0]*xLimUp1R2)*fldSrc[3]+(124.7076581449591*xLimLo1R2-124.7076581449591*xLimUp1R2)*fldSrc[2])*dyLim-203.6467529817258*fldSrc[1]*xLimUp1R3+203.6467529817258*fldSrc[1]*xLimLo1R3)*yLimUpR5+(((110.227038425243*xLimLo0R2*xLimLo[1]-110.227038425243*xLimUp0R2*xLimUp[1])*fldSrc[3]+(90.0*xLimLo[0]*xLimLo[1]-90.0*xLimUp[0]*xLimUp[1])*fldSrc[2])*dyLimR2+(((-220.454076850486*xLimUp[0]*fldSrc[1])-90.0*fldSrc[0])*xLimUp1R2+(220.454076850486*xLimLo[0]*fldSrc[1]+90.0*fldSrc[0])*xLimLo1R2)*dyLim)*yLimUpR4+(((14.14213562373095*xLimLo0R3-14.14213562373095*xLimUp0R3)*fldSrc[3]+(17.32050807568877*xLimLo0R2-17.32050807568877*xLimUp0R2)*fldSrc[2])*dyLimR3+(((-84.85281374238573*xLimUp0R2*fldSrc[1])-69.28203230275508*fldSrc[0]*xLimUp[0])*xLimUp[1]+(84.85281374238573*xLimLo0R2*fldSrc[1]+69.28203230275508*fldSrc[0]*xLimLo[0])*xLimLo[1])*dyLimR2)*yLimUpR3+((12.24744871391589*xLimLo0R3-12.24744871391589*xLimUp0R3)*fldSrc[1]-15.0*fldSrc[0]*xLimUp0R2+15.0*fldSrc[0]*xLimLo0R2)*dyLimR3*yLimUpR2+(293.9387691339815*xLimUp1R3-293.9387691339815*xLimLo1R3)*fldSrc[3]*yLimLoR6+(((305.4701294725887*xLimUp[0]*xLimUp1R2-305.4701294725887*xLimLo[0]*xLimLo1R2)*fldSrc[3]+(124.7076581449591*xLimUp1R2-124.7076581449591*xLimLo1R2)*fldSrc[2])*dyLim+203.6467529817258*fldSrc[1]*xLimUp1R3-203.6467529817258*fldSrc[1]*xLimLo1R3)*yLimLoR5+(((110.227038425243*xLimUp0R2*xLimUp[1]-110.227038425243*xLimLo0R2*xLimLo[1])*fldSrc[3]+(90.0*xLimUp[0]*xLimUp[1]-90.0*xLimLo[0]*xLimLo[1])*fldSrc[2])*dyLimR2+((220.454076850486*xLimUp[0]*fldSrc[1]+90.0*fldSrc[0])*xLimUp1R2+((-220.454076850486*xLimLo[0]*fldSrc[1])-90.0*fldSrc[0])*xLimLo1R2)*dyLim)*yLimLoR4+(((14.14213562373095*xLimUp0R3-14.14213562373095*xLimLo0R3)*fldSrc[3]+(17.32050807568877*xLimUp0R2-17.32050807568877*xLimLo0R2)*fldSrc[2])*dyLimR3+((84.85281374238573*xLimUp0R2*fldSrc[1]+69.28203230275508*fldSrc[0]*xLimUp[0])*xLimUp[1]+((-84.85281374238573*xLimLo0R2*fldSrc[1])-69.28203230275508*fldSrc[0]*xLimLo[0])*xLimLo[1])*dyLimR2)*yLimLoR3+((12.24744871391589*xLimUp0R3-12.24744871391589*xLimLo0R3)*fldSrc[1]+15.0*fldSrc[0]*xLimUp0R2-15.0*fldSrc[0]*xLimLo0R2)*dyLimR3*yLimLoR2))/dyLimR3; 

}

