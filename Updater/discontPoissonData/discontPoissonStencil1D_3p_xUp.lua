local Lin = require "Lib.Linalg"
local function stencilFn(dx, a, b, val)
  local _M = {}

  _M[1] = Lin.Mat(4,4)
  _M[1][1][1] = (2625.0*b-130.0*a)/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[1][1][2] = (3715.248982235241*b-232.0948082142295*a)/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[1][1][3] = (2985.15074996222*b-272.8002932549743*a)/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[1][1][4] = (9167.528292838808*b-1317.584152910166*a)/(3360.0*dx[1]^2*b+224.0*dx[1]^2*a)
  _M[1][2][1] = -(1.0*(2121.762239271874*b+387.9793808954285*a))/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[1][2][2] = -(1.0*(3003.0*b+588.0*a))/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[1][2][3] = -(1.0*(2412.86862468722*b+542.2176684690385*a))/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[1][2][4] = -(1.0*(1058.574985534799*b+293.2848444771737*a))/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[1][3][1] = (301.8691769624717*b-252.6756814574763*a)/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[1][3][2] = (104.5705503476002*b-422.1551847366085*a)/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[1][3][3] = -(1.0*(495.0*b+455.0*a))/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[1][3][4] = -(1.0*(4738.779906262793*b+2041.047525169368*a))/(1680.0*dx[1]^2*b+112.0*dx[1]^2*a)
  _M[1][4][1] = -(1.0*(4167.05831492673*b+761.9763775866021*a))/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[1][4][2] = -(1.0*(5897.774919408165*b+1154.809075128872*a))/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[1][4][3] = -(1.0*(4738.779906262793*b+1064.894360957931*a))/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[1][4][4] = -(1.0*(2079.0*b+576.0*a))/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[2] = Lin.Mat(4,4)
  _M[2][1][1] = -(1.0*(2625.0*b+830.0*a))/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[2][1][2] = (3715.248982235241*b-1881.0071770198*a)/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[2][1][3] = -(1.0*(2985.15074996222*b+1068.8404932449*a))/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[2][1][4] = (9167.528292838808*b-17573.08020809101*a)/(3360.0*dx[1]^2*b+224.0*dx[1]^2*a)
  _M[2][2][1] = (2121.762239271874*b-387.9793808954285*a)/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[2][2][2] = -(1.0*(5691.0*b+2100.0*a))/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[2][2][3] = (2412.86862468722*b-542.2176684690385*a)/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[2][2][4] = -(1.0*(3698.138585829362*b+2346.27875581739*a))/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[2][3][1] = -(1.0*(301.8691769624717*b+605.974421902443*a))/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[2][3][2] = -(1.0*(546.0906518152458*b+1940.364656449916*a))/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[2][3][3] = -(1.0*(8505.0*b+1345.0*a))/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[2][3][4] = -(1.0*(42223.06141198197*b+19079.35730049627*a))/(1680.0*dx[1]^2*b+112.0*dx[1]^2*a)
  _M[2][4][1] = (4167.05831492673*b-761.9763775866021*a)/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[2][4][2] = -(1.0*(11176.90211999729*b+4124.318125460255*a))/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[2][4][3] = (4738.779906262793*b-1064.894360957931*a)/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[2][4][4] = -(1.0*(50463.0*b+7488.0*a))/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[3] = Lin.Vec(4)
  _M[3][1] = 1357.645019878173/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)*val
  _M[3][2] = 1097.371404766865/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)*val
  _M[3][3] = 1214.31462150466/(240.0*dx[1]^2*b+16.0*dx[1]^2*a)*val
  _M[3][4] = 2155.194654781794/(480.0*dx[1]^2*b+32.0*dx[1]^2*a)*val
  return(_M)
end

return(stencilFn)