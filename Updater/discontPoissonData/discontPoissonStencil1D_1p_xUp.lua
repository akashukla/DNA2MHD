local Lin = require "Lib.Linalg"
local function stencilFn(dx, a, b, val)
  local _M = {}

  _M[1] = Lin.Mat(2,2)
  _M[1][1][1] = (324.0*b+24.0*a)/(144.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[1][1][2] = (155.8845726811989*b+6.928203230275509*a)/(72.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[1][2][1] = -(1.0*(394.907584125704*b+152.4204710660612*a))/(144.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[1][2][2] = -(1.0*(174.0*b+76.0*a))/(72.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[2] = Lin.Mat(2,2)
  _M[2][1][1] = -(1.0*(324.0*b+312.0*a))/(144.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[2][1][2] = (155.8845726811989*b-214.7743001385408*a)/(72.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[2][2][1] = (394.907584125704*b-235.5589098293673*a)/(144.0*dx[1]^2*b+32.0*dx[1]^2*a)
  _M[2][2][2] = -(1.0*(846.0*b+524.0*a))/(72.0*dx[1]^2*b+16.0*dx[1]^2*a)
  _M[3] = Lin.Vec(2)
  _M[3][1] = 144.0/(50.91168824543144*dx[1]^2*b+11.31370849898477*dx[1]^2*a)*val
  _M[3][2] = 193.9896904477143/(50.91168824543144*dx[1]^2*b+11.31370849898477*dx[1]^2*a)*val
  return(_M)
end

return(stencilFn)