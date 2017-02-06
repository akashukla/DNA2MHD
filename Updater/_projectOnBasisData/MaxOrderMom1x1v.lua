local _M = { numDensity = {}, momentum = {} } 
 
_M.numDensity[1] = function (f, out, dv, w) 
   out[1] = out[1] + 1.414213562373095*f[1] 
   out[2] = out[2] + 1.414213562373095*f[2] 
end 
_M.numDensity[2] = function (f, out, dv, w) 
   out[1] = out[1] + 1.414213562373095*f[1] 
   out[2] = out[2] + 1.414213562373095*f[2] 
   out[3] = out[3] + 1.414213562373095*f[5] 
end 
_M.numDensity[3] = function (f, out, dv, w) 
   out[1] = out[1] + 1.414213562373095*f[1] 
   out[2] = out[2] + 1.414213562373095*f[2] 
   out[3] = out[3] + 1.414213562373095*f[5] 
   out[4] = out[4] + 1.414213562373095*f[9] 
end 
_M.numDensity[4] = function (f, out, dv, w) 
   out[1] = out[1] + 1.414213562373095*f[1] 
   out[2] = out[2] + 1.414213562373095*f[2] 
   out[3] = out[3] + 1.414213562373095*f[5] 
   out[4] = out[4] + 1.414213562373095*f[9] 
   out[5] = out[5] + 1.414213562373095*f[14] 
end 

 
_M.momentum[1] = function (f, out, dv, w) 
   local dv1 = dv[1] 
   local w1 = w[1] 
   out[1] = out[1] + 1.414213562373095*f[1]*w1+0.408248290463863*f[3]*dv1 
   out[2] = out[2] + 1.414213562373095*f[2]*w1 
end 
_M.momentum[2] = function (f, out, dv, w) 
   local dv1 = dv[1] 
   local w1 = w[1] 
   out[1] = out[1] + 1.414213562373095*f[1]*w1+0.408248290463863*f[3]*dv1 
   out[2] = out[2] + 1.414213562373095*f[2]*w1+0.408248290463863*f[4]*dv1 
   out[3] = out[3] + 1.414213562373095*f[5]*w1 
end 
_M.momentum[3] = function (f, out, dv, w) 
   local dv1 = dv[1] 
   local w1 = w[1] 
   out[1] = out[1] + 1.414213562373095*f[1]*w1+0.408248290463863*f[3]*dv1 
   out[2] = out[2] + 1.414213562373095*f[2]*w1+0.408248290463863*f[4]*dv1 
   out[3] = out[3] + 1.414213562373095*f[5]*w1+0.408248290463863*f[7]*dv1 
   out[4] = out[4] + 1.414213562373095*f[9]*w1 
end 
_M.momentum[4] = function (f, out, dv, w) 
   local dv1 = dv[1] 
   local w1 = w[1] 
   out[1] = out[1] + 1.414213562373095*f[1]*w1+0.408248290463863*f[3]*dv1 
   out[2] = out[2] + 1.414213562373095*f[2]*w1+0.408248290463863*f[4]*dv1 
   out[3] = out[3] + 1.414213562373095*f[5]*w1+0.408248290463863*f[7]*dv1 
   out[4] = out[4] + 1.414213562373095*f[9]*w1+0.408248290463863*f[12]*dv1 
   out[5] = out[5] + 1.414213562373095*f[14]*w1 
end 
return _M 
