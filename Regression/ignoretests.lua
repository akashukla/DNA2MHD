-- File with list of regression tests to ignore: this is mainly used
-- as a stop-gap while some things are being debugged.
return {
   "./dg-maxwell/rt-maxwell-1d-plane-p3.lua",
   "./dg-maxwell/rt-maxwell-2d-plane-p3.lua",
   "./dg-maxwell/rt-maxwell-2d-reflect-p3.lua",
   "./dg-maxwell/rt-maxwell-3d-plane-p3.lua",
   "./gk-sheath/rt-gk-sheath-solovev-R.lua",
   "./mgPoisson/rt-mgPoisson-FEM-wholeSolver-3D.lua",
   "./mom-ecdi1d/rt-5m-ecdi1d.lua",
   "./vm-const-mag/rt-const-mag-1x2v-ms-p3.lua",
   "./vm-const-mag/rt-par-const-mag-1x2v-ms-p1.lua",
   "./vm-lbo/rt-conservation-collisions-1x1v-p3.lua",
   "./vm-lbo/rt-lboRelax-1x1v-p3.lua",
   "./vm-lbo/rt-lboRelax-1x2v-p3.lua",
   "./vm-lbo/rt-lboRelax-2x2v-p3.lua",
   "./vm-mom/rt-1x1v-p3-ser-mom.lua",
   "./vm-mom/rt-1x2v-p3-ser-mom.lua",
   "./vm-mom/rt-1x3v-p3-ser-mom.lua",
   "./vm-mom/rt-2x2v-p3-ser-mom.lua",
   "./vm-mom/rt-2x3v-p3-ser-mom.lua",
   "./vm-two-stream/rt-two-stream-p3.lua",
   "./vm-weibel/rt-weibel-1x2v-p3.lua",
   "./vm-weibel/rt-weibel-2x2v-p3.lua",
   "./mom-axisymmetric/rt-phmaxwell-circular-waveguide-TM.lua",
   "./mom-axisymmetric/rt-phmaxwell-circular-waveguide-TE.lua",
   "./mom-axisymmetric/rt-phmaxwell-axisymmetric.lua",
   "./mom-axisymmetric/rt-5m-dean-axisymmetric.lua",
   "./mom-axisymmetric/rt-5m-mirror.lua",
   "./gk-neutrals/rt-recycleBCs-1x2v-p2.lua",
   "./gk-neutrals/rt-nstx-neut-lowRes-p2.lua",
}
