import numpy as np
import matplotlib.pyplot as plt

lpath = '/scratch/08929/echansen/dna2mhdrun741'
import dna2mhd_utils as dn

te,e = dn.plot_energy(lpath,3,show=False)
e24 =e[:,0]

lpath = '/scratch/08929/echansen/dna2mhdrun742'
te,e = dn.plot_energy(lpath,3,show=False)
e32 = e[:,0]

lpath =	'/scratch/08929/echansen/dna2mhdrun743'
te,e = dn.plot_energy(lpath,3,show=False)
e64 = e[:,0]

epred = np.pi**3 / 2 * np.exp(-0.01 * te)

plt.figure()
plt.plot(te,epred,"k-",label="Predicted")
plt.plot(te,e24,"r+",label="$24^3$")
plt.plot(te,e32,"bo",label="$32^3$")
plt.plot(te,e64,"m*",label="$64^3$")
plt.title("Evolution of Taylor Green Vortex Energy")
plt.xlabel("Time ($\omega_c^{-1}$)")
plt.ylabel("Energy")
plt.legend()
plt.savefig("tgen")
plt.close()

plt.figure()
plt.plot(te,np.abs(e24-epred),"r+",label="$24^3$")
plt.plot(te,np.abs(e32-epred),"bo",label="$32^3$")
plt.plot(te,np.abs(e64-epred),"m*",label="$64^3$")
plt.title("Evolution of Deviation from Taylor Green Vortex Energy")
plt.xlabel("Time ($\omega_c^{-1}$)")
plt.ylabel("Energy Deviation")
plt.legend()
plt.savefig("tgendev")
plt.close()
