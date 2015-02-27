#python plot.py -d Ifront1_1.000.dat

import numpy as np
import pylab as pl
import sys
from optparse import OptionParser

kpc 	= 3.08567758131e21	#kpc to cm
megapc     = 3.08567758131e24	#Mpc to cm

#Parse command line arguments
parser = OptionParser()
parser.add_option('-d', '--data_file', dest='data', help='C2Ray data', default='results/Ifront1_1.000.dat')
(options,args) = parser.parse_args()

megapc = 3.08567758131e24	#Mpc to cm

pl.figure()

#pl.legend(loc='upper left')

spectrum = np.loadtxt('results/spectrum.dat')
p = np.loadtxt(options.data)

r = p[:,0]/megapc #distance in kpc
xHI = p[:,1]
xHII = p[:,2]
T = p[:,3] 
n = p[:,4] 
photo = p[:,5] 
avg = p[:,6] 
xHeIII = p[:,7] 

pl.subplot(2,1,1)
pl.plot(r, T, color='g',label='normal')
#pl.semilogy(r, xHeIII, color='g',label='normal')
#pl.semilogy(r, xHI, color='r',label='normal')
#pl.semilogy(r, photo, color='r',label='normal')
#pl.legend(loc='best', prop={'size':10})
pl.xlabel('$d \; \mathrm{Mpc}$')
pl.ylabel('Temperature $\mathrm{K}$',fontsize=20)
pl.title('$\eta=0.6$')
pl.subplot(2,1,2)
pl.semilogx(spectrum[:,0],spectrum[:,1],label='intrinsic')
pl.semilogx(spectrum[:,0],spectrum[:,2],label='1.5Mpc')
pl.semilogx(spectrum[:,0],spectrum[:,3],label='3.0Mpc')
pl.semilogx(spectrum[:,0],spectrum[:,4],label='4.5Mpc')
pl.semilogx(spectrum[:,0],spectrum[:,5],label='6.0Mpc')
pl.semilogx(spectrum[:,0],spectrum[:,6],label='7.5Mpc')
pl.semilogx(spectrum[:,0],spectrum[:,7],label='10.5Mpc')
pl.semilogx(spectrum[:,0],spectrum[:,8],label='13.5Mpc')
pl.xlabel('frequency')
pl.ylabel('luminosity')
pl.legend()
pl.show()