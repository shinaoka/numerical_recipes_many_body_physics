import numpy
import numpy as np
from numpy.random import *
import irbasis
from irbasis_util.two_point_basis import Basis, sampling_points_matsubara, sampling_points_tau


import scipy
import scipy.integrate as integrate
from matplotlib import pyplot as plt
import cmath
import math
#\matplotlib inline
plt.rc('text', usetex=True)




markers = ['o', 'x',  '+', 'v']
lines = ['-o', '-x',  '-+', '-v']
k = 0


for Lambda in [100.0]:    
    #Lambda = 100.0
    wmax = 1
    beta = Lambda / wmax
    pole = 1
    stat = "F"
    
    #basis = irbasis.load('F', Lambda,"irbasis.h5")
    b = Basis(irbasis.load('F', Lambda), beta)
    dim = b.dim
    whichl = int(b.dim -1)
    
    
    NN = 100000
    y = []
    newl = int(whichl)
    
    
    sl = []
    rhol = []
    gl = []
    llist = []
    for l in range(whichl):
        if l%2 == 1:
            continue
        #print(l)
        
        
        if stat == 'B':
            llist.append(l)
    

            rhol.append(basis.vly(l,pole)/2-basis.vly(l,-pole)/2)          
            
            #rhol.append(np.sqrt(wmax)*b.Vlomega(l,pole/wmax)/2   - np.sqrt(wmax)*b.Vlomega(l,-pole/wmax)/2)
            
            
            
        elif stat == 'F':
            llist.append(l)
            
            
            #if model == "Semicircular":
            rhol.append(scipy.integrate.quad(lambda omega: numpy.sqrt(1-omega**2)*np.sqrt(wmax)*b.Vlomega(l,omega) ,-1,1,limit = 1500)[0] )
            
            
            #rhol.append(np.sqrt(wmax)*b.Vlomega(l,pole/wmax)/2   + np.sqrt(wmax)*b.Vlomega(l,-pole/wmax)/2)
            
            
           #np.sqrt((beta*wmax)/2)*b.Sl  #for F
        #gl[l] =-basis.sl(l) * rhol[l]
        
        #gl[l] =-np.sqrt(2/(beta*wmax))*b.Sl(l)  * rhol[l]
        gl.append(-np.sqrt(2/(beta*wmax))*b.Sl(l)  * rhol[int(l/2)])
        sl.append(b.Sl(l))
    
    #plt.plot(llist,numpy.abs(gl),lines[k+2],label=r"$\beta = {}$".format(int(beta)))
    #plt.plot(llist,numpy.abs(gl),lines[k+2],label=r"$|G_l| {0}\beta = {1}$".format(r"{\hspace{3mm}}",int(beta)))
    plt.plot(llist,numpy.abs(gl),lines[k+2],label=r"$|G_l|$")
    
    #plt.plot(llist,numpy.abs(rhol),lines[k],label=r"$ \beta = {}$".format(int(beta)))
    #plt.plot(llist,numpy.abs(rhol),lines[k],label=r"$|\rho_l |{0} \beta = {1}$".format(r"{\hspace{3mm}}",int(beta)))
    plt.plot(llist,numpy.abs(rhol),lines[k],label=r"$|\rho_l |$")
    
    
    #plt.plot(llist,sl,lines[k],label=r"$ \beta = {0}$".format(int(beta)))
    k+=1
    





plt.yscale("log")
plt.xlim(-0.1,30.1)
plt.ylim(10**-10,10)
#plt.ylim(10**-1,10)

#plt.ylabel(r"$|G_l|$",fontsize=21)
#plt.ylabel(r"$G_l$",fontsize=21)
#plt.yticks([1.0,1e-4,1e-8,1e-12],[r"$1$",r"$10^{-4}$",r"$10^{-8}$",r"$10^{-12}$"])


plt.tick_params(labelsize=21)
#plt.ylabel(r'$S_l$',fontsize = 21)
#plt.ylabel(r'$|\rho_l|$',fontsize = 21)
#plt.ylabel(r'$|G_l - G_l^{lsq}|$',fontsize = 21)
plt.xlabel(r'$l$',fontsize = 21)
plt.legend(frameon=False,fontsize = 21)
plt.tight_layout()

plt.savefig('Glrholmetal100'+'.pdf')
#plt.savefig('Glmetal100'+'.pdf')
#plt.savefig('Sl'+'.pdf')
plt.show()




    
