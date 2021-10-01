import numpy
import scipy
import scipy.integrate as integrate
from matplotlib import pyplot as plt
#import cmath
#import math
#\matplotlib inline
plt.rc('text', usetex=True)

l_list = [10,100,1000,10000,100000,1000000]

#Insu cut =1e-8

plt.figure(1,figsize = (9,4))
plt.xscale("log")
plt.yscale("log")
plt.xlim(8,1.1*1E+6)
plt.ylim(10,10000)


plt.scatter(l_list,[14,26,40,52,60,68],marker="^",label= r"${\rm IR}$")

plt.scatter(l_list,[19, 49, 148, 469, 1358, 4286],marker="o",label=r"${\rm  Legendre}$")
plt.scatter(l_list,[15,43,131,407,1250,3822],marker="x",label=r"${\rm Chebyshev}$")


plt.scatter([10,100,1000],[350,3498,34973],marker="v",label=r"$G(i\omega_n)\leq10^{-8}$")




#plt.title(r"${}$".format(model),fontsize=21)
plt.ylabel(r"$N$",fontsize=21)

plt.tick_params(labelsize=21)
plt.xlabel(r'$\beta$',fontsize = 21)
#plt.legend(fontsize = 15)

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=16)
"""
plt.legend(loc='upper right',
           bbox_to_anchor=(0.2, 0.4, 0.6, .600), 
           borderaxespad=0.,frameon=False,fontsize = 18)

plt.legend(loc='upper right',
           bbox_to_anchor=(1.2, 0.4, 0.6, .600), 
           borderaxespad=0.,fontsize = 18)
"""

plt.tight_layout()


plt.savefig('1e-8'+'.pdf')


plt.show()




    
