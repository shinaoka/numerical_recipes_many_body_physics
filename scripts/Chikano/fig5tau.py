import numpy
import numpy as np
import irbasis
import scipy
import scipy.integrate as integrate
from matplotlib import pyplot as plt


#import cmath
#import math
#\matplotlib inline
plt.rc('text', usetex=True)

mk = ["o","^","v"]
i = 0


def safe_exp(x,beta,stat):
    
    
    if stat == "f":
        return -(1/2)*(numpy.exp(-x)+numpy.exp(x-beta))/(1+numpy.exp(-beta))
    elif stat == "b":
        return -(1/2)*(numpy.exp(-x)+numpy.exp(x-beta))/(1-numpy.exp(-beta))           


def metal_sp(x,omega,stat,beta):
    if stat == "f":
            if omega < 0 :
                return -numpy.sqrt(1-omega**2)*numpy.exp((beta-x)*omega)/(1+numpy.exp(beta*omega))
                    
                    
            else:
                return -numpy.sqrt(1-omega**2)*numpy.exp(-x*omega)/(1+numpy.exp(-beta*omega))
            
    elif stat == "b":
        if 10**-7 < numpy.abs(beta*omega):
            if omega <0:
                return -omega*numpy.sqrt(1-omega**2)*numpy.exp((beta-x)*omega)/(-1+numpy.exp(beta*omega))
            else:
                return -omega*numpy.sqrt(1-omega**2)*numpy.exp(-x*omega)/(1-numpy.exp(-beta*omega))
            
        else:
            if omega <0:
                return -numpy.sqrt(1-omega**2)*numpy.exp((beta-x)*omega)/beta
            else:
                return -numpy.sqrt(1-omega**2)*numpy.exp(-x*omega)/beta
        

def metal(x,beta,stat):
    gt = scipy.integrate.quad(lambda omega:metal_sp(x,omega,stat,beta),-1,1,limit = 500,points = [-1/2,-0.01,0.01,1/2]) 
    
    if stat =="f":
        return 2*gt[0]/numpy.pi #0.392699 
    elif stat =="b":
        return 3*gt[0]/2


for wmax in [1]:
    beta = 100.0
    pole = 1
    Lambda = beta *wmax
    
    stat = "F"
    #basis = irbasis.load(stat, Lambda)#,"irbasis.h5")
    b = irbasis.load('F', Lambda)
    dim = b.dim()
    
    
    #model = "Semicircular"
    model = "Insulator"
    
    whichl =int(dim -1)
    l_list = numpy.arange(whichl+1)
    
    gl = numpy.empty(whichl+1)
    rhol = []
    Sl = []
    mat_u = numpy.zeros([whichl+1,whichl+1])
    
    samp = b.sampling_points_x(dim-1)
    
    sampling_point = beta *(samp +1)/2
    
    gtau = []
    gtausub = []
    
    for l in range(whichl+1):

        if stat == "F":      
            if model == "Semicircular":
                rhol.append(scipy.integrate.quad(lambda omega: 2/numpy.pi *numpy.sqrt(1-omega**2)*numpy.sqrt(1/wmax)*b.vly(l,omega/wmax) ,-1,1,limit = 1500)[0] )
                
                
            #np.sqrt(wmax)*b.Vlomega(l,-pole/wmax)
                
            elif model == "Insulator":    
                rhol.append((b.vly(l,pole))/2+b.vly(l,-pole)/2)     
                

        gl[l] =-numpy.sqrt(beta*wmax/2)*b.sl(l) * rhol[l]
        
        
        len(sampling_point)
        
        for n,tau in enumerate(sampling_point):

            mat_u[n,l] = numpy.sqrt(2/beta)*b.ulx(l,2*tau/beta -1)
          
            
    for tau in sampling_point:
        print(2*tau/beta -1)
        #gtau.append(metal(tau,beta,"f"))
        gtau.append(safe_exp(tau,beta,"f"))
    

    ta = numpy.linspace(0,beta,500)    
    for tau in ta:
        #print(2*tau/beta -1)
        #gtau.append(metal(tau,beta,"f"))
        gtausub.append(safe_exp(tau,beta,"f"))
        
    gtau = numpy.array(gtau)
    Gl_sampling = numpy.linalg.inv(mat_u) @ gtau
    

    print(numpy.abs(gl)-numpy.abs(Gl_sampling))
    
  
    
    plt.figure(2)
    
    lines = ['-+g','-+b','-xr','-*y','om']
    
    #plt.plot(ta,numpy.array([-safe_exp(tau,beta,"f") for tau in ta]),lines[1],mew = 1,label = r"$\rm{Fermion}$")
    
    goftauir = numpy.zeros(len(ta))
    for l in range(whichl+1):
        cttau = 0
        for tau in ta:
            #print(cttau)
            goftauir[cttau] += Gl_sampling[l]*numpy.sqrt(2/beta)*b.ulx(l,2*tau/beta -1)#leg(l,x)#
            
            #goftauir[cttau] += Gl_sampling[l]*b.Ultau(l,tau)
            
            
            cttau += 1
    plt.plot(ta,numpy.abs(-goftauir),mew = 1,label = r"$-G(\tau)$")
    plt.plot(ta,numpy.abs(-goftauir +gtausub),"--",mew = 1,label = r"$|\rm{\Delta G(\tau)}|$")
    plt.scatter(sampling_point,-gtau,marker="x",color = "r",s=74)
        


plt.figure(1)

#plt.yscale("log")
plt.xlim(-5,beta+5)
#plt.ylim(-6,6)
plt.plot(np.linspace(-5,2*beta,100),np.zeros(100),"--k",mew = 1,alpha=0.5,linewidth = 1)


plt.plot(np.linspace(0,Lambda,1000), [b.ulx(dim-1,2*tau/beta -1) for tau in np.linspace(0,Lambda,1000)],color="b") 
plt.scatter(sampling_point, [b.ulx(dim-1,2*tau/beta -1) for tau in sampling_point],marker="x",color="r") 

#[(2*tau/beta -1) for tau in sampling_point]


#plt.title(r"${}$".format(model),fontsize=21)
plt.ylabel(r"$|G_l|$",fontsize=21)
plt.ylabel(r"$U_{39}^{\rm{F}}(\tau)$",fontsize=21)
plt.tick_params(labelsize=21)
plt.xlabel(r'$\tau$',fontsize = 21)
#plt.legend(frameon=False,fontsize = 21)
plt.legend(loc='upper right',
           bbox_to_anchor=(-0.2, -0.05, 0.6, .600), 
           borderaxespad=0.,frameon=False,fontsize = 21)

plt.tight_layout()


plt.savefig('Gltau100'+'.pdf')

#plt.savefig('optLambda10000'+'.pdf')

plt.figure(2)
plt.xlabel(r'$\tau$',fontsize = 21)

plt.yscale("log")
#plt.xticks([0,beta/5,2*beta/5,3*beta/5,4*beta/5,beta])
#plt.xticks([0,beta],[0,r"$\beta$"])
#plt.yticks([10**-1,10**-2,10**0],[r'$10^{-1}$',r'$10^{-2}$',r'$10^{0}$'])
#plt.xlim(-5,beta+5)
plt.xlim(-5,105)
plt.ylim(10**-16,10)
#plt.title(mat,fontsize = 20)
#plt.ylabel(r'$-G^{\rm{F}}(\tau)$',fontsize = 33)
plt.ylabel(r'$-G(\tau)$',fontsize = 23)
#plt.ylabel(r'$-G^{\alpha}(\tau)$',fontsize = 21)
#plt.yscale("log")
plt.tick_params(labelsize=21)
plt.legend(frameon=False,fontsize = 21)
"""
plt.legend(loc='upper right',
           bbox_to_anchor=(0.4, 0.4, 0.6, .600), 
           borderaxespad=0.,frameon=False,fontsize = 15)
"""
plt.tight_layout()
plt.savefig('Goftau100'+'.pdf')

plt.show()




    
