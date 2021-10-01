import numpy
import numpy as np
import irbasis
import scipy
#from func import *
import scipy.integrate as integrate
from matplotlib import pyplot as plt
import cmath
import math

#from irbasis_util.two_point_basis import Basis, sampling_points_matsubara, sampling_points_tau




#\matplotlib inline
plt.rc('text', usetex=True)


Lambda = 100.0
wmax = 1
beta = Lambda / wmax
pole = 1
stat = "F"

b = irbasis.load('F', Lambda)
dim = b.dim()


sampling_point  = b.sampling_points_matsubara(dim-1)





for whichl in [int(dim -1)]:

    x1 = numpy.arange(100)
    x2 = numpy.array(numpy.exp(numpy.linspace(numpy.log(100), numpy.log(1E+7), 500)), dtype=int)
    
    x = numpy.unique(numpy.hstack((x1,x2)))
    #x = 2*np.pi*numpy.unique(numpy.hstack((x1,x2)))/beta
    
    mx = -x
    x = numpy.unique(numpy.hstack((mx,x)))
    
    #unl = basis.compute_unl(x)
    #unl =  b.compute_Unl(x)/np.sqrt(beta)
    unl =  b.compute_unl(x)
    
    #/////////// compute U(iwn) /////////
    #sp = sampling_points_matsubara(b, whichl=b.dim-1)
    sp  = b.sampling_points_matsubara(dim-1)
    
    #unlsamp = basis.compute_unl(sp)
    unlsamp = b.compute_unl(sp)
    #/////////// compute U(iwn) /////////
    if whichl%2 == 0:
        y = numpy.sqrt(beta)*unl[:,whichl].imag
        ysamp = numpy.sqrt(beta)*unlsamp[:,whichl].imag
    else:
        
        y = numpy.sqrt(beta)*unl[:,whichl].real
        ysamp = numpy.sqrt(beta)*unlsamp[:,whichl].real
        
    print(whichl)
        
    fig = plt.figure(figsize=(6,7))

    
    ax1 = fig.add_subplot(2,1,1)
    #plt.plot(x, y)
    
    x_iwn = []
    x_iwn2 = []
    x_sp = []
    
    for n in x:    
        x_iwn.append(  (2*n +1)* numpy.pi /beta )
        
    for n in np.linspace(-1000,1005,100):    
        x_iwn2.append(  (2*n +1)* numpy.pi /beta )
        
    for n in sp:    
        x_sp.append(  (2*n +1)* numpy.pi /beta )
    
    
    
    #ax1.plot(x, y,"-b")
    ax1.plot(x_iwn, y,"-b")
    ax1.plot(x_iwn2,np.zeros(100),"--k",mew = 1,alpha=0.5,linewidth = 1)
    
    
    ax1.scatter(x_sp,ysamp,marker="x",color="r")
    #ax1.set_xticks([],[])
    ax1.set_xlim(-100,100)
    #plt.ylim(1E-7,100)
    ax1.set_ylabel(r"$\hat{U}_{39}^{\rm{F}}(i\omega^{{\rm F}}_n)$",fontsize=21)
    ax1.tick_params(labelsize=21,labelbottom='off')
    
    shx = 1e-1
    ax1.set_xscale("symlog",linthreshx=shx)
    
    #plt.yscale("log")
    #plt.savefig('unl'+'.pdf')  
    #plt.show()
    
    
    zero_med = b.sampling_points_matsubara(whichl=dim-1)

 

    rhol = []
    gl = numpy.empty(whichl +1)
    llist = []
    for l in range(whichl +1):
                
        if stat == 'B':
            llist.append(l)
    

            rhol.append(b.vly(l,pole)/2-b.vly(l,-pole)/2)          
        elif stat == 'F':
            llist.append(l)

            rhol.append(b.vly(l,pole)/2+b.vly(l,-pole)/2) 
            
        gl[l] =-numpy.sqrt(beta*wmax/2)*b.sl(l) * rhol[l]    
    

    
    

    unl2 =b.compute_unl(zero_med)
    
    
    Umat = numpy.zeros([len(zero_med)*2,whichl+1])
 
    Giwn = numpy.zeros(len(zero_med)*2)    
    iwnset=[]
    nset=[]
    for n in range(len(zero_med)):    
        wn = (2*zero_med[n]+1)* numpy.pi /beta 
        iwnset.append(wn)
        #mono pole
        #Giwn[n]=-pole/(wn**2 + pole**2)#real part
        #Giwn[n + len(zero_med)]=-wn/(wn**2 + pole**2)#imag part   
        
        #two poles
        Giwn[n]= 0 #real part
        Giwn[n + len(zero_med)]=-2*wn/(wn**2 + pole**2)#imag part   
        
        
        
    
    
    l2 = whichl +2
    Umat2 =beta/numpy.sqrt(2) * numpy.vstack((unl2[:][0:l2].real , unl2[:][0:l2].imag))
    
    gre = []
    gim = []
    iwnset2 = []
    
    for n in x:    
        wn = (2*n +1)* numpy.pi /beta 
        iwnset2.append(wn)
        #gre.append(-pole/(wn**2 + pole**2))
        #gim.append(-wn/(wn**2 + pole**2))#imag part          
        #two poles
        gre.append(0)
        gim.append(-wn/(wn**2 + pole**2))#imag part          
        
        
    #plt.scatter(x,gim,marker="o",color = "b")
    
    
    
    
    
    rhol = []
    gl2 = numpy.empty(whichl +1)
    llist = []
    for l in range(whichl +1):
                
        if stat == 'B':
            llist.append(l)
    

            #rhol.append(basis.vly(l,pole)/2-basis.vly(l,-pole)/2)          
        elif stat == 'F':
            llist.append(l)

            rhol.append(b.vly(l,pole)/2+b.vly(l,-pole)/2)           
            
            
        
        gl2[l] =-numpy.sqrt(beta*wmax/2)*b.sl(l) * rhol[l]
        
    

    zero_med=sp*1
    

    unl2 =b.compute_unl(zero_med)
    Umat = numpy.zeros([len(zero_med)*2,whichl+1])
 
    Giwn = numpy.zeros(len(zero_med)*2)    
    
    for n in range(len(zero_med)):    
        wn = (2*zero_med[n]+1)* numpy.pi /beta 
        Giwn[n]=0#real part-pole/(wn**2 + pole**2)
        #Giwn[n]=-pole/(wn**2 + pole**2)
        Giwn[n + len(zero_med)]=-wn/(wn**2 + pole**2)#imag part          
    
    
    l2 = whichl +2
    Umat3 =beta/numpy.sqrt(2) * numpy.vstack((unl2[:][0:l2].real , unl2[:][0:l2].imag))
    
#Umat.shape
#unl2[:][:].real.shape

sol = numpy.sqrt(beta*wmax/2)*numpy.linalg.lstsq(Umat3,Giwn,rcond=None)[0]
    
#inv = numpy.linalg.inv(Umat3)

#sol = numpy.sqrt(beta*wmax/2)*inv@Giwn

print(gl2 -sol)

samplingiw= numpy.zeros(len(x))

#unl = basis.compute_unl(x)
unl =b.compute_unl(x)
for l in range(b.dim()-1):
    if l%2 == 0:
        z = numpy.sqrt(beta)*sol[l]*unl[:,l].imag
    else:
        z = numpy.sqrt(beta)*sol[l]*unl[:,l].real
    
    samplingiw = z + samplingiw



ax2 = fig.add_subplot(2,1,2)
    

ax2.plot(x_iwn, np.abs(samplingiw),label=r"$ |{\rm Im}\hat{G}^{{\rm F}}(i\omega^{{\rm F}}_n)| $")
#plt.plot(x, gim,label=r"$ -im $")
ax2.plot(x_iwn, samplingiw-gim,"--",label=r"$| \Delta {\rm Im}\hat{G}^{{\rm F}}(i\omega^{{\rm F}}_n) |$")

print(samplingiw[0],gim[0])

    



ax2.set_xscale("symlog",linthreshx=shx)
ax2.set_yscale("log")
ax2.set_xlim(-100,100)
ax2.set_ylim(1e-15,10)

#plt.scatter(iwnset,Giwn[0:len(zero_med)],marker="x",color = "r",label=r"${\rm Re}[G(i\omega_n)]$")

zero_med = x_sp*1
ax2.scatter(zero_med,np.abs(Giwn[len(zero_med):2*len(zero_med)]),marker="x",color = "r",s=74)

#ax2.scatter(zero_med,np.abs(Giwn[len(zero_med):2*len(zero_med)]),marker="x",color = "r",s=74)
#ax2.scatter(zero_med,-Giwn[len(zero_med):2*len(zero_med)],marker="x",color = "r",s=74)#,label=r"${\rm Sampling}\hspace{3mm} {\rm point}$")#,label=r"$G_l$")

ax2.set_xticks([-10**2,-10**1,-10**0,0,10**0,10**1,10**2])#,[r'$-10^{2}$',r'$0$',r'$10^{2}$'])
ax2.set_ylabel(r"$|{\rm Im}\hat{G}^{{\rm F}}(i\omega^{{\rm F}}_n)|$",fontsize=21)
ax2.tick_params(labelsize=21)#,labelbottom='off')
#plt.ylabel(r'$|G_l|$',fontsize = 21)
#plt.ylabel(r'$|G_l - G_l^{lsq}|$',fontsize = 21)
#plt.xlabel(r'$\omega_n$',fontsize = 21)
plt.xlabel(r'$\omega_n^{{\rm F}}$',fontsize = 21)
#plt.legend(frameon=False,fontsize = 21)
ax2.legend(loc='upper right',
           bbox_to_anchor=(0.4, 0.1, 0.6, .600), 
           borderaxespad=0.,frameon=False,fontsize = 20)

fig.tight_layout()
plt.savefig('subplotiw100'+'.pdf')  
#plt.savefig('Imsamplingpoints100'+'.pdf')  
#plt.show()




    
