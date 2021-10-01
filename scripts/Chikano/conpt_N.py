from __future__ import print_function

import numpy
#import irbasis#irlib
import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from scipy import stats
#from mpmath import *
#from mpmath.calculus.quadrature import GaussLegendre
def legesort(ll1,vc):
    for i in range(len(ll1)):
        for j in range(len(ll1)):
            if ll1[i] < ll1[j]:
                a1 = ll1[i]
                ll1[i] = ll1[j]
                ll1[j] = a1
                
def metal_sp(x,omega,stat,beta):
    if stat == "F":
            if omega < 0 :
                return -numpy.sqrt(1-omega**2)*numpy.exp((beta-x)*omega)/(1+numpy.exp(beta*omega))
                    
                    
            else:
                return -numpy.sqrt(1-omega**2)*numpy.exp(-x*omega)/(1+numpy.exp(-beta*omega))
            
    elif stat == "B":
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
    
    if stat =="F":
        return 2*gt[0]/numpy.pi #0.392699 
    elif stat =="B":
        return 3*gt[0]/2


def leg(l,x):
    c = numpy.zeros(l+1)
    
    c[l]=1
    #return numpy.polynomial.legendre.legval(x,c)
    return scipy.special.eval_legendre(l,x)
def cheu(l,x):
    #c = numpy.zeros(l+1)
    
    #c[l]=1
    #return numpy.polynomial.chebyshev.chebval(x,c)
    return scipy.special.eval_chebyt(l, x)

    


def safe_exp(x,beta,stat):
    #insulating

    if stat == "F":
        #if beta-x < 30:
        
        return -(numpy.exp(x-beta))/(1+numpy.exp(-beta))
        #else:
        #    return 0
    elif stat == "B":
        return (numpy.exp(x-beta))/(numpy.exp(-beta)-1)              
    
    
    
    
    
def power_mesh2(x,y,div,coef):
    s = []
    z = y -x
    for i in range(div):
        s.append(x + z*(coef)**i)
        s.append(y - z*(coef)**i)
    s.sort()
    
    return s

#print(metal(100,100,"b"))
        
def goftau(f,listl,ll,vv,beta,Lambda,cutoff,basis,uni,mat,kind):
    swl= 0
    lcut = 0
    listtau = []
    glcmaxtau = []
    for i in range(len(listl)):
        glcmaxtau.append(beta)
    if kind == "u":
        #for i in numpy.linspace(0,4,3):#[0.5,1,1.2,1.4,1.5,1.6,1.8,2,2.2,2.5,2.8,3,3.3,3.5,4,5,6]:
            #listtau.append((beta-i/beta)/beta)
        #    listtau.append(1-(10)**(-i))
        for i in []:#numpy.linspace(6,12,2):
            listtau.append(1-(10)**(-i))
        listtau.append(1)
    elif kind == "legendre":
        #listtau.append(0.5)
        listtau.append(1)
        
    for taro in listtau:
        gtau = 0
        glctau = 0
        
        goftaul =[]
        game = 0
        
        for l in listl:#[b.dim()-2,b.dim()-1]:                
            l = int(l)
            """
            if mat =="insu":
                if stat == 'b' and l%2==1:
                    
                    continue
                elif stat == 'f' and l%2 ==1:    
                    continue
            """
            if mat =="metal":
                if stat == 'B' and l%2==1:
                    
                    continue
                elif stat == 'F' and l%2 ==1:
                    
                    continue
            
            
            if uni == 0:
                un = numpy.sqrt((2*l+1)/2)
                c = scipy.integrate.quad(lambda tau:  numpy.sqrt(2)*un   *f(tau,beta,stat)*basis(l,2*tau/beta -1) ,0,beta,limit = 600,points=[beta/100000,beta/10000,beta/1000,beta/100,beta/10,beta/1.1,beta/1.001])
                
            elif uni == 1:
                un = 1
                c = scipy.integrate.quad(lambda tau:  numpy.sqrt(2)*un   *f(tau,beta,stat)*basis(l,2*tau/beta -1) ,0,beta,limit = 600,points=[beta/100000,beta/10000,beta/1000,beta/2,beta/1.1,beta/1.001,beta/1.0001,beta/1.00001])
            else:
                un = numpy.sqrt(numpy.pi/2)
                
                c = scipy.integrate.quad(lambda tau:  1/np.sqrt((1-(2*tau/beta -1)**2)) *numpy.sqrt(2)*un *f(tau,beta,stat)*basis(l,2*tau/beta -1) ,0+1,beta-1 ,limit=2600,points=[0.0000005,0.000001,0.00001,0.0001,0.01,beta/10,beta/2,beta/1.001,beta/1.0001,beta/1.00001,beta/1.000001,beta/1.00000005])
                
            
            
            
            
            goftaul.append(numpy.sqrt(2)/(beta/un)*c[0]*basis(l,taro))
            
            gtau += numpy.sqrt(2)/(beta/un) *c[0]*basis(l,taro)
            #if l%100 == 0:
            #    print(l,c[0])
            #print(numpy.sqrt(2)/(beta/un)*c[0]*basis(l,taro))
            
            #print(gtau,c[0])
            
            #vv.append(numpy.abs((-1)**sign /2 - gtau))
            #vv.append(numpy.abs(c[0]))
            #if numpy.abs(c[0]) < cutoff/1000:
            #    break
            
            if beta > 1E+4:
                if numpy.abs(gtau -(-1))< cut :
                    game += 1
                    if game >2.5:
                        ll.append(l)
                        break
            else:
                if numpy.abs(c[0])< cut :
                    game += 1
                    if game >2.5:
                        ll.append(l)
                        break
            """        
            if numpy.abs(c[0])< cut :
                game += 1
                if game >2.5:
                    ll.append(l)
                    break
            """
        
                
    return ll          
            



plt.rc('text', usetex=True)

N = 1000
xvec = numpy.linspace(0,1, N)




## Construct basis
alp = [0.5,0.5,0.5,0.5]


cutmode = 1
cut = 10**-8
#mat ="metal"
mat = "insu" 
legendreon = 1

save = 0
clear = 1
for stat in ["F"]:
    if stat == "F":
        toukei = "Fermion"
        sign = 0
    else:
        toukei = "Boson"
        sign = 1
    
    if clear ==1:
        
        ll1 = []
        vc = []
    
        vv = []
        v10 = []
        
        ll3 = []
    
    for ww in [1.0]:
        #if cutmode ==1:
        lambs = [] 
        cuts=[]
        vc = []
        irvc = []
        ll1=[]
        ll2 = []
        for beta in[10.0,100.0,1000.0,10000.0,100000.0,1000000.0]:#,10000.0,100000.0]:
            Lambda = ww*beta

            
            
            #print("Loading basis functions...")
            
            print(Lambda)
        
            
            be = int(beta)
            
            #print("dim = ", basis.dim)
            ll=[]
            plt.figure(2)
            
            listle = numpy.linspace(0,5000,5001)
            
            if legendreon == 1:
                if mat =="insu":
                    ll1=goftau(safe_exp,listle,ll1,vc,beta,Lambda,cut,leg,0,mat,"legendre")
                    
                elif mat == "metal":
                    ll1=goftau(metal,listle,ll1,vc,beta,Lambda,cut,leg,0,mat,"legendre")
                #vc = numpy.append(vc,beta)
                #goftau(safe_exp,listle,ll1,vc,beta,Lambda,cut,cheu,2,mat,"legendre")
                
                vc = numpy.append(vc,beta)  
print(ll1)
