import numpy as np

beta = 300.0
U = 1.0

giv = lambda v: 0.5*( 1/(1J*v*np.pi/beta - 0.5*U) + 1/(1J*v*np.pi/beta + 0.5*U) )
gtau_single_pole = lambda tau, epsilon: -np.exp(-tau*epsilon)/(1+np.exp(-beta*epsilon))
gtau = lambda taus: 0.5*(gtau_single_pole(taus, 0.5*U) + gtau_single_pole(taus, -0.5*U))