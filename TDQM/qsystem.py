import numpy as np
from numpy import linalg as LA
import cmath
#from fields import Laser

class Hamiltonian:

        def __init__(self,nst,ener,hp):
                self.nst=nst
                self.H0 = ener
                self.W = hp
                self.Ht = np.zeros((nst,nst), dtype=complex)

        def TDHamiltonian(self,laser,t):
                ns = self.nst
                ft = laser.Amplitude(t)
                for i in range(ns):
                      for j in range(ns):
                                delE = self.H0[j]-self.H0[i]
                                ph = np.complex(np.cos(delE*t), np.sin(delE*t))
                                self.Ht[i,j] = ft*self.W[i,j]*ph

                return self.Ht
