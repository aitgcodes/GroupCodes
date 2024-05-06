import numpy as np
from numpy import linalg as LA
import cmath
#from qsystem import Hamiltonian

class Propagator:

        def __init__(self,t0,nt,dt,swexp,swprop):
                self.t0=t0
                self.nt=nt
                self.time=np.zeros(nt, dtype=float)
                self.delt=dt
                self.swexp=swexp
                self.swprop=swprop
                for it in range(self.nt):
                      self.time[it] = self.t0 + it*self.delt	

        def TDExpMethod(self,a,dt):
                (ms, ns) = a.shape
                oper = np.zeros((ns,ns), dtype = 'complex')
                sw = self.swexp
                if sw == "power":
                        fact = [1, 1, 2, 6, 24, 120, 720]
                        nmax = 4
                        for i in range(nmax+1):
                              oper += LA.matrix_power(a,i)*dt**i/fact[i]

       	        return oper

        def TDPropagate(self,ham,laser,it):
                ns = ham.nst
                prg = np.zeros((ns,ns), dtype = 'complex')
                dt = self.delt
                sw = self.swprop
                t = self.time[it]
                zimg = np.complex(0.0,1.0)
                if sw == "exp_mid":
                        tmid = t + dt/2.0
                        a = -zimg*ham.TDHamiltonian(laser,tmid)
                        prg = self.TDExpMethod(a,dt)

                elif sw == "etrs":
                        tp = t + dt	
                        a = -zimg*ham.TDHamiltonian(laser,t)
                        ap = -zimg*ham.TDHamiltonian(laser,tp)
                        Up = self.TDExpMethod(ap,dt/2.0)
                        U = self.TDExpMethod(a,dt/2.0)
                        prg = np.matmul(Up,U)

                elif sw == "exp_sim":
                        a = -zimg*ham.TDHamiltonian(laser,t)
                        prg = self.TDExpMethod(a,dt)

                return prg
