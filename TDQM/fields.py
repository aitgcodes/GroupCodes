import numpy as np
import cmath

class Laser:

        def __init__(self,tlaser,f0,w0,fpol):
                self.f0=f0
                self.w0=w0
                self.fpol=fpol
                self.envlp = "none"
                self.gauss_width = 0.0
                self.gauss_t0 = 0.0
                self.tlaser = tlaser

        def InitEnvelope(self,envlp,params):
                self.envlp = envlp
                if envlp == "gauss":
                        self.gauss_t0 = float(params[0])
                        self.gauss_width = float(params[1])

        def Envelope(self,t,envlp):
                fevlp = 1.0
                if envlp == "gauss":
                        norm = 1.0 #/(np.sqrt(2.0*np.pi)*self.gauss_width)
                        fevlp = norm*np.exp(-(t-self.gauss_t0)**2/2.0/self.gauss_width**2)

                return fevlp

        def Amplitude(self,t):
                ft = 0.0
                if self.tlaser:
                        fvlp = self.Envelope(t,self.envlp)
                        fosc = self.f0*np.sin(self.w0*t)
                        ft = fvlp*fosc

                return ft
