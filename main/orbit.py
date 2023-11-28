import numpy as np

class orbit:
    def __init__(self, porb = None, inc_deg = None, t0 = None, a = None, q = None, f_c = None, f_s = None) :
        
        self.porb, self.inc_deg, self.t0, self.a, self.q, self.f_c, self.f_s = porb, inc_deg, t0, a, q, f_c, f_s
        
        self.ecc, self.omega = None, None
        if self.f_c is not None:
            if self.f_c != 0 or self.f_s != 0:
                self.ecc,self.omega = self.SOLVE_ECC_OMEGA()
            else:
                self.ecc,self.omega = 0.,np.deg2rad(90.)
        else:
            self.ecc,self.omega = 0.,np.deg2rad(90.)
            
        
    def SOLVE_ECC_OMEGA(self):
        
        omega = np.arctan2(self.f_s,self.f_c)
        ecc = self.f_s**2 + self.f_c**2
    
        return ecc,omega