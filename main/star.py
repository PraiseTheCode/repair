class star:
    def __init__(self, radius_a = None, radius_sol = None, mass = None, rotfac = 1, v_equatorial = None, ve_sync = None, 
                 vsini = None, heated = None, shape = "roche", ldc = None, gdc = None, label = None, bfac = None, 
                 Teff = None, logg = None, metallicity = None, orbit = None):
        
        self.radius_a, self.radius_sol, self.mass, self.rotfac, self.v_equatorial, self.vsini, self.heated, self.shape, self.ldc, self.gdc, self.label, self.bfac, self.Teff, self.logg, self.metallicity, self.orbit = radius_a, radius_sol, mass, rotfac, v_equatorial, vsini, heated, shape, ldc, gdc, label, bfac, Teff, logg, metallicity, orbit 
        
        self.radius_R1, self.mass_M1 = None, None
        
        self.sync_fact = None
        
        self.spectrum = None
        self.spectrum_calib = None
        
        self.ldc = None
        self.ldc_spec = None
