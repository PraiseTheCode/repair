import astropy.units as au
import astropy.constants as ac
from scipy.interpolate import interp1d
import numpy as np
import pandas as pd

from orbit import orbit
from star import star

from lcrvcurve import lightcurve
from lcrvcurve import rvcurve

import ellc

import matplotlib.pyplot as plt


class binary:


    def __init__(self, orbit, primary, secondary, sbratio = 1.0, l3 = 0.0, gamma = 0.0, obsSpec = None, obsLC = None, obsRV = None,
                 nbins = 500, phase1 = None, phase2 = None, oversampling = None, width = None, gdc_option = "Claret"):

        self.orbit, self.primary, self.secondary, self.sbratio, self.l3, self.gamma = orbit, primary, secondary, sbratio, l3, gamma
        self.obsSpec, self.obsLC, self.obsRV = obsSpec, obsLC, obsRV

        self.sbratio_spec = None

        self.previous_orbit, self.previous_primary, self.previous_secondary = orbit, primary, secondary

        self.gdc_option = gdc_option

        self.nbins = nbins
        self.phase1 = phase1
        self.phase2 = phase2
        self.oversampling = oversampling
        self.width = width

        a1 = self.orbit.a/(1.+(1/self.orbit.q))
        a2 = self.orbit.a/(1.+self.orbit.q)
        K1 = self.CALCULATE_SEMI_AMPLITUDE(a1)
        K2 = self.CALCULATE_SEMI_AMPLITUDE(a2)



        if self.obsLC is not None:

            self.obsLC.CALC_PHASES(self.orbit.porb, self.orbit.t0)

            if oversampling != None:
                self.obsLC.BIN_LC_NON_UNIFORM(nbins, phase1, phase2, oversampling, width)
            else:
                self.obsLC.BIN_LC(nbins)




    def CALCULATE_SEMI_AMPLITUDE(self,a_i):
        a_i *= au.R_sun
        a_i = a_i.to('km').value
        period = self.orbit.porb * au.day
        period = period.to('second').value
        return 2.*np.pi*a_i*np.sin(np.deg2rad(self.orbit.inc_deg)) / np.sqrt(1.-self.orbit.ecc**2) / period

    def CALCULATE_MASSES(self,k1,k2):
        k1 = k1 * au.km / au.second
        k2 = k2 * au.km / au.second
        period = self.orbit.porb*au.day
        scale = 0.5/(np.pi*ac.G) * (1.-self.orbit.ecc**2)**(1.5) * (k1+k2)**2 * period

        m1 = ((scale * k2).to('Msun')).value
        m2 = ((scale * k1).to('Msun')).value
        sin3i = np.sin(np.deg2rad(self.orbit.inc_deg)) **3
        return m1/sin3i ,m2/sin3i


    def CALCULATE_SEMI_AMPLITUDE_e(self,a_i,porb,incl,ecc):
        a_i *= au.R_sun
        a_i = a_i.to('km').value
        period = porb * au.day
        period = period.to('second').value
        return 2.*np.pi*a_i*np.sin(np.deg2rad(incl)) / np.sqrt(1.-ecc**2) / period

    def CALCULATE_MASSES_e(self,k1,k2,porb,incl,ecc):
        k1 = k1 * au.km / au.second
        k2 = k2 * au.km / au.second
        period = porb*au.day
        scale = 0.5/(np.pi*ac.G) * (1.-ecc**2)**(1.5) * (k1+k2)**2 * period

        m1 = ((scale * k2).to('Msun')).value
        m2 = ((scale * k1).to('Msun')).value
        sin3i = np.sin(np.deg2rad(incl)) **3
        return m1/sin3i ,m2/sin3i


    def CALCULATE_LOGG(self,m,r):
        m *= au.Msun
        r *= au.Rsun
        g = (ac.G * m / (r**2)).to('cm / s^2')
        return np.log10(g.value)

    def CALCULATE_MASS_FROM_LOGG(self,logg,r):
        #m *= au.Msun
        r *= au.Rsun
        g = (10**(logg))*au.cm/au.s**2
        m = g*r**2/ac.G
        r /= au.Rsun
        return m.to(au.Msun) / au.Msun

    def TIMEPERI_TO_TIMETRANS(self, tp, per, ecc, omega):
        try:
            if ecc >= 1:
                return tp
        except ValueError:
            pass

        f = np.pi/2 - omega                                         # true anomaly during transit
        ee = 2 * np.arctan(np.tan(f/2) * np.sqrt((1-ecc)/(1+ecc)))  # eccentric anomaly

        tc = tp + per/(2*np.pi) * (ee - ecc*np.sin(ee))         # time of conjunction

        return tc

    def COMPUTE_MODEL_LC(self):

        times = self.orbit.t0 + np.linspace(0, self.orbit.porb, 1000)

        t0_ellc = self.TIMEPERI_TO_TIMETRANS( self.orbit.t0, self.orbit.porb, self.orbit.ecc, self.orbit.omega)


        lcc = ellc.lc(times, self.primary.radius_a, self.secondary.radius_a, self.sbratio, self.orbit.inc_deg,
                     light_3=self.l3, t_zero=t0_ellc, period=self.orbit.porb,
                     a=self.orbit.a, q=self.orbit.q, f_c=self.orbit.f_c, f_s=self.orbit.f_s,
                     rotfac_1=self.primary.rotfac, rotfac_2=self.secondary.rotfac,
                     vsini_1=self.primary.vsini, vsini_2=self.secondary.vsini, shape_1=self.primary.shape,
                     shape_2=self.secondary.shape, heat_1=self.primary.heated, heat_2=self.secondary.heated,
                     ld_1='mugrid', ldc_1=self.primary.ldc, ld_2='mugrid',
                     ldc_2=self.primary.ldc, gdc_1=self.primary.gdc, gdc_2=self.secondary.gdc,
                     bfac_1=self.primary.bfac, bfac_2=self.secondary.bfac)


        sbratio_spec = self.sbratio
        if self.sbratio_spec is not None:
            sbratio_spec = self.sbratio_spec
        primary_ldc_spec = self.primary.ldc
        if self.primary.ldc_spec is not None:
            primary_ldc_spec = self.primary.ldc_spec
        secondary_ldc_spec = self.secondary.ldc
        if self.secondary.ldc_spec is not None:
            secondary_ldc_spec = self.secondary.ldc_spec

        lcc2 = lcc
        if self.sbratio_spec is not None:
            lcc2 = ellc.lc(times, self.primary.radius_a, self.secondary.radius_a, sbratio_spec, self.orbit.inc_deg,
                     #light_3=self.l3,
                     t_zero=t0_ellc, period=self.orbit.porb,
                     a=self.orbit.a, q=self.orbit.q, f_c=self.orbit.f_c, f_s=self.orbit.f_s,
                     rotfac_1=self.primary.rotfac, rotfac_2=self.secondary.rotfac,
                     vsini_1=self.primary.vsini, vsini_2=self.secondary.vsini, shape_1=self.primary.shape,
                     shape_2=self.secondary.shape, heat_1=self.primary.heated, heat_2=self.secondary.heated,
                     ld_1='mugrid', ldc_1=primary_ldc_spec, ld_2='mugrid',
                     ldc_2=secondary_ldc_spec, gdc_1=self.primary.gdc, gdc_2=self.secondary.gdc,
                     bfac_1=self.primary.bfac, bfac_2=self.secondary.bfac)

        flux1, flux2 = ellc.fluxes(times, self.primary.radius_a, self.secondary.radius_a, sbratio_spec, self.orbit.inc_deg,
                     t_zero=t0_ellc, period=self.orbit.porb,
                     a=self.orbit.a, q=self.orbit.q, f_c=self.orbit.f_c, f_s=self.orbit.f_s,
                     rotfac_1=self.primary.rotfac, rotfac_2=self.secondary.rotfac,
                     vsini_1=self.primary.vsini, vsini_2=self.secondary.vsini, shape_1=self.primary.shape,
                     shape_2=self.secondary.shape, heat_1=self.primary.heated, heat_2=self.secondary.heated,
                     ld_1='mugrid', ldc_1=primary_ldc_spec, ld_2='mugrid',
                     ldc_2=secondary_ldc_spec, gdc_1=self.primary.gdc, gdc_2=self.secondary.gdc,
                     bfac_1=self.primary.bfac, bfac_2=self.secondary.bfac)

        self.flux_ratios_1 = flux1 / lcc2
        self.flux_ratios_2 = flux2 / lcc2


        return times, lcc


    def COMPUTE_MODEL_RV(self):
        times = np.linspace(0, self.orbit.porb, 1000)

        rv1 = self.CALCULATE_MODEL_RV_DIRECT(rv_times=times, component=1)
        rv2 = self.CALCULATE_MODEL_RV_DIRECT(rv_times=times, component=2)

        return times, rv1, rv2

    def CALCULATE_PHASE_FROM_TRUE_ANOMALY(self,T,ecc,omega,pshift):

        E = 2*np.arctan(np.sqrt((1-ecc)/(1+ecc)) * np.tan(T/2))
        M = E - ecc*np.sin(E)
        return (M+omega)/(2.0*np.pi) - 0.25 + pshift


    def CALCULATE_A1(self,a,q):
        return a/(1.+(1/q)) # solar radii

    def CALCULATE_A2(self,a,q):
        return a/(1.+q) # solar radii

    def SOLVE_ECC_OMEGA(self,fc,fs):
        #fs = sqrt(e)*sin(omega)
        #fc = sqrt(e)*cos(omega)

        omega = np.arctan2(fs,fc)
        ecc = fs**2 + fc**2

        return ecc,omega


    def GET_TRUE_ANOMALY(self,M,ecc,itermax=10):

      #-- initial value
      Fn = M + ecc*np.sin(M) + ecc**2/2.*np.sin(2*M)
      #-- iterative solving of the transcendent Kepler's equation
      for i in range(itermax):
        F = Fn
        Mn = F-ecc*np.sin(F)
        Fn = F+(M-Mn)/(1.-ecc*np.cos(F))
        keep = F!=0 #-- take care of zerodivision
        if hasattr(F,'__iter__'):
          if np.all(abs((Fn-F)[keep]/F[keep])<0.00001):
            break
        elif (abs((Fn-F)/F)<0.00001):
          break

      #-- relationship between true anomaly (theta) and eccentric anomalie (Fn)
      true_an = 2.*np.arctan(np.sqrt((1.+ecc)/(1.-ecc))*np.tan(Fn/2.))

      return Fn,true_an

    def GET_RV_CURVE_SINGLE(self,times,period,hjd0,ecc,omega,K,gamma):

      orbital_frequency = 1./period

      M_const = 2. * np.pi * orbital_frequency
      M0 = hjd0 * M_const

      M = times * M_const - M0
      E, true_anomaly = self.GET_TRUE_ANOMALY(M,ecc)

      rv = gamma + K*( ecc*np.cos(omega) + np.cos(true_anomaly+omega) )

      return rv

    def CALCULATE_MODEL_RV_DIRECT(self,rv_times=None,component=1):

        if component == 1:
            a_i = self.CALCULATE_A1(self.orbit.a,self.orbit.q)

        else:
            a_i = self.CALCULATE_A2(self.orbit.a,self.orbit.q)

        t0 = self.orbit.t0
        period = self.orbit.porb
        incl = self.orbit.inc_deg
        ecc = self.orbit.ecc
        omega = self.orbit.omega
        k_i = self.CALCULATE_SEMI_AMPLITUDE(a_i)
        if component == 1:
            rv_syn = self.GET_RV_CURVE_SINGLE(rv_times, period, t0, ecc, omega, k_i, 0.)
        else:
            rv_syn = self.GET_RV_CURVE_SINGLE(rv_times, period, t0, ecc, omega-np.pi, k_i, 0.)


        return rv_syn

    def CALCULATE_MODEL_RV_ellc(self,rv_times):

        t0_ellc = self.TIMEPERI_TO_TIMETRANS( self.orbit.t0, self.orbit.porb, self.orbit.ecc, self.orbit.omega)

        if_fluxw = True
        if self.primary.heated is not None or self.secondary.heated is not None:
            if_fluxw = False

        rv1, rv2 = ellc.rv(rv_times, self.primary.radius_a, self.secondary.radius_a, self.sbratio, self.orbit.inc_deg,
                     t_zero=t0_ellc, period=self.orbit.porb,
                     a=self.orbit.a, q=self.orbit.q, f_c=self.orbit.f_c, f_s=self.orbit.f_s,
                     rotfac_1=self.primary.rotfac, rotfac_2=self.secondary.rotfac,
                     vsini_1=self.primary.vsini, vsini_2=self.secondary.vsini, shape_1=self.primary.shape,
                     shape_2=self.secondary.shape, heat_1=self.primary.heated, heat_2=self.secondary.heated,
                     ld_1='mugrid', ldc_1=self.primary.ldc, ld_2='mugrid',
                     ldc_2=self.primary.ldc, gdc_1=self.primary.gdc, gdc_2=self.secondary.gdc,
                     bfac_1=self.primary.bfac, bfac_2=self.secondary.bfac,
                     flux_weighted=if_fluxw)

        return rv1, rv2


    def CALCULATE_VSINI(self,incl, porb,r_i,a,f_i):
        sini = np.sin(np.deg2rad(incl))
        porb *= 86400.
        r_i *= a * 695700.
        val = 2.*np.pi* r_i * f_i * sini / porb
        return val




    def UPDATE_PARAMS(self, params):

        ecc, omega = 0.0, np.pi
        if params['f_c'] != 0.0 or params['f_s'] != 0.0:
            ecc, omega = self.SOLVE_ECC_OMEGA(params['f_c'], params['f_s'])


        if params['a'] < 0:
            params['a'] = np.abs(params['a'])
        if params['r1'] < 0:
            params['r1'] = np.abs(params['r1'])
        if params['r2'] < 0:
            params['r2'] = np.abs(params['r2'])
        #if params['M1'] < 0:
        #    params['M1'] = np.abs(params['M1'])
        if params['q'] < 0:
            params['q'] = np.abs(params['q'])
        if params['Ve1sini'] < 0:
            params['Ve1sini'] = np.abs(params['Ve1sini'])
        if params['Ve2sini'] < 0:
            params['Ve2sini'] = np.abs(params['Ve2sini'])
        if params['porb'] < 0:
            params['porb'] = np.abs(params['porb'])
        if params['incl'] < 0:
            params['incl'] = np.abs(params['incl'])
        if params['incl'] > 90:
            params['incl'] = 180.0 - params['incl']

        ecc = 0.0
        try:
            ecc,omega = self.SOLVE_ECC_OMEGA(params["f_c"], params["f_s"])
        except:
            ecc,omega = 0.0, None
        a1 = params['a']/(1.+(1/params['q']))
        a2 = params['a']/(1.+params['q'])
        K1 = self.CALCULATE_SEMI_AMPLITUDE_e(a1,params['porb'],params['incl'],ecc)
        K2 = self.CALCULATE_SEMI_AMPLITUDE_e(a2,params['porb'],params['incl'],ecc)
        M1,M2 = self.CALCULATE_MASSES_e(K1, K2,params['porb'],params['incl'],ecc)
        print(K1,K2,M1)

        transformed_params = {
            "m2" : params["q"],
            "M1" : M1,
            "M2" : params["q"]*M1,
            "R1" : params["r1"]*params["a"],
            "R2" : params["r2"]*params["a"],
            "Ve1": params["Ve1sini"] / np.sin(np.deg2rad(params["incl"])),
            "Ve2": params["Ve2sini"] / np.sin(np.deg2rad(params["incl"])),
            "sbratio" : 1.0,
            "logg1" : self.CALCULATE_LOGG(M1, params["r1"]*params["a"]),
            "logg2" : self.CALCULATE_LOGG(params["q"]*M1, params["r2"]*params["a"]),
            "ecc" : ecc,
            "omega" : omega,
            }

        #print("!")
        #print(params["logg1"])
        #print(self.CALCULATE_MASS_FROM_LOGG(params["logg1"], transformed_params["R1"]))
        #print(params["q"])
        #print("!")

        a_star1= params["a"] / (1 + params["q"])
        a_star2= params["q"] * params["a"] / (1 + params["q"])

        if params["q"] > 0:
            R_lobe_1 = params["a"] * (0.49 * params["q"]**(2/3)) / (0.6 * params["q"]**(2/3) + np.log(1 + params["q"]**(1/3)))
            R_lobe_2 = params["a"] * (0.49 * (1/params["q"])**(2/3)) / (0.6 * (1/params["q"])**(2/3) + np.log(1 + (1/params["q"])**(1/3)))
        else:
            # print("Mass ratio q must be greater than 0.")
            R_lobe_1 = 1e100
            R_lobe_2 = 1e100

        R_lobe_1 *= (1-transformed_params["ecc"]-0.01)
        R_lobe_2 *= (1-transformed_params["ecc"]-0.01)

        if transformed_params["R1"] > R_lobe_1:
            print(f'!!! R1 limited from {transformed_params["R1"]} to {R_lobe_1} (R_lobe * (1-ecc) )!!!')
            print(f'(r1 updated from {params["r1"]} to {R_lobe_1 / params["a"]})')
            transformed_params["R1"] = R_lobe_1
            params["r1"] = transformed_params["R1"] / params["a"]
            transformed_params["logg1"] = self.CALCULATE_LOGG(transformed_params["M1"], transformed_params["R1"])

        if transformed_params["R2"] > R_lobe_2:
            print(f'!!! R2 limited from {transformed_params["R2"]} to {R_lobe_2} (R_lobe * (1-ecc) ) !!!')
            print(f'(r2 updated from {params["r2"]} to {R_lobe_2 / params["a"]})')
            params["r2"] = R_lobe_2 / params["a"]
            transformed_params["R2"] = R_lobe_2
            transformed_params["logg2"] = self.CALCULATE_LOGG(transformed_params["M2"], transformed_params["R2"])


        sm = 1.989e30
        sr = 6.957e8
        cG = 6.6743e-11

        vsync1 = 2 * np.pi * transformed_params["R1"] * sr *1e-3  / (params["porb"]*24*60*60)
        vsync2 = 2 * np.pi * transformed_params["R2"] * sr *1e-3 / (params["porb"]*24*60*60)


        rotfac_1 = transformed_params["Ve1"] / vsync1
        rotfac_2 = transformed_params["Ve2"] / vsync2


        breakup1 = np.sqrt(2 * cG * transformed_params["M1"] * sm / 3 / 3 /transformed_params["R1"] / sr ) / 1000.0
        breakup2 = np.sqrt(2 * cG * transformed_params["M2"] * sm / 3 / 3/ transformed_params["R2"] / sr ) / 1000.0

        #print("Ve, sync, breakup, rotfac")

        #print(transformed_params["Ve1"], vsync1, breakup1, rotfac_1)
        #print(transformed_params["Ve2"], vsync2, breakup2, rotfac_2)

        if rotfac_1*vsync1 >= breakup1:
            rotfac_1 = breakup1 / vsync1 / 1.5
            #print("new rotfac1:", rotfac_1)
        if rotfac_2*vsync2 >= breakup2:
            rotfac_2 = breakup2 / vsync2 / 1.5
            #print("new rotfac2:", rotfac_2)

        orb = orbit(porb = params["porb"], inc_deg = params["incl"], t0 = params["t0"], a = params["a"], q = params["q"],
                    f_c = params["f_c"], f_s = params["f_s"])


        primary = star(radius_a = params["r1"], radius_sol = transformed_params["R1"], mass = transformed_params["M1"], v_equatorial = transformed_params["Ve1"],
                       rotfac = rotfac_1,
                       vsini = params["Ve1sini"], Teff = params["Teff1"], heated = params["heat_1"],
                       logg = transformed_params["logg1"], metallicity=params["metallicity1"], gdc = params["gdc_1"], label = 1, orbit = orb)
        secondary = star(radius_a = params["r2"], radius_sol = transformed_params["R2"], mass = transformed_params["M2"], v_equatorial = transformed_params["Ve2"],
                       rotfac = rotfac_2,
                       vsini = params["Ve2sini"], Teff = params["Teff2"], heated = params["heat_2"],
                       logg = transformed_params["logg2"], metallicity=params["metallicity2"], gdc = params["gdc_2"], label = 2, orbit = orb)

        secondary.radius_R1 = params["r2"]/params["r1"]
        secondary.mass_M1 = transformed_params["M2"]/transformed_params["M1"]

        primary.sync_fact = rotfac_1
        secondary.sync_fact = rotfac_2

        primary.ve_sync = vsync1
        secondary.ve_sync = vsync2

        self.previous_orbit, self.previous_primary, self.previous_secondary = self.orbit, self.primary, self.secondary

        self.orbit = orb
        self.primary = primary
        self.secondary = secondary

        self.sbratio = transformed_params["sbratio"]
        self.l3 = params["l3"]
        self.gamma = params["gamma"]



        if self.orbit.t0 != self.previous_orbit.t0 or self.orbit.porb != self.previous_orbit.porb :
            if self.obsLC is not None:
                self.obsLC.CALC_PHASES(self.orbit.porb, self.orbit.t0)

                if self.oversampling != None:
                    self.obsLC.BIN_LC_NON_UNIFORM(self.nbins, self.phase1, self.phase2, self.oversampling, self.width)
                else:
                    self.obsLC.BIN_LC(self.nbins)

            if self.obsSpec is not None:
                self.obsSpec.CALC_PHASES(self.orbit.porb, self.orbit.t0)



    def UPDATE_MODEL_LC(self):

        times, lcc = self.COMPUTE_MODEL_LC()

        mod_lc = lightcurve(times, lcc, np.zeros_like(lcc))
        mod_lc.CALC_PHASES(self.orbit.porb, self.orbit.t0)

        phase = mod_lc.phase
        flux_in_phase = mod_lc.flux_in_phase


        ellcrv1, ellcrv2 = self.CALCULATE_MODEL_RV_ellc(times)


        mod_rv1 = rvcurve(times, ellcrv1, None)
        mod_rv1.CALC_PHASES(self.orbit.porb, self.orbit.t0)

        mod_rv2 = rvcurve(times, ellcrv2, None)
        mod_rv2.CALC_PHASES(self.orbit.porb, self.orbit.t0)


        ii = interp1d(phase, flux_in_phase, fill_value="extrapolate")
        self.modLC = ii(self.obsLC.binned_phase)
        self.modRV1 = self.obsSpec.INTERP_RVS(self.orbit.porb, self.orbit.t0, mod_rv1.phase, mod_rv1.rv_in_phase)
        self.modRV2 = self.obsSpec.INTERP_RVS(self.orbit.porb, self.orbit.t0, mod_rv2.phase, mod_rv2.rv_in_phase)
        self.modRV1_full = mod_rv1.rv_in_phase
        self.modRV2_full = mod_rv2.rv_in_phase
        self.modRV1_full_phase = mod_rv1.phase
        self.modRV2_full_phase = mod_rv2.phase

        self.residLC = (self.obsLC.binned_flux - self.modLC)** 2/ self.modLC
        self.residLCNoC = (self.obsLC.binned_flux - self.modLC)** 2



    def UPDATE_SP_SEPARATION(self, ifplot = False, if_lightfact = False):

        rvs = np.zeros(shape=(2,self.obsSpec.nspec))
        rvs[0] = self.modRV1 + self.gamma
        rvs[1] = self.modRV2 + self.gamma


        if if_lightfact:
            ffi1 = interp1d(np.linspace(0, 1, 1000), self.flux_ratios_1)
            ffi2 = interp1d(np.linspace(0, 1, 1000), self.flux_ratios_2)
            lfs = np.ones(shape=(2,self.obsSpec.nspec))


            for i in range(self.obsSpec.nspec):
                lfs[0,i] = ffi1(self.obsSpec.phases[i])
                lfs[1,i] = ffi2(self.obsSpec.phases[i])

                self.obsSpec.SET_LIGHTFACT(lfs)

        self.obsSpec.SEPARATE(rvs)

        self.wlstep=self.obsSpec.wl_eqlog[1]-self.obsSpec.wl_eqlog[0]


        ntrim = 0

        try:
            ntrim = int(max(self.obsSpec.rvs_bin[0])*2)
            ntrim = np.abs(ntrim)
            residuals = np.empty((len(self.obsSpec.wl_eqlog) - 2*ntrim, self.obsSpec.nspec))
        except:
            ntrim = 0
            residuals = np.empty((len(self.obsSpec.wl_eqlog) - 2*ntrim, self.obsSpec.nspec))


        wls1_all = self.obsSpec.wl_eqlog[:, np.newaxis] + self.obsSpec.rvs_bin[0] * self.wlstep
        wls2_all = self.obsSpec.wl_eqlog[:, np.newaxis] + self.obsSpec.rvs_bin[1] * self.wlstep

        if ifplot:
            plt.plot(self.obsSpec.wl_eqlog, self.obsSpec.fluxes_eqlog[:,0],"-k",lw=0.5)
            plt.plot(self.obsSpec.wl_eqlog, self.obsSpec.mod[0],lw=0.5)
            plt.plot(self.obsSpec.wl_eqlog, self.obsSpec.mod[1],lw=0.5)
            plt.show()


        self.mod1_var =  np.mean(self.obsSpec.fluxes_wl_vars) / np.mean(lfs[0])
        self.mod2_var =  np.mean(self.obsSpec.fluxes_wl_vars) / np.mean(lfs[1])

        #print(f"var LC: {self.obsLC.flux_var}")
        #print(f"var specs: {self.obsSpec.fluxes_wl_vars[0]}, {self.obsSpec.fluxes_wl_vars[1]}, {self.obsSpec.fluxes_wl_vars[2]}, ...")
        #print(f"var sep opt 2: {self.mod1_var}, {self.mod2_var}")

        wls1_all_l = []
        wls2_all_l = []
        for i in range(len(rvs[0])):
            aa = np.e**self.obsSpec.wl_eqlog * (1 + rvs[0][i]/299792.45)
            wls1_all_l.append(aa)
        for i in range(len(rvs[1])):
            aa = np.e**self.obsSpec.wl_eqlog * (1 + rvs[1][i]/299792.45)
            wls2_all_l.append(aa)

        self.summed_mod = np.empty((len(self.obsSpec.wl_eqlog) - 2*ntrim, self.obsSpec.nspec))

        residualsSqNorm = residuals.copy()
        residualsNoSqNorm = residuals.copy()
        for i in range(self.obsSpec.nspec):
            wls1 = wls1_all[:, i]
            wls2 = wls2_all[:, i]


            ii = interp1d(wls1[ntrim:-ntrim], self.obsSpec.mod[0, ntrim:-ntrim], bounds_error=False, fill_value=1.0)
            fl1 = ii(self.obsSpec.wl_eqlog)

            ii = interp1d(wls2[ntrim:-ntrim], self.obsSpec.mod[1, ntrim:-ntrim], bounds_error=False, fill_value=1.0)
            fl2 = ii(self.obsSpec.wl_eqlog)

            summed_flux = fl1*self.obsSpec.lfs[0,i] + fl2*self.obsSpec.lfs[1,i]

            residuals[:, i] = self.obsSpec.fluxes_eqlog[ntrim:-ntrim, i] - summed_flux[ntrim:-ntrim]
            residualsSqNorm[:, i] = ((self.obsSpec.fluxes_eqlog[ntrim:-ntrim, i] - summed_flux[ntrim:-ntrim])**2)/np.abs(summed_flux[ntrim:-ntrim])
            self.summed_mod[:, i] = summed_flux[ntrim:-ntrim]

            residualsNoSqNorm[:, i] = (self.obsSpec.fluxes_eqlog[ntrim:-ntrim, i] - summed_flux[ntrim:-ntrim])/(summed_flux[ntrim:-ntrim])

            if ifplot:
                plt.plot(self.obsSpec.wl_eqlog[ntrim:-ntrim], self.obsSpec.fluxes_eqlog[ntrim:-ntrim, i], "-k",lw=0.5)
                plt.plot(self.obsSpec.wl_eqlog[ntrim:-ntrim], summed_flux[ntrim:-ntrim], "-r", lw=0.5)
                plt.show()


        self.residSpec = residuals
        self.residSpecSqNorm = residualsSqNorm

        self.residSpecNoSqNorm = residualsNoSqNorm
