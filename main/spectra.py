import numpy as np
from scipy.interpolate import interp1d
from scipy.linalg import svd
from scipy.optimize import curve_fit
import fnmatch
import os
import subprocess
import shutil
from astropy.io import fits
from scipy.fft import fft

class specSeries:
    def __init__(self, path = None, nspec = None, bjd = None, wls = None, fluxes = None, phases = None):
        
        self.path, self.nspec, self.bjd, self.wls, self.fluxes_wl, self.phases = path, nspec, bjd, wls, fluxes, phases
        
        self.wl_eqlog, self.fluxes_eqlog = None, None
        self.rvs_kms, self.rvs_bin = None, None
        self.lfs, self.sig = None, None
        
        self.k, self.m, self.n = 2, self.nspec, None
        
        self.dftobs, self.dftmod, self.dftres = None, None, None
        self.mod, self.resid = None, None
        
    def LOAD_HERMES(self, path):
        
        files = fnmatch.filter(os.listdir(path), '*_cf.fits')

        bjd_ar = []
        for nrm_file in files:
            fits_file = nrm_file.replace('.ndat', '.fits')
            with fits.open(path+fits_file) as hdul:
                bjd_ar.append([nrm_file,  hdul[0].header['BJD']])
             
        bjd_ar=np.array(bjd_ar)  

        bjd_ar = bjd_ar[np.argsort(bjd_ar[:, 1])]
        self.bjd = np.array( [ float(bjd_ar[i][1]) for i in range(len(bjd_ar)) ] )

        files = [ str(bjd_ar[i,0]).replace('.fits', '.ndat') for i in range(len(bjd_ar))]


        spectra =  [np.loadtxt(path+f) for f in files]
        
        wls = []
        fluxes = []
        for i in range(len(files)):
            wls.append(np.array(spectra[i][:,0]))
            fluxes.append(np.array(spectra[i][:,1]))
        self.wls = np.array(wls)
        self.fluxes_wl = np.array(fluxes)
        
        self.path = path 
        self.nspec = len(spectra)
        
        
        self.WL_TO_EQLOGWL()
        
        self.k, self.m, self.n = 2, self.nspec, len(self.wl_eqlog)
        
        self.lfs = np.ones(shape=(2,self.m))
        self.sig = np.ones(shape = (self.m))
        
        for i in range(len(self.lfs)):
            for j in range(len(self.lfs[i])):
                self.lfs[i][j]/=2.0
        
        self.dftobs = self.TO_FREQ(self.fluxes_eqlog)
        
        
        
    def LOAD_TESTSUIT(self, path, period, t0):
        files = fnmatch.filter(os.listdir(path), 'comb_sn*')

        bjd_ar = []
        for nrm_file in files:
            bjd_ar.append([nrm_file,  t0 + float(nrm_file[14:-4])*period] )
             
        bjd_ar=np.array(bjd_ar)  

        bjd_ar = bjd_ar[np.argsort(bjd_ar[:, 1])]
        self.bjd = np.array( [ float(bjd_ar[i][1]) for i in range(len(bjd_ar)) ] )
        
        #print(self.bjd)
        
        #import matplotlib.pyplot as plt
        #plt.plot(self.bjd, np.ones_like(self.bjd),"ok")
        #plt.show()
        
        #print(files)

        files = [ str(bjd_ar[i,0]) for i in range(len(bjd_ar))]
        
        #print(files)

        spectra =  [np.loadtxt(path+f) for f in files]
        
        self.CALC_PHASES(0.8, 2.0)
        #print(self.phases)
        
        wls = []
        fluxes = []
        for i in range(len(files)):
            wls.append(np.array(spectra[i][:,0]))
            fluxes.append(np.array(spectra[i][:,1]))
        self.wls = np.array(wls)
        self.fluxes_wl = np.array(fluxes)
        
        self.path = path 
        self.nspec = len(spectra)
        
        
        self.WL_TO_EQLOGWL()
        
        self.k, self.m, self.n = 2, self.nspec, len(self.wl_eqlog)
        
        self.lfs = np.ones(shape=(2,self.m))
        self.sig = np.ones(shape = (self.m))
        
        for i in range(len(self.lfs)):
            for j in range(len(self.lfs[i])):
                self.lfs[i][j]/=2.0
        
        self.dftobs = self.TO_FREQ(self.fluxes_eqlog)
        
        
    def CALC_PHASES(self, porb, t_zero):
        
        phase = ((self.bjd - t_zero) / porb) % 1
        phase[phase > 1] -= 1
        phase[phase < 0] += 1
        
        self.phases = phase
        
    
    def INTERP_RVS(self, porb, t_zero, rvs_phase, rvs):
        
        self.CALC_PHASES(porb, t_zero)
        interpolator = interp1d(rvs_phase, rvs, fill_value='extrapolate')
        rv = interpolator(self.phases)
        return rv
        
        
    def SEPARATE(self, rvs):
        
        self.SET_RV(rvs)
        
        self.dftmod, self.dftres = self.FD3SEP_COMPLEX()
        self.mod = self.TO_WL(self.dftmod)
        self.resid = self.TO_WL(self.dftres)
        
        self.EQLOGWL_TO_WL()
    
    
    def WL_TO_EQLOGWL(self):
        
        log_wavelength_min = np.max([np.log(wwl.min()) for wwl in self.wls])
        log_wavelength_max = np.min([np.log(wwl.max()) for wwl in self.wls])
        
        step = np.mean(np.diff(self.wls[0])) # / 1.5 
        self.wlstep = step
        nnn = int((np.e**log_wavelength_max-np.e**log_wavelength_min)/step) # Number of bins for resampling
        
        log_wavelength_step = (log_wavelength_max - log_wavelength_min) / nnn
        
        common_log_wavelength = np.arange(log_wavelength_min, log_wavelength_max, log_wavelength_step)
        
        self.fluxes_eqlog = np.zeros((nnn, len(self.fluxes_wl)))
        
        for i, spec in enumerate(self.fluxes_wl):
            log_wavelength = np.log(self.wls[i])
            f = interp1d(log_wavelength, spec, fill_value="extrapolate")
            self.fluxes_eqlog[:, i] = f(common_log_wavelength)
        
        self.wl_eqlog = common_log_wavelength
        
        
    
    def EQLOGWL_TO_WL(self):
        
        self.wl_eqlog_to_wl = np.e**self.wl_eqlog
        
        
    def TRIM(self, wl1, wl2):
        start_wavelength = np.log(wl1)
        end_wavelength = np.log(wl2)
        
        start_index = (np.abs(self.wl_eqlog - start_wavelength)).argmin()
        end_index = (np.abs(self.wl_eqlog - end_wavelength)).argmin()
        
        self.wl_eqlog = self.wl_eqlog[start_index:end_index + 1]
        self.fluxes_eqlog = self.fluxes_eqlog[start_index:end_index + 1, :]
        
        self.dftobs = self.TO_FREQ(self.fluxes_eqlog)
        
        self.n = len(self.wl_eqlog)
    
    def TRIM_LOG(self, logwl1, logwl2):
        start_wavelength = logwl1
        end_wavelength = logwl2
        
        start_index = (np.abs(self.wl_eqlog - start_wavelength)).argmin()
        end_index = (np.abs(self.wl_eqlog - end_wavelength)).argmin()
        
        self.wl_eqlog = self.wl_eqlog[start_index:end_index + 1]
        self.fluxes_eqlog = self.fluxes_eqlog[start_index:end_index + 1, :]
        
        self.dftobs = self.TO_FREQ(self.fluxes_eqlog)
        
        self.n = len(self.wl_eqlog)
        
    def SET_RV(self, rvs_kms):
        self.rvs_kms = rvs_kms
        self.KMS_TO_BINS(self.rvs_kms)
        
    def KMS_TO_BINS(self, rvs):
        rvstep = 299792.458 * ( - 1 + np.exp((self.wl_eqlog[-1]-self.wl_eqlog[0]) / float(len(self.wl_eqlog))) )
        
        self.rvs_bin = rvs / rvstep
        
    def SET_LIGHTFACT(self, lfm):
        self.lfs = lfm
        
    def SET_WEIGHTS(self, weights):
        self.sig = weights
        
    def TO_WL(self, array2d):
        N = self.n
        reconstructed_signals = np.zeros((len(array2d), N))
        for k in range(len(array2d)):
            reconstructed_signals[k, :] = np.fft.irfft(array2d[k, :], n=N)
        return reconstructed_signals
    
    def TO_FREQ(self, array2d):
        array2dT = array2d.T
    
        array2dT_dft = np.fft.rfft(array2dT, axis=1) 
        
        return array2dT_dft
    
    
    def dft_fwd(self, m, n, mxin):
        # Create an output array of the same shape as the input, but complex type
        mxout = np.empty((m, n), dtype=complex)
    
        # Normalization factor
        a = 1.0 / np.sqrt(n)
    
        for j in range(m):
            # Apply the normalization factor and compute the DFT
            mxout[j] = a * fft(mxin[j])
    
            # Adjust the zero-frequency component
            mxout[j][0] = mxout[j][1].real
            mxout[j][1] = 0
    
            # If n is even, set the Nyquist frequency to zero
            if n % 2 == 0:
                mxout[j][n // 2] = 0
    
        return mxout

    def dft_bck(self, m, n, mxin):
        # Create the output array with the correct shape for time-domain signals
        mxout = np.zeros((m, n))
        a = 1.0 / np.sqrt(n)
    
        for j in range(m):
            # Perform the actual inverse DFT here using your custom implementation
            # Placeholder using numpy's irfft for demonstration:
            mxout[j] = np.fft.irfft(mxin[j] / a, n=n)
    
        return mxout
    
    def TO_FREQ_b(self, array2d):
        M, N = array2d.shape
        array2dT = array2d.T
        
        array2dT_dft = self.dft_fwd(N, M, array2dT)
        
        return array2dT_dft
    
    def TO_WL_b(self, array2d):
        K, N = array2d.shape
        reconstructed_signals = self.dft_bck(K, self.n, array2d)
        return reconstructed_signals
    
    
    
    def FD3SEP_COMPLEX(self):
        K, M, N = self.k, self.m, self.n
        #print(K,M,N)
        lfm, rvm, sig = self.lfs, self.rvs_bin, self.sig
        #print(lfm.shape)
        #print(rvm.shape)
        #print(sig.shape)
        
        #print(self.dftobs.shape)
        
        
        A = np.zeros((2 * M, 2 * K))
        U = np.zeros_like(A)
        V = np.zeros((2 * K, 2 * K))
        S = np.zeros(2 * K)
        b = np.zeros(2 * M)
        x = np.zeros(2 * K)
        dftmod = np.zeros((K, N // 2 + 1), dtype=complex)
        dftres = np.zeros((M, N // 2 + 1), dtype=complex)
    
        for n in range(N // 2 + 1):
            q = 2.0 * np.pi * n / N
            for j in range(M):
                s = sig[j]
                b[2 * j] = self.dftobs[j, n].real / s
                b[2 * j + 1] = self.dftobs[j, n].imag / s
    
                for k in range(K):
                    v = rvm[k, j]
                    fv = np.floor(v)
                    rez = lfm[k, j] * (((fv + 1) - v) * np.cos(fv * q) + (v - fv) * np.cos((fv + 1) * q))
                    imz = -lfm[k, j] * (((fv + 1) - v) * np.sin(fv * q) + (v - fv) * np.sin((fv + 1) * q))
                    A[2 * j, 2 * k] = rez / s
                    A[2 * j, 2 * k + 1] = -imz / s
                    A[2 * j + 1, 2 * k] = imz / s
                    A[2 * j + 1, 2 * k + 1] = rez / s
    
            U, S, Vt = svd(A, full_matrices=False)
            V = Vt.T
    
            # Solve for model parameters
            # Avoiding small singular values
            cutoff = 1e-9 * S[0]
            S_inv = np.zeros_like(S)
            for i in range(len(S)):
                if S[i] > cutoff:
                    S_inv[i] = 1. / S[i]
    
            x = V @ (np.diag(S_inv) @ (U.T @ b))
    
            for k in range(K):
                dftmod[k, n] = x[2 * k] + 1j * x[2 * k + 1]
    
            # Compute s2 and dftres
            for j in range(M):
                s = sig[j]
                bc_real = bc_imag = 0.0
                for k in range(K):
                    bc_real += (A[2 * j, 2 * k] * x[2 * k] - A[2 * j, 2 * k + 1] * x[2 * k + 1])
                    bc_imag += (A[2 * j + 1, 2 * k] * x[2 * k] + A[2 * j + 1, 2 * k + 1] * x[2 * k + 1])
                dftres[j, n] = (b[2 * j] - bc_real + 1j * (b[2 * j + 1] - bc_imag)) * s
    
        return dftmod, dftres
        


class synthetic:
    
    def __init__(self, wd_synthv = "synthV_Imu/", config_synthv = "SynthV.config", exec_synthv = "SynthV", 
                 if_convolve = False, wd_convolve = "convolve/", config_convolve = "convolve.conf", exec_convolve = "convolve", 
                 db_atmos_models = "Kurucz/", 
                 if_imu = False, if_lines=True,
                 pref = ""):
        
        self.wd_synthv = wd_synthv
        self.config_synthv = config_synthv
        self.exec_synthv = exec_synthv
        
        self.if_convolve = if_convolve
        self.wd_convolve = wd_convolve
        self.config_convolve = config_convolve
        self.exec_convolve = exec_convolve
        
        self.db_atmos_models = db_atmos_models
        
        self.if_imu = if_imu
        self.if_lines = if_lines
    
        self.pref = pref
        
        self.grid_logg = None
    
    def FIND_NEAREST_MODEL(self, Teff, logg, met): 
        
        atm_path = self.db_atmos_models
        dirs_metdir = [d for d in os.listdir(atm_path) if os.path.isdir(os.path.join(atm_path, d))]
        
        
        
        pref = "lp"
        if met < 0:
            pref = "lm"
        
        met_dir = ""
        diff_met = 1e5
        
        for d in dirs_metdir:
            if pref in d:
                loc_diff_met = abs(abs(met) - float(d[2:4]) / 10.0)
                if loc_diff_met < diff_met:
                    diff_met = loc_diff_met
                    met_dir = d
    
        files_atm = [f for f in os.listdir(os.path.join(atm_path, met_dir)) if f.endswith(".mod")]
        
        nearest = ""
        diff_teff = 1e5
        diff_logg = 1e5
        
        indTeff = 7
        lenTeff = 5
        indlogg = 13
        lenlogg = 4
        factTeff = 1.0
        factlogg = 100.0
        
        for f in files_atm:
            loc_diff_teff = abs(Teff - float(f[indTeff:indTeff + lenTeff]) / factTeff)
            loc_diff_logg = abs(logg - float(f[indlogg:indlogg + lenlogg]) / factlogg)
    
            if loc_diff_teff < diff_teff and loc_diff_logg <= diff_logg:
                diff_teff = loc_diff_teff
                diff_logg = loc_diff_logg
                nearest = f
                
                self.grid_logg = float(f[indlogg:indlogg + lenlogg]) / factlogg
                
    
        #print(nearest)
    
        src = os.path.join(atm_path, met_dir, nearest)
        dst = os.path.join(self.pref, self.wd_synthv, nearest)
        shutil.copy(src, dst)
    
        return nearest

    
        
    def REWRITE_CONFIG_SYNTHV(self, atmModel, vmicro, wl1, wl2, wlstep):
        
        with open(os.path.join(self.pref, self.wd_synthv, self.config_synthv), 'w') as f:
            if self.if_lines:
                f.write("VALD2012.lns" + "\n")
            else:
                f.write("skip ! VALD2012.lns" + "\n")
            f.write("skip molec_sel.lns" + "\n")
            f.write(atmModel + "\n")
            f.write(str(wl1) + "\t" + str(wl2) + "\t" + str(wlstep) + "\n")
            f.write(str(vmicro) + "\n")
            f.write("tempSynth" + ".syn" + "\n")
            f.write("solar" + "\n")
            f.write("skip! Name of file with stratification data" + "\n")
            f.write("skip vt.dat! Name of file with Vturb vs.depths(km / sec)" + "\n")
            f.write("skip! Name of file with Vrad vs.depths(km / sec)")
            
            
    
    def REWRITE_CONFIG_CONVOLVE(self, vsini, vmacro, R):
        with open(os.path.join(self.pref, self.wd_convolve, self.config_convolve), 'w') as f:
            f.write("tempSynth.bsn" + "\n")
            f.write(str(vsini) + "\t" + "!Vsini (km/sec)" + "\n")
            f.write(str(vmacro) + "\t" + "!Vmacro (km/sec)" + "\n")
            f.write(str(R) + "\t" + "! R==l/dl " + "\n")
            f.write("tempSynth.rgs" + "\n")
            
            
            
    def FIND_QUAD_LIMB_DARK_COEFF(self, mus, Imus):
        
        def limb_darkening(mu, u1, u2):
            return 1 - u1 * (1 - mu) - u2 * (1 - mu)**2
        
        I_max = Imus[0] 
        normalized_intensities = Imus / I_max
        
        params, params_covariance = curve_fit(limb_darkening, mus, normalized_intensities, bounds=(0, [1., 1.]))
        u1, u2 = params
        
        return u1, u2
    
    
    def COMPUTE(self, Teff, logg, met = 0.0, vmicro = 2.0, vsini = 0.0, vmacro = 0.0, R = 50000, wl1 = 4000, wl2 = 5000, wlstep = 0.05):
        
        atmModel = self.FIND_NEAREST_MODEL(Teff, logg, met)
        
        self.REWRITE_CONFIG_SYNTHV(atmModel, vmicro, wl1, wl2, wlstep)
        
        process = subprocess.Popen(
            os.path.join(self.pref, self.wd_synthv, self.exec_synthv),
            cwd=self.wd_synthv,
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        process.wait()
        #if process.returncode != 0:
        #    raise Exception(f"SynthV process failed with return code {process.returncode}")
        
        
        if self.if_convolve:
            
            src = os.path.join(self.pref, self.wd_synthv, "tempSynth.bsn")
            dst = os.path.join(self.pref, self.wd_convolve, "tempSynth.bsn")
            shutil.copy(src, dst)
            
            self.REWRITE_CONFIG_CONVOLVE(vsini, vmacro, R)
            
            process = subprocess.Popen(
                os.path.join(self.pref, self.wd_convolve, self.exec_convolve),
                cwd=self.wd_convolve,
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
            process.wait()
            
         
        synth = np.loadtxt(os.path.join(self.pref, self.wd_synthv, "tempSynth.syn"))
        if self.if_convolve:
            synth = np.loadtxt(os.path.join(self.pref, self.wd_convolve, "tempSynth.rgs"))
            
        #imu = None
        #mu = None
        #if self.if_imu:
        mu = [0.9636,  0.8864,  0.8018,  0.7071,  0.5976,  0.4629,  0.2673]
        imu = np.loadtxt(os.path.join(self.pref, self.wd_synthv, "tempSynth..imu"), skiprows = 1)
        
        
        os.remove(os.path.join(self.pref, self.wd_synthv, atmModel))
        for filename in os.listdir(os.path.join(self.pref, self.wd_synthv)):
            if "tempSynth" in filename:
                os.remove(os.path.join(self.pref, self.wd_synthv, filename))
        if self.if_convolve:
            os.remove(os.path.join(self.pref, self.wd_convolve, "tempSynth.bsn"))
            os.remove(os.path.join(self.pref, self.wd_convolve, "tempSynth.rgs"))
            
            
        return synth, mu, imu, self.grid_logg
    
    def LDC_QUAD_WL(self, Teff, logg, met = 0.0, wl1 = 5000, wl2 = 9000, wlstep = 100, vmicro=2.0):
        
        atmModel = self.FIND_NEAREST_MODEL(Teff, logg, met)
        
        self.REWRITE_CONFIG_SYNTHV(atmModel, vmicro, wl1, wl2, wlstep)
        
        process = subprocess.Popen(
            os.path.join(self.pref, self.wd_synthv, self.exec_synthv),
            cwd=self.wd_synthv,
            shell=True
        )
        process.wait()
        
        
        mu = [0.9636,  0.8864,  0.8018,  0.7071,  0.5976,  0.4629,  0.2673]
        imu = np.loadtxt(os.path.join(self.pref, self.wd_synthv, "tempSynth..imu"), skiprows = 1)
        
        ldcs = []
        for i in range(len(imu[:,0])):
            u1, u2 = self.FIND_QUAD_LIMB_DARK_COEFF(mu, imu[i][8:])
            ldcs.append([imu[i][0], u1, u2])
        ldcs = np.array(ldcs)
        
        os.remove(os.path.join(self.pref, self.wd_synthv, atmModel))
        for filename in os.listdir(os.path.join(self.pref, self.wd_synthv)):
            if "tempSynth" in filename:
                os.remove(os.path.join(self.pref, self.wd_synthv, filename))
        
        return ldcs
        

'''        
if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
            
    sss = synthetic(db_atmos_models = "/Users/nadezhda/Projects/gasp/gasp/bin/Debug/Kurucz/", if_convolve=True, if_imu=False)
    
    s_prim, m_prim, im_prim = sss.COMPUTE(7000, 4.0, vsini=200.0, vmacro = 50.0, wl1=4000, wl2=5000, wlstep = 0.05)
    #s_sec, m_sec, im_sec = sss.COMPUTE(5000, 4.3, vsini=10.0, vmacro = 5.0, wl1=4000, wl2=10000, wlstep = 1)
    
    
    ldcs = sss.LDC_QUAD_WL(7000, 4.0, wl1 = 4000, wl2 = 10000, wlstep = 100)
    plt.plot(ldcs[:,0],ldcs[:,1],"ok")
    plt.plot(ldcs[:,0],ldcs[:,2],"or")
    plt.show()
    
    
    plt.plot(s_prim[:,0], s_prim[:,1], "-k")
    #plt.plot(s_sec[:,0], s_sec[:,1], "-r")
    plt.show()
    
    plt.plot(s_prim[:,0], s_prim[:,2], "-k")
    plt.plot(s_prim[:,0], s_prim[:,3], "-k")
    #plt.plot(s_sec[:,0], s_sec[:,2], "-r")
    #plt.plot(s_sec[:,0], s_sec[:,3], "-r")
    plt.show()
    
    #plt.plot(s_prim[:,0], (s_prim[:,2]+s_sec[:,2])/(s_prim[:,3]+s_sec[:,3]), "-k")
    #plt.show()
'''
    
    