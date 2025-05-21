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
import pandas as pd

from scipy.interpolate import griddata

import scipy.interpolate as sci
import six.moves as smo
import matplotlib.pyplot as plt

import sqlite3

import os

class specSeries:
    def __init__(self, path = None, nspec = None, bjd = None, wls = None, fluxes = None, phases = None):

        self.path, self.nspec, self.bjd, self.wls, self.fluxes_wl, self.phases = path, nspec, bjd, wls, fluxes, phases
        self.fluxes_wl_vars = None

        self.wl_eqlog, self.fluxes_eqlog = None, None
        self.rvs_kms, self.rvs_bin = None, None
        self.lfs, self.sig = None, None

        self.k, self.m, self.n = 2, self.nspec, None

        self.dftobs, self.dftmod, self.dftres = None, None, None
        self.mod, self.resid = None, None

    def LOAD_HERMES(self, path, wl1, wl2):

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
        self.fluxes_wl_vars = []
        for i in range(len(spectra)):
            locfl, locwl = [], []
            for j in range(len(spectra[i][:,0])):
                if wl1 <= spectra[i][j,0] <= wl2:
                    locwl.append(spectra[i][j,0])
                    locfl.append(spectra[i][j,1])

            wls.append(np.array(locwl))
            fluxes.append(np.array(locfl))

            s = pd.Series(locfl)
            rolling_variance = s.rolling(30).var()
            self.fluxes_wl_vars.append(np.mean(rolling_variance))


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


    def LOAD_LIST(self, path, wl1, wl2):

        ll = fnmatch.filter(os.listdir(path), '*.list*')[0]

        fnames = []
        bjd = []

        with open(path+ll, 'r') as file:
            lines = file.readlines()

        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                fnames.append(parts[0])
                try:
                    bjd_value = float(parts[1])
                    bjd.append(bjd_value)
                except ValueError:
                    print(f"Could not convert {parts[1]} to float. Skipping this entry.")


        files = fnmatch.filter(os.listdir(path), '*_cf.fits')


        self.bjd = np.array( bjd )

        files = fnames

        spectra =  [np.loadtxt(path+f) for f in files]

        wls = []
        fluxes = []
        self.fluxes_wl_vars = []
        for i in range(len(spectra)):
            locfl, locwl = [], []
            for j in range(len(spectra[i][:,0])):
                if wl1 <= spectra[i][j,0] <= wl2:
                    locwl.append(spectra[i][j,0])
                    locfl.append(spectra[i][j,1])

            wls.append(np.array(locwl))
            fluxes.append(np.array(locfl))

            s = pd.Series(locfl)
            rolling_variance = s.rolling(30).var()
            self.fluxes_wl_vars.append(np.mean(rolling_variance))

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



    def LOAD_TESTSUIT(self, path, period, t0, wl1, wl2):
        files = fnmatch.filter(os.listdir(path), 'comb_sn*')

        bjd_ar = []
        for nrm_file in files:
            # Extract phase part from filename, assuming it starts at index 14 and goes until the '.dat'
            phase_value = float(nrm_file[14:-4])  # Phase part extracted from the filename
            bjd = t0 + phase_value * period       # Calculate the BJD based on period and t0
            bjd_ar.append((nrm_file, bjd))        # Store filename and BJD as a tuple

        # Sort by BJD in ascending order
        bjd_ar = sorted(bjd_ar, key=lambda x: x[1])

        # Extract sorted filenames and BJDs
        sorted_files = [str(item[0]) for item in bjd_ar]
        self.bjd = np.array([float(item[1]) for item in bjd_ar])

        print("Sorted BJDs:", self.bjd)
        print("Sorted Files:", sorted_files)

        # Load sorted spectra files
        spectra = [np.loadtxt(os.path.join(path, f)) for f in sorted_files]



        self.CALC_PHASES(period,t0)
        print(self.phases)

        wls = []
        fluxes = []
        self.fluxes_wl_vars = []
        for i in range(len(spectra)):
            locfl, locwl = [], []
            for j in range(len(spectra[i][:,0])):
                if wl1 <= spectra[i][j,0] <= wl2:
                    locwl.append(spectra[i][j,0])
                    locfl.append(spectra[i][j,1])

            wls.append(np.array(locwl))
            fluxes.append(np.array(locfl))

            s = pd.Series(locfl)
            rolling_variance = s.rolling(30).var()
            self.fluxes_wl_vars.append(np.mean(rolling_variance))


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
            f = interp1d(log_wavelength, spec, bounds_error=False, fill_value=1.0)
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
                 atmos_mode = "lin_interp",
                 if_pyastr = False,
                 if_convolve = False, wd_convolve = "convolve/", config_convolve = "convolve.conf", exec_convolve = "convolve",
                 db_atmos_models = "Kurucz/",
                 if_imu = False, if_lines=True,
                 abund_tables = "abundances/",
                 pref = "/Users/nadezhda/Dropbox/Binaries/MergeMeth/multi-objective/"):

        self.wd_synthv = wd_synthv
        self.config_synthv = config_synthv
        self.exec_synthv = exec_synthv

        self.atmos_mode = atmos_mode

        self.if_pyastr = if_pyastr
        self.if_convolve = if_convolve
        self.wd_convolve = wd_convolve
        self.config_convolve = config_convolve
        self.exec_convolve = exec_convolve

        self.db_atmos_models = db_atmos_models

        self.if_imu = if_imu
        self.if_lines = if_lines

        self.pref = pref

        self.grid_logg = None

        self.abund_tables = abund_tables

        if self.atmos_mode == "lin_interp":
            self.atmosphere = atmosphere(db_atmos_models)

        self.initialize_database()

    def initialize_database(self):
        #self.db_connection = sqlite3.connect('/STER/spectraspitter/spectra_get_sp.db', check_same_thread=False)
        self.db_connection = sqlite3.connect(self.pref+'spectra_cache.db', check_same_thread=False)
        with self.db_connection:
            self.db_connection.execute('''
                CREATE TABLE IF NOT EXISTS spectra (
                    Teff REAL,
                    logg REAL,
                    metallicity REAL,
                    vsini REAL,
                    spectrum BLOB,
                    imu BLOB,
                    PRIMARY KEY (Teff, logg, metallicity, vsini)
                )
            ''')

    def save_spectrum_to_db(self, Teff, logg, metallicity, vsini, spectrum, imu):
        with self.db_connection:
            self.db_connection.execute('''
                INSERT OR REPLACE INTO spectra (Teff, logg, metallicity, vsini, spectrum, imu)
                VALUES (?, ?, ?, ?, ?, ?)
            ''', (Teff, logg, metallicity, vsini, spectrum.tobytes(), imu.tobytes()))

    def get_spectrum_from_db(self, Teff, logg, metallicity, vsini):

        cursor = self.db_connection.cursor()
        query = '''
            SELECT spectrum, imu, Teff, logg, metallicity, vsini FROM spectra
            WHERE Teff BETWEEN ? AND ?
            AND logg BETWEEN ? AND ?
            AND metallicity BETWEEN ? AND ?
            AND vsini BETWEEN ? AND ?
        '''
        params = (
            Teff - 20.0, Teff + 20.0,
            logg - 0.01, logg + 0.01,
            metallicity - 0.005, metallicity + 0.005,
            vsini - 1.0, vsini + 1.0
        )
        cursor.execute(query, params)
        rows = cursor.fetchall()

        for row in rows:
            bin_synth, bin_imu, db_Teff, db_logg, db_metallicity, db_vsini = row
            if (abs(db_Teff - Teff) <= 20.0 and
                abs(db_logg - logg) <= 0.01 and
                abs(db_metallicity - metallicity) <= 0.005 and
                abs(db_vsini - vsini) <= 1.0):

                synth = np.frombuffer(bin_synth).copy()
                imu = np.frombuffer(bin_imu).copy()
                return synth,imu

        return None, None


    def FIND_NEAREST_MODEL(self, Teff, logg, met):

        atm_path = self.db_atmos_models
        dirs_metdir = [d for d in os.listdir(atm_path) if os.path.isdir(os.path.join(atm_path, d))]

        mpref = "lp"
        if met < -1e-3:
            mpref = "lm"

        met_dir = ""
        diff_met = 1e5

        for d in dirs_metdir:
            if mpref in d:
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



    def REWRITE_CONFIG_SYNTHV(self, atmModel, fabund, vmicro, wl1, wl2, wlstep):

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
            f.write(self.abund_tables + fabund + "\n")
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


    def FIND_LIN_LIMB_DARK_COEFF(self, mus, Imus):

        def limb_darkening(mu, u1):
            return 1 - u1 * (1 - mu)

        I_max = Imus[0]
        normalized_intensities = Imus / I_max

        params, params_covariance = curve_fit(limb_darkening, mus, normalized_intensities, bounds=(0, [1.]))
        u1 = params

        return u1

    def COMPUTE(self, Teff, logg, met, vmicro = 2.0, vsini = 0.0, vmacro = 0.0, R =  85000, wl1 = 4000, wl2 = 5000, wlstep = 0.05):

         #uncomment for cached spectra
        try:
            cached_synth, cached_imu = self.get_spectrum_from_db(Teff, logg, met, vsini)
            if cached_synth is not None:
                # Decompose the cached_spectrum back into the expected format
                # Assuming cached_spectrum is an array where the first column is wavelength and the second is intensity
                synth = cached_synth.reshape(-1, 6)
                # ... the rest of your code for processing the cached spectrum ...
                mu = [0.9636,  0.8864,  0.8018,  0.7071,  0.5976,  0.4629,  0.2673]
                imu = cached_imu.reshape(-1, 15)
                #print(imu)
                #print(synth)
                #print("spectrum recycled!!! ", Teff)
                return synth, mu, imu
        except Exception as e:
            print(e)



        atmModel = None
        #print(self.atmos_mode)

        if self.atmos_mode == "nearest":
            atmModel = self.FIND_NEAREST_MODEL(Teff, logg, met)

        elif self.atmos_mode == "lin_interp":
            try:
                atmModel = f'teff{np.round(Teff,1)}_logg{np.round(logg,2)}_met{np.round(met,2)}.mod'
                #print(atmModel)
                #print(Teff, logg, met)

                interp_atm = self.atmosphere.interpolate_model(Teff, logg, met, os.path.join(self.pref, self.wd_synthv, atmModel))
            except Exception as e:
                print(e)
                #atmModel = self.FIND_NEAREST_MODEL(Teff, logg, met)
                #print("!!! interpolation failed. using nearest atmosphere !!!")


        amet = str(int(np.abs(np.round(met,2))*100)).zfill(5)
        asign = "lp"
        if met < -1e-3:
            asign = "lm"
        fabund = f'{asign}{amet}.abn'

        self.REWRITE_CONFIG_SYNTHV(atmModel, fabund, vmicro, wl1, wl2, wlstep)


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


        #print(os.path.join(self.pref, self.wd_synthv, "tempSynth..imu"))

        #imu = None
        #mu = None
        #if self.if_imu:
        mu = [0.9636,  0.8864,  0.8018,  0.7071,  0.5976,  0.4629,  0.2673]
        imu = np.loadtxt(os.path.join(self.pref, self.wd_synthv, "tempSynth..imu"), skiprows = 1)



        if self.if_convolve:
            synth = np.loadtxt(os.path.join(self.pref, self.wd_convolve, "tempSynth.rgs"))
        elif self.if_pyastr:
            synth = np.loadtxt(os.path.join(self.pref, self.wd_synthv, "tempSynth.syn"))

            if 0 < vsini <= 400:

                ldcs = []
                #print(imu)
                for i in range(0,len(imu[:,0]),100):
                    #print(imu[i][8:])
                    u1 = self.FIND_LIN_LIMB_DARK_COEFF(mu, imu[i][8:])
                    ldcs.append(u1)
                ldcs = np.array(ldcs)
                #print(ldcs)
                linLDC = sum(ldcs)/len(ldcs)
                #print(linLDC)
                #print(vsini)

                synth[:,1] = self.fastRotBroad(synth[:,0], synth[:,1], linLDC, vsini)

            if 0 < R <= 200000:
                #print(R)

                synth[:,1] = self.instrBroadGaussFast(synth[:,0], synth[:,1], R,edgeHandling="firstlast")


        else:
            synth = np.loadtxt(os.path.join(self.pref, self.wd_synthv, "tempSynth.syn"))




        try:
            os.remove(os.path.join(self.pref, self.wd_synthv, atmModel))
        except OSError as e:
            print(f"Error: {e.strerror}")
        for filename in os.listdir(os.path.join(self.pref, self.wd_synthv)):
            if "tempSynth" in filename:
                os.remove(os.path.join(self.pref, self.wd_synthv, filename))
        if self.if_convolve:
            os.remove(os.path.join(self.pref, self.wd_convolve, "tempSynth.bsn"))
            os.remove(os.path.join(self.pref, self.wd_convolve, "tempSynth.rgs"))


        self.save_spectrum_to_db(Teff, logg, met, vsini, synth, imu)

        #print(imu)
        #print(synth)
        #print("Sp saved to db")

        return synth, mu, imu

    def LDC_QUAD_WL(self, Teff, logg, met, wl1 = 4000, wl2 = 9000, wlstep = 100, vmicro=2.0):

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



    def broadGaussFast(self, x, y, sigma, edgeHandling=None, maxsig=None):

        dxs = x[1:] - x[0:-1]

        if maxsig is None:
            lx = len(x)
        else:
            lx = int(((sigma * maxsig) / dxs[0]) * 2.0) + 1

        nx = (np.arange(lx, dtype=int) - sum(divmod(lx, 2)) + 1) * dxs[0]


        def evaluate_h_sig(x, h, sig):

            y = (
                h * np.exp(-((0.0 - x) ** 2) / (2.0 * sig**2))
            )
            return y


        e = evaluate_h_sig( nx, 1.0 / np.sqrt(2 * np.pi * sigma ** 2), sigma)
        e /= np.sum(e)
        if edgeHandling == "firstlast":
            nf = len(y)
            y = np.concatenate((np.ones(nf) * y[0], y, np.ones(nf) * y[-1]))
            result = np.convolve(y, e, mode="same")[nf:-nf]
        else:
            result = np.convolve(y, e, mode="same")

        return result

    def instrBroadGaussFast(self, wvl, flux, resolution, edgeHandling=None, fullout=False, maxsig=None, equid=False):

        wvlorig = wvl
        wvl = np.linspace(wvl[0], wvl[-1], 2*len(wvl))
        flux = sci.interp1d(wvlorig, flux)(wvl)

        meanWvl = np.mean(wvl)
        fwhm = 1.0 / float(resolution) * meanWvl
        sigma = fwhm / (2.0 * np.sqrt(2. * np.log(2.)))

        result = self.broadGaussFast(
            wvl, flux, sigma, edgeHandling=edgeHandling, maxsig=maxsig)

        result = sci.interp1d(wvl, result)(wvlorig)

        if not fullout:
            return result
        else:
            return (result, fwhm)

    def fastRotBroad(self, wvl, flux, epsilon, vsini, effWvl=None):

        dwl = wvl[1] - wvl[0]

        if effWvl is None:
            effWvl = np.mean(wvl)

        gdl = _Gdl(vsini, epsilon)

        # The number of bins needed to create the broadening kernel
        binnHalf = int(np.floor(((vsini / 299792.458) * effWvl / dwl))) + 1
        gwvl = (np.arange(4*binnHalf) - 2*binnHalf) * dwl + effWvl
        # Create the broadening kernel
        dl = gwvl - effWvl
        g = gdl.gdl(dl, effWvl, dwl)
        # Remove the zero entries
        indi = np.where(g > 0.0)[0]
        g = g[indi]

        result = np.convolve(flux, g, mode="same") * dwl
        return result



    def rotBroad(self, wvl, flux, epsilon, vsini, edgeHandling="firstlast"):

        dwl = wvl[1] - wvl[0]

        validIndices = None

        if edgeHandling == "firstlast":
            # Number of bins additionally needed at the edges
            binnu = int(np.floor(((vsini / 299792.458) * max(wvl)) / dwl)) + 1
            # Defined 'valid' indices to be returned
            validIndices = np.arange(len(flux)) + binnu
            # Adapt flux array
            front = np.ones(binnu) * flux[0]
            end = np.ones(binnu) * flux[-1]
            flux = np.concatenate((front, flux, end))
            # Adapt wavelength array
            front = (wvl[0] - (np.arange(binnu) + 1) * dwl)[::-1]
            end = wvl[-1] + (np.arange(binnu) + 1) * dwl
            wvl = np.concatenate((front, wvl, end))
        else:
            validIndices = np.arange(len(flux))

        result = np.zeros(len(flux))
        gdl = _Gdl(vsini, epsilon)

        for i in smo.range(len(flux)):
            dl = wvl[i] - wvl
            g = gdl.gdl(dl, wvl[i], dwl)
            result[i] = np.sum(flux * g)
        result *= dwl

        return result[validIndices]


class _Gdl:

    def __init__(self, vsini, epsilon):

        self.vc = vsini / 299792.458
        self.eps = epsilon

    def gdl(self, dl, refwvl, dwl):

        self.dlmax = self.vc * refwvl
        self.c1 = 2.*(1. - self.eps) / \
            (np.pi * self.dlmax * (1. - self.eps/3.))
        self.c2 = self.eps / (2. * self.dlmax * (1. - self.eps/3.))
        result = np.zeros(len(dl))
        x = dl/self.dlmax
        indi = np.where(np.abs(x) < 1.0)[0]
        result[indi] = self.c1 * \
            np.sqrt(1. - x[indi]**2) + self.c2*(1. - x[indi]**2)

        result /= (np.sum(result) * dwl)
        return result


class atmosphere:

    def __init__(self, path):
        self.path = path

        self.teff_grid, self.logg_grid, self.metallicity_grid = self.collect_model_parameters()


    def parse_model_filename(self, filename):
        parts = filename.split('_')
        metallicity_sign = -1 if parts[0][1] == 'm' else 1
        metallicity = metallicity_sign * float(parts[0][3:]) / 10
        teff = int(parts[1])
        logg = float(parts[2]) / 100
        return teff, logg, metallicity


    def find_linear_neighbors(self, teff, logg, metallicity):
        """
        Find the 8 corner neighbors forming a cube around the desired (teff, logg, metallicity) point.
        """
        def find_bounds(value, grid):
            grid = np.array(grid)
            idx = np.searchsorted(grid, value)
            idx = np.clip(idx, 1, len(grid) - 1)
            lower_idx = idx - 1
            upper_idx = idx
            return lower_idx, upper_idx

        # Get the bounding indices in each dimension
        teff_lower_idx, teff_upper_idx = find_bounds(teff, self.teff_grid)
        logg_lower_idx, logg_upper_idx = find_bounds(logg, self.logg_grid)
        metallicity_lower_idx, metallicity_upper_idx = find_bounds(metallicity, self.metallicity_grid)

        # Create a list of the 8 neighbors
        neighbors = []
        for t_idx in [teff_lower_idx, teff_upper_idx]:
            for g_idx in [logg_lower_idx, logg_upper_idx]:
                for m_idx in [metallicity_lower_idx, metallicity_upper_idx]:
                    neighbor = (self.teff_grid[t_idx], self.logg_grid[g_idx], self.metallicity_grid[m_idx])
                    neighbors.append(neighbor)

        return neighbors


    def find_quadratic_neighbors(self, teff, logg, metallicity):
        """
        Find the neighbors in a 3x3x3 grid around the desired (teff, logg, metallicity) point.
        """
        # Helper function to find the closest indices in a grid
        def find_bounds(value, grid):
            grid = np.array(grid)
            idx = np.searchsorted(grid, value)
            idx = np.clip(idx, 1, len(grid) - 2)
            return idx - 1, idx, idx + 1

        # Get the bounding indices in each dimension
        teff_idxs = find_bounds(teff, self.teff_grid)
        logg_idxs = find_bounds(logg, self.logg_grid)
        metallicity_idxs = find_bounds(metallicity, self.metallicity_grid)

        # Create a list of the 27 neighbors
        neighbors = []
        for t_idx in teff_idxs:
            for g_idx in logg_idxs:
                for m_idx in metallicity_idxs:
                    neighbors.append((self.teff_grid[t_idx], self.logg_grid[g_idx], self.metallicity_grid[m_idx]))

        return neighbors

    def collect_model_parameters(self):

        '''
        base_directory = self.path

        teff_values = []
        logg_values = []
        metallicity_values = []


        for subdir, _, files in os.walk(base_directory):
            for file in files:
                if file.endswith('.mod'):
                    teff, logg, metallicity = self.parse_model_filename(file)
                    teff_values.append(teff)
                    logg_values.append(logg)
                    metallicity_values.append(metallicity)

        teff_values = np.unique(teff_values)
        logg_values = np.unique(logg_values)
        metallicity_values = np.unique(metallicity_values)
        '''

        teff_values = np.concatenate((np.arange(2500,9901,100), np.arange(10000,30001,250), np.arange(30000,34001,500)), axis=None)
        logg_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                       1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                       2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
                       3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0,
                       4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0]
        metallicity_values = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

        #print(teff_values, logg_values, metallicity_values)
        #input()

        return teff_values, logg_values, metallicity_values

    def get_neighbours(self, desired_teff, desired_logg, desired_metallicity):

        #base_directory = self.path

        neighbors = self.find_linear_neighbors(desired_teff, desired_logg, desired_metallicity)

        neighbor_files = []
        #for subdir, _, files in os.walk(base_directory):
        #    for neighbor in neighbors:
        #        neighbor_teff, neighbor_logg, neighbor_metallicity = neighbor
        #        for file in files:
        #            if file.endswith('.mod'):
        #                teff, logg, metallicity = self.parse_model_filename(file)
        #                if (teff == neighbor_teff and logg == neighbor_logg and
        #                        metallicity == neighbor_metallicity):
        #                    neighbor_files.append(os.path.join(subdir, file))
        for neighbor in neighbors:
            neighbor_teff, neighbor_logg, neighbor_metallicity = neighbor
            ateff = str(int(np.round(neighbor_teff,1))).zfill(5)
            alogg = str(int(np.round(neighbor_logg*100.0,0))).zfill(4)
            amet = str(int(np.abs(np.round(neighbor_metallicity,1))*10)).zfill(4)
            amet_dir = str(int(np.abs(neighbor_metallicity)*10)).zfill(2)
            asign = "lp"
            if neighbor_metallicity < -1e-3:
                asign = "lm"
            conv = "off"
            if neighbor_teff < 9000:
                conv = "on"
            fn = self.path+f'{asign}{amet_dir}k2/{asign}{amet}_{ateff}_{alogg}_0020_{conv}.mod'
            #print(neighbor)
            #print(fn)
            neighbor_files.append(fn)


        return neighbor_files

    def read_atmosphere_file(self,filename):

        start_reading = False


        # RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB

        Rhox = []
        T = []
        P = []
        Xne = []
        Abross = []
        Accrad = []
        Vturb = []

        num_layers = 0

        with open(filename, 'r') as file:
            for line in file:
                # Check for the start of the numerical data
                if line.startswith('READ DECK6'):
                    start_reading = True
                    try:
                        num_layers = int(line[11:13])
                    except Exception as e:
                        print(e)
                    continue

                # Read the numerical data
                if start_reading:
                    values = line.split()
                    #print(values)
                    if len(values) >= 6:
                        Rhox.append(float(values[0]))
                        T.append(float(values[1]))
                        P.append(float(values[2]))
                        Xne.append(float(values[3]))
                        Abross.append(float(values[4]))
                        Accrad.append(float(values[5]))
                        Vturb.append(float(values[6]))
                    else:
                        break


        #print(num_layers)

        Rhox = np.array(Rhox)
        T = np.array(T)
        P = np.array(P)
        Xne = np.array(Xne)
        Abross = np.array(Abross)
        Accrad = np.array(Accrad)
        Vturb = np.array(Vturb)

        return Rhox, T, P, Xne, Abross, Accrad, Vturb




    def interpolate_model(self, desired_teff, desired_logg, desired_metallicity, output_filename):

        neighbor_files = self.get_neighbours(desired_teff, desired_logg, desired_metallicity)
        #print(neighbor_files)
        if len(neighbor_files) != 8:
            return

        min_layers = min(len(self.read_atmosphere_file(file)[0]) for file in neighbor_files)
        num_layers = min_layers
        #print(num_layers)
        #input()

        points = []  # Teff, logg, metallicity
        values = []  # Parameter values for interpolation
        interpolated_layers = np.zeros((num_layers, 7))

        for file in neighbor_files:
            teff, logg, metallicity = self.parse_model_filename(os.path.basename(file))
            atmos_params = self.read_atmosphere_file(file)
            for layer_params in zip(*atmos_params):
                points.append((teff, logg, metallicity))
                values.append(layer_params)
                #print((teff, logg, metallicity))
                #print(layer_params)
            #input()

        for layer_idx in range(num_layers):
            # Collect the data for this layer from each neighbor file
            points = []  # Teff, logg, metallicity
            values = []  # Parameter values for this layer

            for file in neighbor_files:
                teff, logg, metallicity = self.parse_model_filename(os.path.basename(file))
                atmos_params = self.read_atmosphere_file(file)
                # Append the parameters for the current layer
                points.append((teff, logg, metallicity))
                # Ensure that we only append the parameters if we're not out of bounds for this file
                if layer_idx < len(atmos_params[0]):
                    values.append([param[layer_idx] for param in atmos_params])

            # Convert lists to NumPy arrays
            points = np.array(points)
            values = np.array(values)

            # Perform the interpolation for this layer
            for param_index in range(len(values[0])):  # Loop over each parameter
                param_values = values[:, param_index]
                interpolated_value = griddata(points, param_values, (desired_teff, desired_logg, desired_metallicity), method='linear')
                interpolated_layers[layer_idx, param_index] = interpolated_value if interpolated_value.size > 0 else np.nan


        with open(output_filename, 'w') as file:
            # Write the header
            file.write(f"TEFF   {desired_teff:.1f}  GRAVITY {desired_logg:.5f} LTE\n")
            file.write("TITLE ---\n")

            file.write(f"READ DECK6 {len(interpolated_layers)} RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB\n")

            for layer_values in interpolated_layers:
                file.write(f" {layer_values[0]:.8E} {layer_values[1]:8.1f} {layer_values[2]:.3E} {layer_values[3]:.3E} "
                           f"{layer_values[4]:.5E} {layer_values[5]:.5E} {layer_values[6]:.3E}\n")

            file.write("BEGIN\n")




        return interpolated_layers
