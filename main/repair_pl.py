import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import chebyshev
import math
import random

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit, minimize
from scipy.integrate import simps

from orbit import orbit
from star import star
from binary import binary

from lcrvcurve import lightcurve

from spectra import specSeries
from spectra import synthetic

from deap import base, creator, tools, algorithms

from multiprocessing import Pool, current_process, TimeoutError

import pickle, os, json, shutil

import seaborn as sns

import pandas as pd

from scipy.stats.qmc import Sobol
from scipy.spatial.distance import euclidean

global pars_opt, param_bounds, mainconfig

import logging

logging.raiseExceptions = False


sns.set_style("white")
sns.set_style("ticks")
sns.set_palette("rocket")


with open('config.json', 'r') as config_file:
    mainconfig = json.load(config_file)


dlevel = logging.WARNING
if mainconfig["loglevel"] == "debug":
    dlevel = logging.DEBUG
elif mainconfig["loglevel"] == "info":
    dlevel = logging.INFO
elif mainconfig["loglevel"] == "warning":
    dlevel = logging.WARNING
elif mainconfig["loglevel"] == "error":
    dlevel = logging.ERROR
elif mainconfig["loglevel"] == "critical":
    dlevel = logging.CRITICAL

logging.basicConfig(filename=f'{mainconfig["saveto"]}/main.log', level=dlevel, format='%(asctime)s - %(levelname)s - %(message)s')

synthini_s = []
synthini_simu = []
for idir in range(mainconfig["nthreads"]):
    #if os.path.exists(f'SynthV_core{idir}'):
   #     shutil.rmtree(f'SynthV_core{idir}', ignore_errors=True)
   # os.mkdir(f'SynthV_core{idir}')
   # shutil.copytree(mainconfig["synthVpath"], f'SynthV_core{idir}')

    #if os.path.exists(f'convolve_core{idir}'):
    #    shutil.rmtree(f'convolve_core{idir}', ignore_errors=True)
   # os.mkdir(f'convolve_core{idir}')
   # shutil.copytree(mainconfig["convolvepath"], f'convolve_core{idir}')


    synthini_s.append(synthetic(wd_synthv=f'SynthV_core{idir}', wd_convolve=f'convolve_core{idir}', db_atmos_models = mainconfig["atmmodels"], if_pyastr=True, if_convolve=False, if_imu=False, atmos_mode="lin_interp", abund_tables=mainconfig["abunds"], pref = os.getcwd()))
    synthini_simu.append(synthetic(wd_synthv=f'SynthV_core{idir}', db_atmos_models = mainconfig["atmmodels"], if_convolve=False, if_imu=True, if_lines = False, atmos_mode="lin_interp", abund_tables=mainconfig["abunds"],  pref = os.getcwd()))


sobol_engine = Sobol(d=len(mainconfig["params_opt"]), scramble=True)

specs = specSeries()
if mainconfig["obsmode"] == "test":
    specs.LOAD_TESTSUIT(mainconfig["path_specs"], period = mainconfig["params_init"]["porb"], t0 = mainconfig["params_init"]["t0"], wl1 = mainconfig["wl11"], wl2 = mainconfig["wl21"])
elif mainconfig["obsmode"] == "obs-list":
    specs.LOAD_LIST(mainconfig["path_specs"], wl1 = mainconfig["wl11"], wl2 = mainconfig["wl21"])
else:
    specs.LOAD_HERMES(mainconfig["path_specs"], wl1 = mainconfig["wl11"], wl2 = mainconfig["wl21"])

specs.TRIM(mainconfig["wl11"]+10, mainconfig["wl21"]-10)



lcf = np.loadtxt(mainconfig["path_lc"])
if mainconfig["obsmode"] == "test":
    lcs = lightcurve(lcf[:,0], lcf[:,1], np.zeros_like(lcf[:,0]))
else:
    lcs = lightcurve(lcf[:,0], lcf[:,1], lcf[:,2])

passband = np.loadtxt("tess_transmission.txt", skiprows=6,delimiter=",")
passband[:,0] *= 10.0

if "custom_passband" in mainconfig:
    passband = np.loadtxt(mainconfig["custom_passband"], delimiter="\t")

mu_new = np.arange(0, 1.001, 0.025)

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException

def integrate_passband(wavelengths, flux, ifTess = True):

    if ifTess:
        filter_wavelengths = passband[:,0]
        filter_response =  passband[:,1]
    else:
        filter_wavelengths = np.linspace(3000.0, 11000.0, 1000)
        filter_response = np.ones_like(filter_wavelengths)


    min_wavelength = max(min(filter_wavelengths), min(wavelengths))
    max_wavelength = min(max(filter_wavelengths), max(wavelengths))


    if min_wavelength > max_wavelength:
        return 0.0

    common_wavelengths = np.linspace(min_wavelength, max_wavelength, 1000)

    filter_interp = interp1d(filter_wavelengths, filter_response, bounds_error=False, fill_value=0.0)
    spectrum_interp = interp1d(wavelengths, flux, fill_value="extrapolate")

    filter_values = filter_interp(common_wavelengths)
    spectrum_values = spectrum_interp(common_wavelengths)

    numerator = simps(spectrum_values * filter_values, common_wavelengths)
    denominator = simps(filter_values, common_wavelengths)

    if np.abs(numerator-denominator) < 0.0001 : return numerator

    return numerator / denominator if denominator != 0 else 0.0


def FIND_QUAD_LIMB_DARK_COEFF(mus, Imus):

    def limb_darkening(mu, u1, u2):
        return 1 - u1 * (1 - mu) - u2 * (1 - mu)**2

    I_max = Imus[0]
    normalized_intensities = Imus / I_max

    params, params_covariance = curve_fit(limb_darkening, mus, normalized_intensities, bounds=(0, [1., 1.]))
    u1, u2 = params

    return u1, u2

def relative_difference(q1,q2):
    return np.abs(q1-q2)/q1

def interpolate_intensity(mu, Imu):

    f = interp1d(mu, Imu, kind='quadratic', fill_value="extrapolate")

    Imu_new = f(mu_new)
    Imu_new /= Imu_new[-1]

    return Imu_new


def fetch_spec_new(obs_wl, obs_flux, synth_wl, synth_flux):

    synth_flux_resampled = np.interp(obs_wl, synth_wl, synth_flux)

    obs_flux_corrected = 1 - ((1 - obs_flux) - (1 - obs_flux).min()) / (1 - obs_flux).ptp()
    synth_flux_corrected = 1 - ((1 - synth_flux_resampled) - (1 - synth_flux_resampled).min()) / (1 - synth_flux_resampled).ptp()

    residual = obs_flux_corrected - synth_flux_corrected
    def fit_func(wl, *coeffs):
        return chebyshev.chebval(wl, coeffs)

    degree = 5
    initial_coeffs = [0] * (degree + 1)

    popt, _ = curve_fit(fit_func, obs_wl, residual, p0=initial_coeffs, maxfev=5000)

    obs_flux_final = obs_flux_corrected - fit_func(obs_wl, *popt)

    obs_flux_final = 1-(1-obs_flux_final)*(np.mean(1-synth_flux) / np.mean(1-obs_flux_final) )

    return  obs_flux_final



def update_system(mainconfig, system, params, synthV_instance, synthV_imu_instance, ifplot_main = False, if_save_mods = False, if_sbratio_spec = True, titleprefix=f'', obj_type = 'chi2'):

    system.UPDATE_PARAMS(params)

    logging.info(f'Primary Teff {system.primary.Teff}, logg {system.primary.logg}, vsini {system.primary.vsini}')
    logging.info(f'Secondary Teff {system.secondary.Teff}, logg {system.secondary.logg}, vsini {system.secondary.vsini}')


    b_prim, b_sec = 1.0, 1.0

    try:

        system.primary.spectrum, m_prim, im_prim = synthV_instance.COMPUTE(system.primary.Teff, system.primary.logg, system.primary.metallicity, R=mainconfig["spres"], vsini=system.primary.vsini, vmacro = mainconfig["vmacro1"], wl1=mainconfig["wl11"], wl2=mainconfig["wl21"], wlstep = mainconfig["wlstep1"])
        nosp, mu, imu = synthV_imu_instance.COMPUTE(system.primary.Teff, system.primary.logg, system.primary.metallicity, R=0, wl1=passband[0][0]-100, wl2=passband[-1][0]+100,  wlstep = 5)
        nosp[:,2:4] *= 1e5
        system.primary.spectrum_calib = nosp

        if if_save_mods:
            with open(f"primary_{np.round(system.primary.Teff,1)}_{np.round(system.primary.logg,1)}.spectrum", "w")  as out:
                for i in range(len(system.primary.spectrum[:,0])):
                    out.write(str(system.primary.spectrum[i][0]) + "\t" + str(system.primary.spectrum[i][1])+ "\t" + str(system.primary.spectrum[i][2])+ "\t" + str(system.primary.spectrum[i][3]) + "\n" )


        b_prim = integrate_passband(system.primary.spectrum_calib[:,0], system.primary.spectrum_calib[:,2])

        if if_sbratio_spec:
            b_prim_spec = integrate_passband(system.primary.spectrum[:,0], system.primary.spectrum[:,2], ifTess=False)

        integrated_Imu = []

        for i in range(7):
            integrated_Imu.append(integrate_passband(imu[:,0], imu[:,1+i]))


        mugrid = interpolate_intensity(mu,  integrated_Imu)

        system.primary.ldc = mugrid

        if if_sbratio_spec:
            integrated_Imu = []
            for i in range(7):
                integrated_Imu.append(integrate_passband(im_prim[:,0], im_prim[:,1+i], ifTess=False))

            mugrid = interpolate_intensity(m_prim,  integrated_Imu)

            system.primary.ldc_spec = mugrid


    except Exception as e:
        print("!!exception Synth 1!!")
        print(e)
        if obj_type == "logP":
            system.obj_synth1 = -1e10
            system.obj_synth2 = -1e10
            system.obj_lc = -1e10
            system.obj_dis = -1e10
        else:
            system.obj_synth1 = 1e3
            system.obj_synth2 = 1e3
            system.obj_lc = 1e3
            system.obj_dis = 1e3
        return system

    try:
        system.secondary.spectrum, m_sec, im_sec = synthV_instance.COMPUTE(system.secondary.Teff, system.secondary.logg, system.secondary.metallicity, R=mainconfig["spres"], vsini=system.secondary.vsini, vmacro = mainconfig["vmacro2"], wl1=mainconfig["wl12"], wl2=mainconfig["wl22"], wlstep = mainconfig["wlstep2"])
        nosp, mu, imu = synthV_imu_instance.COMPUTE(system.secondary.Teff, system.secondary.logg, system.secondary.metallicity, R=0, wl1=passband[0][0]-100, wl2=passband[-1][0]+100,  wlstep = 5)
        nosp[:,2:4] *= 1e5
        system.secondary.spectrum_calib = nosp

        if if_save_mods:
            with open(f"secondary_{np.round(system.secondary.Teff,1)}_{np.round(system.secondary.logg,1)}.spectrum", "w")  as out:
                for i in range(len(system.secondary.spectrum[:,0])):
                    out.write(str(system.secondary.spectrum[i][0]) + "\t" + str(system.secondary.spectrum[i][1])+ "\t" + str(system.secondary.spectrum[i][2])+ "\t" + str(system.secondary.spectrum[i][3]) + "\n" )


        b_sec = integrate_passband(system.secondary.spectrum_calib[:,0], system.secondary.spectrum_calib[:,2])

        if if_sbratio_spec:
            b_sec_spec = integrate_passband(system.secondary.spectrum[:,0],  system.secondary.spectrum[:,2], ifTess=False)

        integrated_Imu = []
        for i in range(7):
            integrated_Imu.append(integrate_passband(imu[:,0], imu[:,1+i]))

        mugrid = interpolate_intensity(np.array(mu),  np.array(integrated_Imu))

        system.secondary.ldc = mugrid

        if if_sbratio_spec:
            integrated_Imu = []

            for i in range(7):
                integrated_Imu.append(integrate_passband(im_sec[:,0], im_sec[:,1+i], ifTess=False))

            mugrid = interpolate_intensity(m_sec,  integrated_Imu)

            system.secondary.ldc_spec = mugrid



    except Exception as e:
        print("!!exception Synth 2!!")
        print(e)
        if obj_type == "logP":
            system.obj_synth1 = -1e10
            system.obj_synth2 = -1e10
            system.obj_lc = -1e10
            system.obj_dis = -1e10
        else:
            system.obj_synth1 = 1e3
            system.obj_synth2 = 1e3
            system.obj_lc = 1e3
            system.obj_dis = 1e3
        return system


    system.sbratio = b_sec / b_prim
    logging.info("sbratio: ", system.sbratio)
    system.sbratio_spec = b_sec_spec / b_prim_spec
    logging.info("sbratio_spec: ", system.sbratio_spec)

    try:
        system.UPDATE_MODEL_LC()
    except Exception as e:
        logging.warning("!!exception LC!!")
        logging.warning(e)
        if obj_type == "logP":
            system.obj_synth1 = -1e10
            system.obj_synth2 = -1e10
            system.obj_lc = -1e10
            system.obj_dis = -1e10
        else:
            system.obj_synth1 = 1e3
            system.obj_synth2 = 1e3
            system.obj_lc = 1e3
            system.obj_dis = 1e3
        return system

    if obj_type == "logP":
        residuals_lc = system.residLCNoC
        residuals_lc = residuals_lc[np.isfinite(residuals_lc)]
        logP_lc = -0.5*float(len(residuals_lc))*np.log(2*np.pi*system.obsLC.flux_var) - 1.0/2.0/system.obsLC.flux_var*np.sum(residuals_lc)
        #print(f'LC: logP {logP_lc}, varM {0.5*float(len(residuals_lc))*np.log(2*np.pi*system.obsLC.binned_flux_var)}, sumSqResid {np.sum(residuals_lc)}')
        system.obj_lc = logP_lc
    else:
        residuals_lc = system.residLC
        residuals_lc = residuals_lc[np.isfinite(residuals_lc)]
        schi2lc = np.sum(residuals_lc)

        if schi2lc < 1e-6:
            schi2lc = 1e3
        system.obj_lc = schi2lc

    try:
        system.UPDATE_SP_SEPARATION(ifplot=False, if_lightfact = True)
    except Exception as e:
        logging.error("!!exception DIS!!")
        logging.error(e)
        if obj_type == "logP":
            system.obj_synth1 = -1e10
            system.obj_synth2 = -1e10
            system.obj_lc = -1e10
            system.obj_dis = -1e10
        else:
            system.obj_synth1 = 1e3
            system.obj_synth2 = 1e3
            system.obj_lc = 1e3
            system.obj_dis = 1e3
        return system

    if obj_type == "logP":
        residuals_sp = system.residSpec
        logP_dis = 0
        for i in range(len(residuals_sp[0])):
            aa = np.sum(residuals_sp[:,i]**2 / system.obsSpec.fluxes_wl_vars[i] + np.log(2*np.pi*system.obsSpec.fluxes_wl_vars[i]))
            #print(f'DIS {i}: logP {-0.5*aa}, varM {0.5*float(len(residuals_sp[:,i]))*np.log(2*np.pi*system.obsSpec.fluxes_wl_vars[i])}, sumSqResid {np.sum(residuals_sp[:,i]**2)}')
            logP_dis += aa
        logP_dis *= -0.5

        system.obj_dis = logP_dis
    else:
        residuals_sp = system.residSpec
        residuals_spSqN = system.residSpecSqNorm
        chi2_dis = 0
        for i in range(len(residuals_sp[0])):
            aa = np.sum(residuals_spSqN[:,i])
            chi2_dis += aa

        system.obj_dis = chi2_dis

    if obj_type == "logP":

        ffi = interp1d(system.primary.spectrum[:,0], system.primary.spectrum[:,1], bounds_error=False, fill_value=1.0)
        mm = ffi(system.obsSpec.wl_eqlog_to_wl)
        ob1 = fetch_spec_new(system.obsSpec.wl_eqlog_to_wl, system.obsSpec.mod[0], system.obsSpec.wl_eqlog_to_wl, mm)
        residuals_sp1= (np.array(ob1) - np.array(mm))**2
        logP_syn1 = -0.5*float(len(residuals_sp1))*np.log(2*np.pi*system.mod1_var) - 1.0/2.0/system.mod1_var*np.sum(residuals_sp1)

        ffi = interp1d(system.secondary.spectrum[:,0], system.secondary.spectrum[:,1], bounds_error=False, fill_value=1.0)
        mm = ffi(system.obsSpec.wl_eqlog_to_wl)
        ob2 = fetch_spec_new(system.obsSpec.wl_eqlog_to_wl, system.obsSpec.mod[1], system.obsSpec.wl_eqlog_to_wl, mm)
        residuals_sp2= (np.array(ob2) - np.array(mm))**2
        logP_syn2 = -0.5*float(len(residuals_sp2))*np.log(2*np.pi*system.mod2_var) - 1.0/2.0/system.mod2_var*np.sum(residuals_sp2)


        system.obj_synth1 = logP_syn1
        system.obj_synth2 = logP_syn2

        #print(f'Syn1: logP {logP_syn1}, varM {0.5*float(len(residuals_sp1))*np.log(2*np.pi*system.mod1_var)}, sumSqResid {np.sum(residuals_sp1)}')
        #print(f'Syn2: logP {logP_syn2}, varM {0.5*float(len(residuals_sp2))*np.log(2*np.pi*system.mod2_var)}, sumSqResid {np.sum(residuals_sp2)}')

        print(system.obj_lc, system.obj_dis, system.obj_synth1, system.obj_synth2)

        logging.debug(f'objs: {system.obj_lc} {system.obj_dis} {system.obj_synth1} {system.obj_synth2}')
    else:
        ffi = interp1d(system.primary.spectrum[:,0], system.primary.spectrum[:,1], bounds_error=False, fill_value=1.0)
        mm = ffi(system.obsSpec.wl_eqlog_to_wl)
        ob1 = fetch_spec_new(system.obsSpec.wl_eqlog_to_wl, system.obsSpec.mod[0], system.obsSpec.wl_eqlog_to_wl, mm)
        residuals_sp1= (np.array(ob1) - np.array(mm))**2/ np.array(mm)
        schi2synth1 = np.sum(residuals_sp1)
        ffi = interp1d(system.secondary.spectrum[:,0], system.secondary.spectrum[:,1], bounds_error=False, fill_value=1.0)
        mm = ffi(system.obsSpec.wl_eqlog_to_wl)
        ob2 = fetch_spec_new(system.obsSpec.wl_eqlog_to_wl, system.obsSpec.mod[1], system.obsSpec.wl_eqlog_to_wl, mm)
        residuals_sp2= (np.array(ob2) - np.array(mm))**2/ np.array(mm)
        schi2synth2 = np.sum(residuals_sp2)

        system.obj_synth1 = schi2synth1
        system.obj_synth2 = schi2synth2

        print(system.obj_lc, system.obj_dis, system.obj_synth1, system.obj_synth2)

        logging.debug(f'objs: {system.obj_lc} {system.obj_dis} {system.obj_synth1} {system.obj_synth2}')

    if ifplot_main:

        def plot_spectrum(ax, wl_range, ob, system_spectrum, system_wl_eqlog):
            indices = np.where((system_wl_eqlog >= wl_range[0]) & (system_wl_eqlog <= wl_range[1]))
            ax.plot(system_wl_eqlog[indices], ob[indices], lw=0.5, c='cornflowerblue')

            system_spectrum = np.array(system_spectrum)
            spectrum_indices = np.where((system_spectrum[:,0] >= wl_range[0]) & (system_spectrum[:,0] <= wl_range[1]))
            ax.plot(system_spectrum[spectrum_indices,0].flatten(), system_spectrum[spectrum_indices,1].flatten(), lw=0.5, c="crimson")

            ax.set_xlim(wl_range)

        fig = plt.figure()
        fig.set_figheight(12)
        fig.set_figwidth(9)

        ax_lc = plt.subplot2grid(shape=(5, 2), loc=(0, 0), colspan=2)

        ax_dis = plt.subplot2grid(shape=(5, 2), loc=(1, 0), colspan=2)

        ax1 = plt.subplot2grid(shape=(5, 2), loc=(2, 0), colspan=1)
        ax2 = plt.subplot2grid(shape=(5, 2), loc=(2, 1), colspan=1)
        ax3 = plt.subplot2grid(shape=(5, 2), loc=(3, 0), colspan=1)
        ax4 = plt.subplot2grid(shape=(5, 2), loc=(3, 1), colspan=1)
        ax5 = plt.subplot2grid(shape=(5, 2), loc=(4, 0), colspan=1)
        ax6 = plt.subplot2grid(shape=(5, 2), loc=(4, 1), colspan=1)

        ax_lc.plot(system.obsLC.binned_phase-1, system.obsLC.binned_flux,marker = 'o',markersize=2,c='cornflowerblue')
        ax_lc.plot(system.obsLC.binned_phase-1, system.modLC,lw=1, c="crimson")
        ax_lc.plot(system.obsLC.binned_phase, system.obsLC.binned_flux,marker = 'o',markersize=2,c='cornflowerblue')
        ax_lc.plot(system.obsLC.binned_phase, system.modLC,lw=1, c="crimson")
        ax_lc.plot(system.obsLC.binned_phase+1, system.obsLC.binned_flux,marker = 'o',markersize=2,c='cornflowerblue')
        ax_lc.plot(system.obsLC.binned_phase+1, system.modLC,lw=1, c="crimson")
        ax_lc.set_xlim([-0.4,1.4])

        shift = (np.average(np.abs(residuals_sp[:,0]-residuals_sp[:,-1])))*1.5
        for i in range(len(residuals_sp[0])):
            c="cornflowerblue"
            if i%2 == 0:
                c="crimson"
            ax_dis.plot(residuals_sp[:,i] - shift*i, lw=0.5, c=c, alpha=0.5)

        ax_dis.set_ylim([-(len(residuals_sp[0])+1)*shift, shift])
        ax_dis.set_xlim([100,len(residuals_sp[:,0])-100])

        plot_spectrum(ax1, [mainconfig["wl11"], mainconfig["wl21"]], ob1, system.primary.spectrum, system.obsSpec.wl_eqlog_to_wl)
        plot_spectrum(ax3, [4300,4380], ob1, system.primary.spectrum, system.obsSpec.wl_eqlog_to_wl)
        plot_spectrum(ax5, [4450,4490], ob1, system.primary.spectrum, system.obsSpec.wl_eqlog_to_wl)

        plot_spectrum(ax2, [mainconfig["wl12"], mainconfig["wl22"]], ob2, system.secondary.spectrum, system.obsSpec.wl_eqlog_to_wl)
        plot_spectrum(ax4, [4300,4380], ob2, system.secondary.spectrum, system.obsSpec.wl_eqlog_to_wl)
        plot_spectrum(ax6, [4450,4490], ob2, system.secondary.spectrum, system.obsSpec.wl_eqlog_to_wl)

        plt.suptitle(f'{titleprefix}a_{np.round(np.abs(system.orbit.a),1)}, r1_{np.round(np.abs(system.primary.radius_a),2)}, r2_{np.round(np.abs(system.secondary.radius_a),2)}, M1_{np.round(system.primary.mass,2)} (logg1_{system.primary.logg:.2f}), q_{np.round(system.orbit.q,2)}, incl_{np.round(system.orbit.inc_deg,1)}, Teff1_{np.round(system.primary.Teff,1)}, Teff2_{np.round(system.secondary.Teff,1)}, met1_{np.round(system.primary.metallicity,2)} ')

        plt.savefig(f'{mainconfig["saveto"]}/{titleprefix}Teff1_{system.primary.Teff:.1f}_Teff2_{system.secondary.Teff:.1f}_logg1_{system.primary.logg:.2f}_logg2_{system.secondary.logg:.2f}.png',dpi=150)
        plt.show()

    if if_save_mods:
        with open(f"T{np.round(system.orbit.porb,1)}_M1{np.round(system.primary.mass,1)}_inc{system.orbit.inc_deg}.lc", "w")  as out:
            for i in range(len(system.obsLC.binned_phase)):
                out.write(str(system.obsLC.binned_phase[i]) + "\t" + str(system.modLC[i]) + "\n" )
        with open(f"T{np.round(system.orbit.porb,1)}_M1{np.round(system.primary.mass,1)}_inc{system.orbit.inc_deg}.rv", "w")  as out:
            for i in range(len(system.obsSpec.phases)):
                out.write(str(system.obsSpec.phases[i]) + "\t" + str(system.modRV1[i])+ "\t" + str(system.modRV2[i]) + "\n" )
        with open(f"T{np.round(system.orbit.porb,1)}_M1{np.round(system.primary.mass,1)}_inc{system.orbit.inc_deg}.params", "w")  as out:
            for pp in params:
                out.write(pp + "\t" + str(params[pp]) + "\n")
    return system



def dict_to_array(param_dict):
    return np.array(list(param_dict.values()))

def array_to_dict(param_array, template_dict):
    return {key: val for key, val in zip(template_dict, param_array)}

def objective(mainconfig, params_array, process_id):

    objective_type = mainconfig.get('objective_type', 'chi2') # chi2 = default value if key is not specified in config file

    synthV_instance = synthini_s[process_id - 1]
    synthV_imu_instance = synthini_simu[process_id - 1]

    transformed_params = {
        "sbratio" : 1.0,
        }

    orb = orbit(porb = mainconfig["params_init"]["porb"], inc_deg = mainconfig["params_init"]["incl"], t0 = mainconfig["params_init"]["t0"], a = mainconfig["params_init"]["a"], q = mainconfig["params_init"]["q"], f_c = mainconfig["params_init"]["f_c"], f_s = mainconfig["params_init"]["f_s"])

    system_local = binary(orbit = orb, primary = star(), secondary = star(), sbratio = transformed_params["sbratio"], l3 = mainconfig["params_init"]["l3"], gamma = mainconfig["params_init"]["gamma"],
                    obsLC = lcs, obsSpec = specs, nbins = mainconfig["nbins"],  gdc_option = "opt", phase1 = mainconfig["phase1"], phase2 = mainconfig["phase2"], oversampling=mainconfig["oversampling"], width=mainconfig["width"])

    system_local.UPDATE_PARAMS(mainconfig["params_init"])

    params_to_opt = array_to_dict(params_array, mainconfig["params_opt"])
    params_all = mainconfig["params_init"].copy()
    params_all.update(params_to_opt)

    system_local = update_system(mainconfig, system_local, params_all, synthV_instance, synthV_imu_instance, obj_type = objective_type)

    logging.debug(f'id {process_id}, objs {system_local.obj_lc + system_local.obj_dis + system_local.obj_synth1 + system_local.obj_synth2}')

    return system_local, system_local.obj_lc + system_local.obj_dis + system_local.obj_synth1 + system_local.obj_synth2,





def get_first_pareto_front(population):
    pareto_fronts = tools.sortNondominated(population, len(population), first_front_only=True)
    return pareto_fronts[0]

def calculate_weights_from_pareto_front(population):
    first_pareto_front = get_first_pareto_front(population)

    objectives = np.array([ind.fitness.values for ind in first_pareto_front])

    min_values = np.min(objectives, axis=0)
    max_values = np.max(objectives, axis=0)

    ranges = max_values - min_values
    inverse_ranges = 1 / np.where(ranges > 0, ranges, 1)

    weights = inverse_ranges / np.sum(inverse_ranges)

    return  weights

def weighted_sum_method(population, weights, n=1):
    if n == 1:
        objectives = np.array([ind.fitness.values for ind in population])
        weighted_sums = np.dot(objectives, weights)
        objective_type = mainconfig.get('objective_type', 'chi2') # chi2 = default value if key is not specified in config file
        if objective_type == 'logP':
            best_solution_index = np.argmax(weighted_sums)
        else:
            best_solution_index = np.argmin(weighted_sums)
        return population[best_solution_index], weighted_sums
    else:
        objectives = np.array([ind.fitness.values for ind in population])
        weighted_sums = np.dot(objectives, weights)
        sorted_indices = np.argsort(weighted_sums)
        if objective_type == 'logP':
            sorted_indices = sorted_indices.reverse()
        best_n_indices = sorted_indices[:n]
        best_n_solutions = [population[i] for i in best_n_indices]
        best_n_weighted_sums = weighted_sums[best_n_indices]

        return best_n_solutions, best_n_weighted_sums


def random_value(bound):
    return np.random.uniform(bound[0], bound[1])


def create_individual_random():
    aa = [random_value(mainconfig["params_bounds"][p]) for p in mainconfig["params_opt"]]
    return aa



def get_logg_range(Teff): # LLmodels atmosphere models grid ranges
    logg_max = 5.0
    logg_min = 4.0
    if 2500 <= Teff <= 5000:
        logg_min = 0.1
    elif 5000 < Teff <= 5500:
        logg_min = 0.5
    elif 5500 < Teff <= 10000:
        logg_min = 2.5
    elif 10000 < Teff <= 11750:
        logg_min = 3.0
    elif 11750 < Teff <= 14000:
        logg_min = 2.0
    elif 14000 < Teff <= 27500:
        logg_min = 3.0
    elif 27500 < Teff <= 30000:
        logg_min = 3.25
    elif 30000 < Teff <= 33000:
        logg_min = 3.5
    return logg_min, logg_max


def estimate_radius_from_teff_logg(teff, logg): # solely for initial grid sampling in physical ranges - atmospheres grid range

    def estimate_mass_from_teff(teff): # solely for selecting range of mass-lum relation
        if teff < 6500:
            # Reference: a general Teff-M relation where M ∝ T_eff^0.5 for stars < 1.5 M_sun
            mass = (teff / 5777)**0.5
        # For M > 1.5 Msun (using the Eker et al. (2018) quadratic relation for high-mass stars)
        elif 6500 <= teff <= 31000:
            log_teff = np.log10(teff)
            # Solve for logM from the quadratic relation: -0.170 * (logM)^2 + 0.888 * logM + 3.671 - logTeff = 0
            a = -0.170
            b = 0.888
            c = 3.671 - log_teff
            discriminant = b**2 - 4 * a * c
            if discriminant < 0:
                return np.nan
            log_m = (-b + np.sqrt(discriminant)) / (2 * a)
            mass = 10**log_m
        else:
            log_teff = np.log10(31000)
            a = -0.170
            b = 0.888
            c = 3.671 - log_teff
            discriminant = b**2 - 4 * a * c
            if discriminant < 0:
                return np.nan
            log_m = (-b + np.sqrt(discriminant)) / (2 * a)
            mass = 10**log_m

        return mass

    def get_k_alpha(M):
        if M < 0.43 :
            k, alpha = 0.23, 2.3  # L ∝ 0.23 * M^2.3
        elif M < 2:
            k, alpha = 1.0, 4.0  # L ∝ M^4.0
        elif M < 55:
            k, alpha = 1.4, 3.5  # L ∝ 1.4 * M^3.5
        else:
            k, alpha = 32000.0, 1.0  # L ∝ 32000 * M

        k_prime = 1 / k**(1 / alpha)
        return k_prime, alpha

    Msun = 1.989e33  # solar mass in grams
    Teff_sun = 5778  # solar temperature in Kelvin
    G = 6.67430e-8  # gravitational constant in cm^3 g^-1 s^-2
    sigma = 5.670374419e-5  # Stefan-Boltzmann constant in erg cm^-2 s^-1 K^-4
    Lsun = 3.828e33  # Solar luminosity in erg/s
    Rsun_cm = 6.96e10  # Solar radius in cm

    M = estimate_mass_from_teff(teff)
    k, alpha = get_k_alpha(M)

    g = 10**logg
    k_in_grams = k * Msun

    if alpha != 1:
        term = G * M * k_in_grams * (4 * np.pi * sigma * teff**4 / Lsun) ** (1 / alpha)
        exponent = 1 / (2 - 2 / alpha)
        R = (term / g) ** exponent
    else:
        R = np.sqrt(G * M*Msun / g)
    R_solarR = R / Rsun_cm
    return R_solarR


def r1_dynamic_bounds(params):
    logg_min, logg_max = get_logg_range(params['Teff1'])
    rmin_atm = estimate_radius_from_teff_logg(params['Teff1'], logg_max) / params['a']
    rmax_atm = estimate_radius_from_teff_logg(params['Teff1'], logg_min) / params['a']

    rmax_roche = 0.49 * params['q']**(2.0/3.0) / ( 0.6 * params['q']**(2.0/3.0) + np.log(1+params['q']**(1.0/3.0)))
    ecc = params['f_s']**2 + params['f_c']**2
    rmax_roche *= (1-ecc)

    print(f'atm rmin, rmax: {rmin_atm}, {rmax_atm}')
    print(f'roche rmax: {rmax_roche}')

    rmin = max([rmin_atm, mainconfig['params_bounds']['r1'][0]])
    rmax = min([rmax_atm, rmax_roche, mainconfig['params_bounds']['r1'][1]])

    return [rmin, rmax]

def r2_dynamic_bounds(params):
    logg_min, logg_max = get_logg_range(params['Teff2'])
    rmin_atm = estimate_radius_from_teff_logg(params['Teff2'], logg_max) / params['a']
    rmax_atm = estimate_radius_from_teff_logg(params['Teff2'], logg_min) / params['a']

    rmax_roche = 0.49 * (1/params['q'])**(2.0/3.0) / ( 0.6 * (1/params['q'])**(2.0/3.0) + np.log(1+(1/params['q'])**(1.0/3.0)))
    ecc = params['f_s']**2 + params['f_c']**2
    rmax_roche *= (1-ecc)

    print(f'atm rmin, rmax: {rmin_atm}, {rmax_atm}')
    print(f'roche rmax: {rmax_roche}')

    rmin = max([rmin_atm, mainconfig['params_bounds']['r2'][0]])
    rmax = min([rmax_atm, rmax_roche, mainconfig['params_bounds']['r2'][1]])

    return [rmin, rmax]

def incl_dynamic_bounds(params):

    sin_i_min = params['r1'] + params['r2']
    if sin_i_min >= 1:
        i_min = 89.0
    else:
        i_min = math.degrees(math.asin(sin_i_min))

    incl_min = max([i_min, mainconfig['params_bounds']['incl'][0] ])
    incl_max = min([90.0, mainconfig['params_bounds']['incl'][1] ])

    return [incl_min, incl_max]

def scale_sample(sample, bounds):
    return bounds[0] + (bounds[1] - bounds[0]) * sample

def scale_sample_a(sample, a_bounds): # more uniform sampling in mass
    a_min, a_max = a_bounds
    a_scaled = ((a_max**3 - a_min**3) * sample + a_min**3) ** (1 / 3)
    return a_scaled

def create_individual():

    sample = sobol_engine.random()[0]
    params_sampled = {}

    for i, param in enumerate(mainconfig["params_opt"]):
        if param not in ["r1", "r2", "incl"]:
            bounds = mainconfig["params_bounds"][param]
            if param == "a":
                sampled_value = scale_sample_a(sample[i], bounds)
                params_sampled[param] = sampled_value
            else:
                sampled_value = scale_sample(sample[i], bounds)
                params_sampled[param] = sampled_value

    for param, default_value in mainconfig["params_init"].items():
        if param not in params_sampled:
            params_sampled[param] = default_value

    for param in ["r1", "r2", "incl"]:
        if param in mainconfig["params_opt"]:
            i = mainconfig["params_opt"].index(param)
            if param == "r1":
                bounds = r1_dynamic_bounds(params_sampled)
            elif param == "r2":
                bounds = r2_dynamic_bounds(params_sampled)
            elif param == "incl":
                bounds = incl_dynamic_bounds(params_sampled)

            print(param)
            print(bounds)

            sampled_value = scale_sample(sample[i], bounds)
            params_sampled[param] = sampled_value

    individual = [params_sampled[param] for param in mainconfig["params_opt"]]

    return individual



def create_individual_uni():
    sample = sobol_engine.random()[0]
    individual = []
    for i, param in enumerate(mainconfig["params_opt"]):
        bounds = mainconfig["params_bounds"][param]
        individual.append(scale_sample(sample[i], bounds))
    return individual




def evaluate(individual, mainconfig=None, ind_id=None):

    system_local, obj_local = objective(mainconfig, individual, ind_id) if ind_id is not None else objective(mainconfig, individual, 1)

    updated_params = {
        'r1' : system_local.primary.radius_a,
        'r2' : system_local.secondary.radius_a,
        'q' : system_local.secondary.mass_M1,
        'm2' : system_local.secondary.mass_M1,
        'Teff1' : system_local.primary.Teff,
        'Teff2' : system_local.secondary.Teff,
        'metallicity1' : system_local.primary.metallicity,
        'metallicity2' : system_local.secondary.metallicity,
        'Ve1' : system_local.primary.v_equatorial,
        'Ve2' : system_local.secondary.v_equatorial,
        'Ve1sini' : system_local.primary.vsini,
        'Ve2sini' : system_local.secondary.vsini,
        'heat_1' : system_local.primary.heated,
        'heat_2' : system_local.secondary.heated,
        'gdc_1' : system_local.primary.gdc,
        'gdc_2' : system_local.secondary.gdc,

        'bfac_1' : system_local.primary.bfac,
        'bfac_2' : system_local.secondary.bfac,

        'incl' : system_local.orbit.inc_deg,
        'porb' : system_local.orbit.porb,
        't0' : system_local.orbit.t0,
        'a' : system_local.orbit.a,
        'f_c' : system_local.orbit.f_c,
        'f_s' : system_local.orbit.f_s,

        'l3' : system_local.l3,
        'gamma' : system_local.gamma,
    }


    for i, param in enumerate(mainconfig["params_opt"]):
        individual[i] = updated_params[param]

    return system_local.obj_lc, system_local.obj_dis, system_local.obj_synth1, system_local.obj_synth2,


def evaluate_chunk(chunk):
    async_results = [pool.apply_async(toolbox.evaluate, args=(ind, mainconfig, ind_id)) for ind_id, ind in enumerate(chunk, start=1)]

    objective_type = mainconfig.get('objective_type', 'chi2') # chi2 = default value if key is not specified in config file

    results = []
    for async_result in async_results:
        try:
            result = async_result.get(timeout=120)
            if mainconfig["ifsecondspec"] == 1 :
                results.append(result)
            else:
                results.append(result[:-1])
        except TimeoutError:
            logging.warning("!!! TIMEOUT for one set of parameters !!!")
            if objective_type == "logP":
                if mainconfig["ifsecondspec"] == 1:
                    results.append((-1.0e10, -1.0e10, -1.0e10, -1.0e10))
                else:
                    results.append((-1.0e10, -1.0e10, -1.0e10))
            else:
                if mainconfig["ifsecondspec"] == 1:
                    results.append((1.0e3, 1.0e3, 1.0e3, 1.0e3))
                else:
                    results.append((1.0e3, 1.0e3, 1.0e3))
        except Exception as e:
            logging.error(e)
            if objective_type == "logP":
                if mainconfig["ifsecondspec"] == 1:
                    results.append((-1.0e10, -1.0e10, -1.0e10, -1.0e10))
                else:
                    results.append((-1.0e10, -1.0e10, -1.0e10))
            else:
                if mainconfig["ifsecondspec"] == 1:
                    results.append((1.0e3, 1.0e3, 1.0e3, 1.0e3))
                else:
                    results.append((1.0e3, 1.0e3, 1.0e3))

    return results



def custom_map(evaluate_function, individuals):
    chunk_size = mainconfig["nthreads"]
    chunks = [individuals[i:i + chunk_size] for i in range(0, len(individuals), chunk_size)]

    all_results = []
    for chunk in chunks:
        all_results.extend(evaluate_chunk(chunk))

    return all_results


def save_checkpoint(generation, population, hall_of_fame, logbook, filename="checkpoint.pkl"):
    with open(filename, "wb") as cp_file:
        pickle.dump({
            "generation": generation,
            "population": population,
            "halloffame": hall_of_fame,
            "logbook": logbook
        }, cp_file)

def load_checkpoint(filename="checkpoint.pkl"):
    if not os.path.isfile(filename):
        return None
    with open(filename, "rb") as cp_file:
        content = pickle.load(cp_file)
    return content



def run_ga_checkpoints(ngen=20, popsize=50, cxpb=0.5, mutpb=0.2, mutpb_stagnated=0.5, chunk_size=8, checkpoint_interval=2, verbose = True, plot_every_best = False, mass_extinction_threshold=10):

    checkpoint_file = "checkpoint_ngen{}_popsize{}_initPorb{}_initA{}.pkl".format(ngen, popsize, mainconfig["params_init"]["porb"], mainconfig["params_init"]["a"])

    hof = tools.ParetoFront()
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.median, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)
    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])


    df = pd.DataFrame([mainconfig["params_init"]])
    df.to_csv(f'{mainconfig["saveto"]}/init_state_initPorb{mainconfig["params_init"]["porb"]}_initA{mainconfig["params_init"]["a"]}_ngen{ngen}_popsize{popsize}.dat', sep="\t", index = False)

    if os.path.isfile(checkpoint_file):
        with open(checkpoint_file, "rb") as cp_file:
            cp = pickle.load(cp_file)
            pop, start_gen, hof, logbook, rndstate = cp["population"], cp["generation"], cp["halloffame"], cp["logbook"], cp["rndstate"]
            random.setstate(rndstate)
    else:
        start_gen = 0
        pop = toolbox.population(n = int(popsize * 1.1))
        chunks = [pop[i:i + chunk_size] for i in range(0, len(pop), chunk_size)]
        fitness_results = []
        for chunk in chunks:
            fitness_results.extend(evaluate_chunk(chunk))
        for ind, fit in zip(pop, fitness_results):
            ind.fitness.values = fit

    def has_pareto_front_stagnated(previous_front, current_front, tolerance=0.01):
        if len(current_front) != len(previous_front):
            return False

        changes = [euclidean(ind1.fitness.values, ind2.fitness.values)
                   for ind1, ind2 in zip(previous_front, current_front)]
        average_change = sum(changes) / len(changes)
        return average_change < tolerance

    previous_hof = None
    previous_best = 1e5
    stagnation_counter = 0

    def calculate_weighted_sum(individual, weights):
        return sum(w * val for w, val in zip(weights, individual.fitness.values))

    ref_points = tools.uniform_reference_points(4, 5)
    weights = [1e-5,1e-5,1e-5,1e-5]

    for gen in range(start_gen, ngen):

        logging.info(f'GENERATION {gen} STARTED')

        if gen % checkpoint_interval == 0:
            cp = dict(population=pop, generation=gen, halloffame=hof, logbook=logbook, rndstate=random.getstate())
            with open(checkpoint_file, "wb") as cp_file:
                pickle.dump(cp, cp_file)

            last_gen_parameters = []
            last_gen_fitness = []
            for ind in pop:
                params_to_optimize = dict(zip([param for param in mainconfig["params_opt"] if param in mainconfig["params_bounds"]], ind))
                last_gen_parameters.append(params_to_optimize)
                last_gen_fitness.append(ind.fitness.values)

            df_params = pd.DataFrame(last_gen_parameters)
            df_fitness = pd.DataFrame(last_gen_fitness, columns=[f'obj{i+1}' for i in range(len(last_gen_fitness[0]))])
            df = pd.concat([df_params, df_fitness], axis=1)
            df.to_csv(f'{mainconfig["saveto"]}/gen{gen}_initPorb{mainconfig["params_init"]["porb"]}_initA{mainconfig["params_init"]["a"]}_ngen{ngen}_popsize{popsize}.dat', sep="\t")


        if stagnation_counter <= 5 :
            offspring = algorithms.varAnd(pop, toolbox, cxpb, mutpb)
        elif 5 < stagnation_counter < mass_extinction_threshold :
            offspring = algorithms.varAnd(pop, toolbox, cxpb, mutpb_stagnated)
            logging.info(f"Stagnation: mutpb increased")
        else:
            # Trigger mass extinction
            num_replace = int(len(pop) * 0.5)
            new_individuals = toolbox.population(n=num_replace)
            pop[-num_replace:] = new_individuals
            offspring = algorithms.varAnd(pop, toolbox, cxpb, mutpb)
            stagnation_counter = 0
            logging.info(f"Stagnation: mass extinction triggered")

        fits = toolbox.map(toolbox.evaluate, offspring)
        for fit, ind in zip(fits, offspring):
            ind.fitness.values = fit

        pop[:] = tools.selNSGA2(offspring + pop, popsize)

        record = stats.compile(pop)
        logbook.record(gen=gen, nevals=len(offspring), **record)

        weights = calculate_weights_from_pareto_front(pop)
        logging.info(f"Calculated Weights: {weights}")
        best_solution, weighted_sums = weighted_sum_method(pop, weights)
        logging.info(f"Objectives: {best_solution.fitness.values}")

        hof.update(pop)
        #hof = sorted(hof, key=lambda ind: calculate_weighted_sum(ind, weights))

        best = calculate_weighted_sum(best_solution, weights)

        if previous_hof and has_pareto_front_stagnated(previous_hof, hof):
            stagnation_counter += 1
        elif np.abs(best-previous_best) < 1e-8:
            stagnation_counter += 1
        else:
            stagnation_counter = 0

        previous_hof = list(hof)
        previous_best = calculate_weighted_sum(hof[0], weights)


        best_params = dict(zip([param for param in mainconfig["params_opt"] if param in mainconfig["params_bounds"]], best_solution))
        full_best_params = mainconfig["params_init"].copy()
        full_best_params.update(best_params)

        logging.info("Best parameters:")
        for key, value in full_best_params.items():
            logging.info(f"{key}: {value}")

        if plot_every_best:
            transformed_params = {
                "sbratio" : 1.0,
                }

            orb = orbit(porb = mainconfig["params_init"]["porb"], inc_deg = mainconfig["params_init"]["incl"], t0 = mainconfig["params_init"]["t0"], a = mainconfig["params_init"]["a"], q = mainconfig["params_init"]["q"], f_c = mainconfig["params_init"]["f_c"], f_s = mainconfig["params_init"]["f_s"])

            system_obs = binary(orbit = orb, primary = star(), secondary = star(), sbratio = transformed_params["sbratio"], l3 = mainconfig["params_init"]["l3"], gamma = mainconfig["params_init"]["gamma"],
                            obsLC = lcs, obsSpec = specs, nbins = mainconfig["nbins"],  gdc_option = "opt", phase1 = mainconfig["phase1"], phase2 = mainconfig["phase2"], oversampling=mainconfig["oversampling"], width=mainconfig["width"])

            system_obs.UPDATE_PARAMS(mainconfig["params_init"])

            process_id = current_process()._identity[0] if current_process()._identity else 1

            synthV_instance = synthini_s[process_id - 1]
            synthV_imu_instance = synthini_simu[process_id - 1]

            objective_type = mainconfig.get('objective_type', 'chi2') # chi2 = default value if key is not specified in config file

            system_obs = update_system(mainconfig, system_obs, full_best_params, synthV_instance, synthV_imu_instance, ifplot_main=True, titleprefix=f'gen{gen}___', obj_type = objective_type)

        logging.info(f'GENERATION {gen} FINISHED')

    return pop, logbook, hof



def run_grid(popsize=50, chunk_size=8, verbose = True, plot_every_best = False, file_out = "default.file"):

    npars = len(mainconfig["params_opt"])
    pop = toolbox.population(n=popsize)
    chunks = [pop[i:i + chunk_size] for i in range(0, len(pop), chunk_size)]
    nchunk = 0
    for chunk in chunks:
        fitness_results = evaluate_chunk(chunk)

        with open(file_out,"a") as out:
            for ia in range(chunk_size):
                fstr = ""
                for ib in range(npars):
                    fstr += f'{chunk[ia][ib]}\t'
                for ib in range(4):
                    fstr += f'{fitness_results[ia][ib]}\t'
                fstr += '\n'
                out.write(fstr)
        nchunk += 1
        print(f'{nchunk*chunk_size} out of {popsize}')




def main():

    objective_type = mainconfig.get('objective_type', 'chi2') # chi2 = default value if key is not specified in config file
    print(objective_type)

    transformed_params = {
        "sbratio" : 1.0,
        }

    orb = orbit(porb = mainconfig["params_init"]["porb"], inc_deg = mainconfig["params_init"]["incl"], t0 = mainconfig["params_init"]["t0"], a = mainconfig["params_init"]["a"], q = mainconfig["params_init"]["q"], f_c = mainconfig["params_init"]["f_c"], f_s = mainconfig["params_init"]["f_s"])

    system_obs = binary(orbit = orb, primary = star(), secondary = star(), sbratio = transformed_params["sbratio"], l3 = mainconfig["params_init"]["l3"], gamma = mainconfig["params_init"]["gamma"],
                    obsLC = lcs, obsSpec = specs, nbins = mainconfig["nbins"],  gdc_option = "opt", phase1 = mainconfig["phase1"], phase2 = mainconfig["phase2"], oversampling=mainconfig["oversampling"], width=mainconfig["width"])

    system_obs.UPDATE_PARAMS(mainconfig["params_init"])

    global pool
    pool = Pool(processes=int(mainconfig["nthreads"]))

    process_id = current_process()._identity[0] if current_process()._identity else 1

    synthV_instance = synthini_s[process_id - 1]
    synthV_imu_instance = synthini_simu[process_id - 1]

    system_obs = update_system(mainconfig,system_obs, mainconfig["params_init"], synthV_instance, synthV_imu_instance, ifplot_main=True, if_save_mods = False, obj_type = objective_type)

    if objective_type == 'logP':
        objweights=(1.0, 1.0, 1.0)
        if mainconfig["ifsecondspec"] == 1:
            objweights=(1.0, 1.0, 1.0, 1.0)
        creator.create("FitnessMax", base.Fitness, weights=objweights)
        creator.create("Individual", list, fitness=creator.FitnessMax)
    else:
        objweights=(-1.0, -1.0, -1.0)
        if mainconfig["ifsecondspec"] == 1:
            objweights=(-1.0, -1.0, -1.0, -1.0)
        creator.create("FitnessMin", base.Fitness, weights=objweights)
        creator.create("Individual", list, fitness=creator.FitnessMin)

    global toolbox
    toolbox = base.Toolbox()
    toolbox.register("individual", tools.initIterate, creator.Individual, create_individual)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("evaluate", evaluate)
    toolbox.register("mate", tools.cxBlend, alpha=mainconfig["cxblend_a"])
    toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=1, indpb=0.1)
    toolbox.register("map", custom_map)
    toolbox.register("select", tools.selNSGA2)


    pop, stats, hof = run_ga_checkpoints(
        ngen = mainconfig["ngen"],
        popsize = mainconfig["popsize"] ,
        cxpb= mainconfig["cxpb"],
        mutpb= mainconfig["mutpb"],
        mutpb_stagnated=mainconfig["mutpb_stagnated"],
        chunk_size= mainconfig["nthreads"],
        checkpoint_interval=mainconfig["checkpoint_interval"],
        verbose = False,
        plot_every_best=False,
        mass_extinction_threshold=mainconfig["mass_extinction_threshold"]
        )


    pareto_weights = calculate_weights_from_pareto_front(pop)

    best_solution, weighted_sums = weighted_sum_method(pop, pareto_weights)

    evaluate(best_solution, mainconfig=mainconfig)

    best_params = dict(zip([param for param in mainconfig["params_opt"] if param in mainconfig["params_bounds"]], best_solution))
    full_best_params = mainconfig["params_init"].copy()
    full_best_params.update(best_params)

    logging.critical("Best parameters:")
    for key, value in full_best_params.items():
        logging.critical(f"{key}: {value}")

    system_obs = update_system(mainconfig,system_obs, full_best_params,
                               synthV_instance, synthV_imu_instance, ifplot_main=True, obj_type = objective_type)
    print("!")
    print(system_obs.primary.mass)
    print("!")

    pool.close()
    pool.join()

if __name__ == '__main__':
    main()
