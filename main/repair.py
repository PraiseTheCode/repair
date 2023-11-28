import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import chebyshev

import random

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
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
    if os.path.exists(f'SynthV_core{idir}'):
        shutil.rmtree(f'SynthV_core{idir}', ignore_errors=True)
    #os.mkdir(f'SynthV_core{idir}')
    shutil.copytree(mainconfig["synthVpath"], f'SynthV_core{idir}')
    
    if os.path.exists(f'convolve_core{idir}'):
        shutil.rmtree(f'convolve_core{idir}', ignore_errors=True)
    #os.mkdir(f'convolve_core{idir}')
    shutil.copytree(mainconfig["convolvepath"], f'convolve_core{idir}')
    
    synthini_s.append(synthetic(wd_synthv=f'SynthV_core{idir}', wd_convolve=f'convolve_core{idir}', db_atmos_models = mainconfig["atmmodels"], if_convolve=True, if_imu=False, pref = os.getcwd()))
    synthini_simu.append(synthetic(wd_synthv=f'SynthV_core{idir}', db_atmos_models = mainconfig["atmmodels"], if_convolve=False, if_imu=True, if_lines = False, pref = os.getcwd()))



sobol_engine = Sobol(d=len(mainconfig["params_opt"]), scramble=True)

specs = specSeries()
if mainconfig["obsmode"] == "test":
    specs.LOAD_TESTSUIT(mainconfig["path_specs"], period = mainconfig["params_init"]["porb"], t0 = mainconfig["params_init"]["t0"])
else:
    specs.LOAD_HERMES(mainconfig["path_specs"])


specs.TRIM(mainconfig["wl11"]+10, mainconfig["wl21"]-10)



lcf = np.loadtxt(mainconfig["path_lc"])
if mainconfig["obsmode"] == "test":
    lcs = lightcurve(lcf[:,0], lcf[:,1], np.zeros_like(lcf[:,0]))
else:
    lcs = lightcurve(lcf[:,0], lcf[:,1], lcf[:,2])

passband = np.loadtxt("tess_transmission.txt", skiprows=6,delimiter=",")
passband[:,0] *= 10.0

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
    spectrum_interp = interp1d(wavelengths, flux, bounds_error=False, fill_value=0.0)

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



def update_system(system, params, synthV_instance, synthV_imu_instance, ifplot_main = False, if_save_mods = False, if_sbratio_spec = True):
    
    system.UPDATE_PARAMS(params)
    
    logging.info(f'Primary Teff {system.primary.Teff}, logg {system.primary.logg}, vsini {system.primary.vsini}')
    logging.info(f'Secondary Teff {system.secondary.Teff}, logg {system.secondary.logg}, vsini {system.secondary.vsini}')


    b_prim, b_sec = 1.0, 1.0
    
    try:
        system.primary.spectrum, m_prim, im_prim, system.primary.grid_logg  = synthV_instance.COMPUTE(system.primary.Teff, system.primary.logg, vsini=system.primary.vsini, vmacro = mainconfig["vmacro1"], wl1=mainconfig["wl11"], wl2=mainconfig["wl21"], wlstep = mainconfig["wlstep1"])
        nosp, mu, imu, grid_logg = synthV_imu_instance.COMPUTE(system.primary.Teff, system.primary.logg, wl1=5000, wl2=11000, wlstep = 5)
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
        logging.error("!!exception Synth 1!!")
        logging.error(e)
        system.chi2_synth1 = 1e3
        system.chi2_synth2 = 1e3
        system.chi2_lc = 1e3
        system.chi2_dis = 1e3
        return system
    
    try:
        system.secondary.spectrum, m_sec, im_sec, system.secondary.grid_logg  = synthV_instance.COMPUTE(system.secondary.Teff, system.secondary.logg, vsini=system.secondary.vsini, vmacro = mainconfig["vmacro2"], wl1=mainconfig["wl12"], wl2=mainconfig["wl22"], wlstep = mainconfig["wlstep2"])
        nosp, mu, imu, grid_logg = synthV_imu_instance.COMPUTE(system.secondary.Teff, system.secondary.logg, wl1=5000, wl2=11000, wlstep = 5)
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
        logging.error("!!exception Synth 2!!")
        logging.error(e)
        system.chi2_synth1 = 1e3
        system.chi2_synth2 = 1e3
        system.chi2_lc = 1e3
        system.chi2_dis = 1e3
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
        system.chi2_lc = 1e3
        system.chi2_dis = 1e3
        system.chi2_synth1 = 1e3
        system.chi2_synth2 = 1e3
        return system
    
    residuals_lc = system.residLC
    residuals_lc = residuals_lc[np.isfinite(residuals_lc)]
    schi2lc = np.sum(residuals_lc**2) / float(len(residuals_lc))
    
    
    system.chi2_lc = schi2lc #* 2e4
    
    
    if system.primary.sync_fact > 1.0 :
        system.chi2_lc *= np.abs(1.0 - system.primary.sync_fact)**0.5
    if system.secondary.sync_fact > 1.0 :
        system.chi2_lc *= np.abs(1.0 - system.secondary.sync_fact)**0.5
    logging.debug(f'sync1 {system.primary.sync_fact}')
    logging.debug(f'sync2 {system.secondary.sync_fact}')
    logging.debug(f'corr1 {np.abs(1.0 - system.primary.sync_fact)**0.5}')
    
    
    try:
        system.UPDATE_SP_SEPARATION(ifplot=False, if_lightfact = True)
    except Exception as e:
        logging.error("!!exception DIS!!")
        logging.error(e)
        system.chi2_lc = 1e3
        system.chi2_dis = 1e3
        system.chi2_synth1 = 1e3
        system.chi2_synth2 = 1e3
        return system
    
        
    residuals_sp = system.residSpec
    chi2_dis = 0
    for i in range(len(residuals_sp[0])):
        aa = np.sum(residuals_sp[:,i]**2 ) / float(len(residuals_sp[:,i]))
        chi2_dis += aa
        
    system.chi2_dis = chi2_dis
    
    schi2synth = 0
    ffi = interp1d(system.primary.spectrum[:,0], system.primary.spectrum[:,1], fill_value="extrapolate")
    mm = ffi(system.obsSpec.wl_eqlog_to_wl)
    ob1 = fetch_spec_new(system.obsSpec.wl_eqlog_to_wl, system.obsSpec.mod[0] + 0.5, system.obsSpec.wl_eqlog_to_wl, mm) #, 500, 3)
    schi2synth1 = sum( (ob1 - mm)**2 ) / float(len(mm))
    ffi = interp1d(system.secondary.spectrum[:,0], system.secondary.spectrum[:,1], fill_value="extrapolate")
    mm = ffi(system.obsSpec.wl_eqlog_to_wl)
    ob2 = fetch_spec_new(system.obsSpec.wl_eqlog_to_wl, system.obsSpec.mod[1] + 0.5, system.obsSpec.wl_eqlog_to_wl, mm) #, 500, 3)
    schi2synth2 = sum( (ob2 - mm)**2 ) / float(len(mm))
    #schi2synth = schi2synth1 + schi2synth2 
    
    #system.chi2_synth = schi2synth #* 3e3
    system.chi2_synth1 = schi2synth1
    system.chi2_synth2 = schi2synth2
    
    if np.abs(system.primary.grid_logg - system.primary.logg) > 0.5:
        system.chi2_synth1 *= 100.0
    if np.abs(system.secondary.grid_logg - system.secondary.logg) > 0.5:
        system.chi2_synth2 *= 100.0
    logging.debug(f'logg1 {system.primary.logg}, grid {system.primary.grid_logg}')
    logging.debug(f'logg2 {system.secondary.logg}, grid {system.secondary.grid_logg}')
        
    logging.debug(f'chi2s: {system.chi2_lc} {system.chi2_dis} {system.chi2_synth1} {system.chi2_synth2}')
    
    if ifplot_main:
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
        
        shift = (np.average(np.abs(residuals_sp[:,0])) + np.average(np.abs(residuals_sp[:,-1]))) / 5 #*2 #5
        for i in range(len(residuals_sp[0])):
            ax_dis.plot(residuals_sp[:,i] - shift*i, lw=0.5, c="crimson", alpha=0.3)
            
        ax_dis.set_ylim([-(len(residuals_sp[0])+1)*shift, shift])
        ax_dis.set_xlim([100,len(residuals_sp[:,0])-100])
        
        ax1.plot(system.obsSpec.wl_eqlog_to_wl, ob1, lw=1,c='cornflowerblue')
        ax1.plot(system.primary.spectrum[:,0], system.primary.spectrum[:,1], lw=1, c="crimson")
        ax2.plot(system.obsSpec.wl_eqlog_to_wl, ob2, lw=1,c='cornflowerblue')
        ax2.plot(system.secondary.spectrum[:,0], system.secondary.spectrum[:,1], lw=1, c="crimson")
        ax1.set_xlim([mainconfig["wl11"], mainconfig["wl21"]])
        ax2.set_xlim([mainconfig["wl11"],mainconfig["wl21"]])
        ax3.plot(system.obsSpec.wl_eqlog_to_wl, ob1, lw=1,c='cornflowerblue')
        ax3.plot(system.primary.spectrum[:,0], system.primary.spectrum[:,1],  lw=1, c="crimson")
        ax4.plot(system.obsSpec.wl_eqlog_to_wl, ob2, lw=1,c='cornflowerblue')
        ax4.plot(system.secondary.spectrum[:,0], system.secondary.spectrum[:,1], lw=1, c="crimson")
        ax3.set_xlim([4300,4380])
        ax4.set_xlim([4300,4380])
        ax5.plot(system.obsSpec.wl_eqlog_to_wl, ob1, lw=1,c='cornflowerblue')
        ax5.plot(system.primary.spectrum[:,0], system.primary.spectrum[:,1], lw=1, c="crimson")
        ax6.plot(system.obsSpec.wl_eqlog_to_wl, ob2, lw=1,c='cornflowerblue')
        ax6.plot(system.secondary.spectrum[:,0], system.secondary.spectrum[:,1], lw=1, c="crimson")
        ax5.set_xlim([4450,4500])
        ax6.set_xlim([4450,4500])
        
        ax_lc.text(1.25,0.95,str(np.round(np.log10(system.chi2_lc)+10.0,3)))
        ax_dis.text(4200,0.0-shift,str(np.round(np.log10(system.chi2_dis)+10.0,3)))
        
        ax1.text(4450,0.8,str(np.round(np.log10(system.chi2_synth1)+10.0,3)))
        ax2.text(4450,0.8,str(np.round(np.log10(system.chi2_synth2)+10.0,3)))
        ax3.text(4370,0.8,str(np.round(np.log10(system.chi2_synth1)+10.0,3)))
        ax4.text(4370,0.8,str(np.round(np.log10(system.chi2_synth2)+10.0,3)))
        ax5.text(4490,0.8,str(np.round(np.log10(system.chi2_synth1)+10.0,3)))
        ax6.text(4490,0.8,str(np.round(np.log10(system.chi2_synth2)+10.0,3)))
        
        
        plt.savefig(f'{mainconfig["saveto"]}/Teff1_{system.primary.Teff:.1f}_Teff2_{system.secondary.Teff:.1f}_logg1_{system.primary.logg:.2f}_logg2_{system.secondary.logg:.2f}.png',dpi=150)
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

def chi2(params_array, process_id):
    
    
    synthV_instance = synthini_s[process_id - 1]
    synthV_imu_instance = synthini_simu[process_id - 1]
    
    transformed_params = {
        "q" : mainconfig["params_init"]['m2'], 
        "a" : mainconfig["params_init"]['asr'] * mainconfig["params_init"]['R1'] * (1 + mainconfig["params_init"]['r2']), 
        "sbratio" : 1.0,
        }

    orb = orbit(porb = mainconfig["params_init"]["porb"], inc_deg = mainconfig["params_init"]["incl"], t0 = mainconfig["params_init"]["t0"], a = transformed_params["a"], q = transformed_params["q"], f_c = mainconfig["params_init"]["f_c"], f_s = mainconfig["params_init"]["f_s"])
    
    system_local = binary(orbit = orb, primary = star(), secondary = star(), sbratio = transformed_params["sbratio"], l3 = mainconfig["params_init"]["l3"], gamma = mainconfig["params_init"]["gamma"],
                    obsLC = lcs, obsSpec = specs, nbins = mainconfig["nbins"],  gdc_option = "opt", phase1 = mainconfig["phase1"], phase2 = mainconfig["phase2"], oversampling=mainconfig["oversampling"], width=mainconfig["width"]) 
    
    system_local.UPDATE_PARAMS(mainconfig["params_init"])
    
    params_to_opt = array_to_dict(params_array, mainconfig["params_opt"]) 
    params_all = mainconfig["params_init"].copy()
    params_all.update(params_to_opt)
    
    system_local = update_system(system_local, params_all, synthV_instance, synthV_imu_instance)
    
    logging.debug(f'id {process_id}, chi2 {np.sqrt(system_local.chi2_lc**2 + system_local.chi2_dis**2 + system_local.chi2_synth1**2 + system_local.chi2_synth2**2)}')
    
   #return system_local, (np.log10(system_local.chi2_lc)+10.0 + np.log10(system_local.chi2_dis)+10.0 + np.log10(system_local.chi2_synth1)+10.0 + np.log10(system_local.chi2_synth2)+10.0),
    return system_local, np.sqrt(system_local.chi2_lc**2 + system_local.chi2_dis**2 + system_local.chi2_synth1**2 + system_local.chi2_synth2**2),




def plot_pareto_front(population, best_solution=None):
    objectives = np.array([ind.fitness.values for ind in population])

    filter_mask = (objectives < 0.1e3).all(axis=1)
    #filter_mask = (objectives < 10.0).all(axis=1)
    filtered_objectives = objectives[filter_mask]

    if best_solution is not None:
        best_objectives = np.array(best_solution.fitness.values)
        

    fig, axs = plt.subplots(1, 4, figsize=(12, 4))

    # LC vs Dis
    axs[0].scatter(filtered_objectives[:, 0], filtered_objectives[:, 1], alpha=0.6)
    if best_solution is not None:
        axs[0].scatter(best_objectives[0], best_objectives[1], color='red', marker='x', s=100)
    axs[0].set_xlabel('LC')
    axs[0].set_ylabel('Dis')
    axs[0].set_title(r'LC $\chi^2$ vs Dis $\chi^2$')

    # LC vs Synth1
    axs[1].scatter(filtered_objectives[:, 0], filtered_objectives[:, 2], alpha=0.6)
    if best_solution is not None:
        axs[1].scatter(best_objectives[0], best_objectives[2], color='red', marker='x', s=100)
    axs[1].set_xlabel('LC')
    axs[1].set_ylabel('Synth1')
    axs[1].set_title(r'LC $\chi^2$ vs Synth1 $\chi^2$')

    # Dis vs Synth1
    axs[2].scatter(filtered_objectives[:, 1], filtered_objectives[:, 2], alpha=0.6)
    if best_solution is not None:
        axs[2].scatter(best_objectives[1], best_objectives[2], color='red', marker='x', s=100)
    axs[2].set_xlabel('Dis')
    axs[2].set_ylabel('Synth1')
    axs[2].set_title(r'Dis $\chi^2$ vs Synth1 $\chi^2$')

    # Synth1 vs Synth2
    axs[3].scatter(filtered_objectives[:, 2], filtered_objectives[:, 3], alpha=0.6)
    if best_solution is not None:
        axs[3].scatter(best_objectives[2], best_objectives[3], color='red', marker='x', s=100)
    axs[3].set_xlabel('Synth1')
    axs[3].set_ylabel('Synth2')
    axs[3].set_title(r'Synth1 $\chi^2$ vs Synth2 $\chi^2$')

    plt.tight_layout()
    plt.show()


def calculate_weights_from_pareto_front(population):
    objectives = np.array([ind.fitness.values for ind in population])
    min_values = np.min(objectives, axis=0)
    max_values = np.max(objectives, axis=0)
    ranges = max_values - min_values
    inverse_ranges = 1 / np.where(ranges > 0, ranges, 1)  
    weights = inverse_ranges / np.sum(inverse_ranges)
    
    #weights = np.ones_like(weights) 
    
    return weights

def weighted_sum_method(population, weights):
    objectives = np.array([ind.fitness.values for ind in population])
    weighted_sums = np.dot(objectives, weights)
    best_solution_index = np.argmin(weighted_sums)
    return population[best_solution_index], weighted_sums

def random_value(bound):
    return np.random.uniform(bound[0], bound[1])


def create_individual_random():
    aa = [random_value(mainconfig["params_bounds"][p]) for p in mainconfig["params_opt"]]
    return aa

def scale_sample(sample, bounds):
    return bounds[0] + (bounds[1] - bounds[0]) * sample

def create_individual():
    sample = sobol_engine.random()[0]  
    individual = []
    for i, param in enumerate(mainconfig["params_opt"]):
        bounds = mainconfig["params_bounds"][param]
        individual.append(scale_sample(sample[i], bounds))
    return individual




def evaluate(individual, ind_id=None):
    
    system_local, chi2_local = chi2(individual, ind_id) if ind_id is not None else chi2(individual, 1)
    
    updated_params = {
        'R1' : system_local.primary.radius_sol, 
        'r2' : system_local.secondary.radius_R1, 
        'M1' : system_local.primary.mass, 
        'm2' : system_local.secondary.mass_M1,
        'Teff1' : system_local.primary.Teff, 
        'Teff2' : system_local.secondary.Teff,
        'Ve1' : system_local.primary.v_equatorial, 
        'Ve2' : system_local.secondary.v_equatorial,
        'heat_1' : system_local.primary.heated, 
        'heat_2' : system_local.secondary.heated,
        'gdc_1' : system_local.primary.gdc, 
        'gdc_2' : system_local.secondary.gdc,
        
        'bfac_1' : system_local.primary.bfac, 
        'bfac_2' : system_local.secondary.bfac,
        
        'incl' : system_local.orbit.inc_deg, 
        'porb' : system_local.orbit.porb, 
        't0' : system_local.orbit.t0,
        'asr' : system_local.orbit.a / (system_local.primary.radius_sol + system_local.secondary.radius_sol),
        'f_c' : system_local.orbit.f_c, 
        'f_s' : system_local.orbit.f_s,
         
        'l3' : system_local.l3,
        'gamma' : system_local.gamma,
    }
    
    
    for i, param in enumerate(mainconfig["params_opt"]):
        individual[i] = updated_params[param]
  
        return system_local.chi2_lc, system_local.chi2_dis, system_local.chi2_synth1, system_local.chi2_synth2,
    

def evaluate_chunk(chunk):
    async_results = [pool.apply_async(toolbox.evaluate, args=(ind, ind_id)) for ind_id, ind in enumerate(chunk, start=1)]

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
            if mainconfig["ifsecondspec"] == 1:
                results.append((1.0e3, 1.0e3, 1.0e3, 1.0e3))
            else:
                results.append((1.0e3, 1.0e3, 1.0e3))
        except Exception as e:
            logging.error(e)
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



def run_ga_checkpoints(ngen=20, popsize=50, cxpb=0.5, mutpb=0.2, chunk_size=8, checkpoint_interval=2, verbose = True, plot_every_best = False):
    
    checkpoint_file = "checkpoint_ngen{}_popsize{}_initPorb{}_initMp{}.pkl".format(ngen, popsize, mainconfig["params_init"]["porb"], mainconfig["params_init"]["M1"])

    hof = tools.ParetoFront()
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.median, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)
    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])
    

    df = pd.DataFrame([mainconfig["params_init"]])
    df.to_csv(f'{mainconfig["saveto"]}/init_state_initPorb{mainconfig["params_init"]["porb"]}_initMp{mainconfig["params_init"]["M1"]}_ngen{ngen}_popsize{popsize}.dat', sep="\t", index = False)

    if os.path.isfile(checkpoint_file):
        with open(checkpoint_file, "rb") as cp_file:
            cp = pickle.load(cp_file)
            pop, start_gen, hof, logbook, rndstate = cp["population"], cp["generation"], cp["halloffame"], cp["logbook"], cp["rndstate"]
            random.setstate(rndstate)
    else:
        start_gen = 0
        pop = toolbox.population(n=popsize)
        chunks = [pop[i:i + chunk_size] for i in range(0, len(pop), chunk_size)]
        fitness_results = []
        for chunk in chunks:
            fitness_results.extend(evaluate_chunk(chunk))
        for ind, fit in zip(pop, fitness_results):
            ind.fitness.values = fit


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
            df.to_csv(f'{mainconfig["saveto"]}/gen{gen}_initPorb{mainconfig["params_init"]["porb"]}_initMp{mainconfig["params_init"]["M1"]}_ngen{ngen}_popsize{popsize}.dat', sep="\t")

        
        offspring = algorithms.varAnd(pop, toolbox, cxpb, mutpb)
        fits = toolbox.map(toolbox.evaluate, offspring)
        for fit, ind in zip(fits, offspring):
            ind.fitness.values = fit

        pop[:] = tools.selNSGA2(offspring + pop, popsize)

        hof.update(pop)

        record = stats.compile(pop)
        logbook.record(gen=gen, nevals=len(offspring), **record)
        
        weights = calculate_weights_from_pareto_front(pop)
        logging.info(f"Calculated Weights: {weights}")
        
        best_solution, weighted_sums = weighted_sum_method(pop, weights)
        logging.info(f"Objectives: {best_solution.fitness.values}")

        best_params = dict(zip([param for param in mainconfig["params_opt"] if param in mainconfig["params_bounds"]], best_solution))
        full_best_params = mainconfig["params_init"].copy()
        full_best_params.update(best_params)

        logging.info("Best parameters:")
        for key, value in full_best_params.items():
            logging.info(f"{key}: {value}")
        
        if plot_every_best:
            transformed_params = {
                "q" : mainconfig["params_init"]['m2'], 
                "a" : mainconfig["params_init"]['asr'] * mainconfig["params_init"]['R1'] * (1 + mainconfig["params_init"]['r2']), 
                "sbratio" : 1.0,
                }
            
            orb = orbit(porb = mainconfig["params_init"]["porb"], inc_deg = mainconfig["params_init"]["incl"], t0 = mainconfig["params_init"]["t0"], a = transformed_params["a"], q = transformed_params["q"], f_c = mainconfig["params_init"]["f_c"], f_s = mainconfig["params_init"]["f_s"])
            
            system_obs = binary(orbit = orb, primary = star(), secondary = star(), sbratio = transformed_params["sbratio"], l3 = mainconfig["params_init"]["l3"], gamma = mainconfig["params_init"]["gamma"],
                            obsLC = lcs, obsSpec = specs, nbins = mainconfig["nbins"],  gdc_option = "opt", phase1 = mainconfig["phase1"], phase2 = mainconfig["phase2"], oversampling=mainconfig["oversampling"], width=mainconfig["width"]) 
            
            system_obs.UPDATE_PARAMS(mainconfig["params_init"])
            
            process_id = current_process()._identity[0] if current_process()._identity else 1
            
            synthV_instance = synthini_s[process_id - 1]
            synthV_imu_instance = synthini_simu[process_id - 1]
            
            system_obs = update_system(system_obs, full_best_params, synthV_instance, synthV_imu_instance, ifplot_main=True)

        logging.info(f'GENERATION {gen} FINISHED')

    return pop, logbook, hof


def run_ga(ngen=20, popsize=50, cxpb=0.5, mutpb=0.2, chunk_size=8):
    pop = toolbox.population(n=popsize)
    hof = tools.ParetoFront()
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.median, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)

    chunks = [pop[i:i + chunk_size] for i in range(0, len(pop), chunk_size)]
    fitness_results = []

    for chunk in chunks:
        fitness_results.extend(evaluate_chunk(chunk))

    algorithms.eaMuPlusLambda(pop, toolbox, mu=popsize, lambda_=popsize, cxpb=cxpb, mutpb=mutpb,
                              ngen=ngen, stats=stats, halloffame=hof, verbose=True)

    return pop, stats, hof





def main():
    
    
    transformed_params = {
        "q" : mainconfig["params_init"]['m2'], 
        "a" : mainconfig["params_init"]['asr'] * mainconfig["params_init"]['R1'] * (1 + mainconfig["params_init"]['r2']), 
        "sbratio" : 1.0,
        }
    
    orb = orbit(porb = mainconfig["params_init"]["porb"], inc_deg = mainconfig["params_init"]["incl"], t0 = mainconfig["params_init"]["t0"], a = transformed_params["a"], q = transformed_params["q"], f_c = mainconfig["params_init"]["f_c"], f_s = mainconfig["params_init"]["f_s"])
    
    system_obs = binary(orbit = orb, primary = star(), secondary = star(), sbratio = transformed_params["sbratio"], l3 = mainconfig["params_init"]["l3"], gamma = mainconfig["params_init"]["gamma"],
                    obsLC = lcs, obsSpec = specs, nbins = mainconfig["nbins"],  gdc_option = "opt", phase1 = mainconfig["phase1"], phase2 = mainconfig["phase2"], oversampling=mainconfig["oversampling"], width=mainconfig["width"]) 
    
    system_obs.UPDATE_PARAMS(mainconfig["params_init"])
    
    global pool 
    pool = Pool(processes=int(mainconfig["nthreads"]))
    
    process_id = current_process()._identity[0] if current_process()._identity else 1
    
    synthV_instance = synthini_s[process_id - 1]
    synthV_imu_instance = synthini_simu[process_id - 1]
    
    system_obs = update_system(system_obs, mainconfig["params_init"], synthV_instance, synthV_imu_instance, ifplot_main=True, if_save_mods = False)
    
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
    toolbox.register("mate", tools.cxBlend, alpha=0.5)
    toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=1, indpb=0.1)
    toolbox.register("map", custom_map)
    toolbox.register("select", tools.selNSGA2)
    
    
    pop, stats, hof = run_ga_checkpoints(
        ngen = mainconfig["ngen"], 
        popsize = mainconfig["popsize"] , 
        cxpb= mainconfig["cxpb"], 
        mutpb= mainconfig["mutpb"],
        chunk_size= mainconfig["nthreads"],
        checkpoint_interval=mainconfig["checkpoint_interval"],
        verbose = False,
        plot_every_best=True,
        )
    
    pool.close()
    pool.join()

    weights = calculate_weights_from_pareto_front(pop)
    
    best_solution, weighted_sums = weighted_sum_method(pop, weights)
    
    evaluate(best_solution)

    best_params = dict(zip([param for param in mainconfig["params_opt"] if param in mainconfig["params_bounds"]], best_solution))
    full_best_params = mainconfig["params_init"].copy()
    full_best_params.update(best_params)

    logging.critical("Best parameters:")
    for key, value in full_best_params.items():
        logging.critical(f"{key}: {value}")
        
    system_obs = update_system(system_obs, full_best_params, synthV_instance, synthV_imu_instance, ifplot_main=True)

    

if __name__ == '__main__':
    main()
    
    