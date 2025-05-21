import sys, json, os
import argparse

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
from scipy.optimize import approx_fprime
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.interpolate import LinearNDInterpolator

from PIL import Image
import seaborn as sns
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib.transforms as transforms
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from matplotlib.patches import Ellipse
from PyQt5 import QtWidgets, uic, QtGui, QtCore
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar
import mplcursors
from matplotlib.colors import Normalize

from functools import partial

import repair_local as repair
from lcrvcurve import lightcurve
from spectra import specSeries
from spectra import synthetic
from orbit import orbit
from star import star
from binary import binary

import astropy.units as au
import astropy.constants as ac

from deap import base, creator, tools

import multiprocessing as mp


sun_teff = 5777
sol_sp = 5777**4 / 10**(4.43)
G = 6.67430e-8


msun_mesa = 1.9884098706980504E+033
rsun_mesa = 6.9570000000000000E+010
lsun_mesa = 3.8280000000000003E+033


def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])

    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    sigma_x = np.sqrt(cov[0, 0])
    sigma_y = np.sqrt(cov[1, 1])

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse), pearson, sigma_x, sigma_y


def CALCULATE_SEMI_AMPLITUDE(a_i, porb, incl, ecc):
    a_i *= au.R_sun
    a_i = a_i.to('km').value
    period = porb * au.day
    period = period.to('second').value
    return 2.*np.pi*a_i*np.sin(np.deg2rad(incl)) / np.sqrt(1.-ecc**2) / period

def CALCULATE_MASSES(k1,k2,porb,incl,ecc):
    k1 = k1 * au.km / au.second
    k2 = k2 * au.km / au.second
    period = porb*au.day
    scale = 0.5/(np.pi*ac.G) * (1.-ecc**2)**(1.5) * (k1+k2)**2 * period

    m1 = ((scale * k2).to('Msun')).value
    m2 = ((scale * k1).to('Msun')).value
    sin3i = np.sin(np.deg2rad(incl)) **3
    return m1/sin3i ,m2/sin3i

def SOLVE_ECC_OMEGA(fc,fs):
    #fs = sqrt(e)*sin(omega)
    #fc = sqrt(e)*cos(omega)

    omega = np.arctan2(fs,fc)
    ecc = fs**2 + fc**2

    return ecc,omega

def deap_population_to_df(population, labels, if_add_columns = False, weights = [1,1,1,1]):

    columns = labels + ['M1','gen','obj1', 'obj2', 'obj3', 'obj4']

    data = []

    for individual in population:
        row = individual[:]  # Get the individual's data
        row.extend(individual.fitness.values)  # Add the fitness values
        data.append(row)

    df = pd.DataFrame(data, columns=columns)

    if if_add_columns:
        df['sum_obj'] =  df['obj1']+df['obj2']+df['obj3']+df['obj4']
        df['wsum_obj'] =  df['obj1']*weights[0]+df['obj2']*weights[1]+df['obj3']*weights[2]+df['obj4']*weights[3]
        df['rating'] = range(len(df) + 1, 1, -1)
        df['Unnamed: 0'] = range(len(df) + 1, 1, -1)

    return df

def df_to_deap_population(df, labels, sort_by_fitness=False):


    labels = labels + ['M1','gen'] #'obj1', 'obj2', 'obj3', 'obj4']

    population = []
    for _, row in df.iterrows():
        #print(row)
        #print(labels)
        #print(row[labels].tolist())
        individual = creator.Individual(row[labels].tolist())
        individual.fitness.values = (row['obj1'], row['obj2'], row['obj3'], row['obj4'])
        population.append(individual)

    if sort_by_fitness:
        population.sort(key=lambda ind: ind.fitness.values)


    return population


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


def load_isochrones(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Find the index where the data starts
    header_end_index = None
    for i, line in enumerate(lines):
        if line.startswith("# EEP"):
            header_end_index = i + 1  # Data starts after this line
            break

    if header_end_index is None:
        raise ValueError("Could not find the end of the header.")

    data = pd.read_csv(file_path, skiprows=header_end_index, delim_whitespace=True, comment='#')

    # Adding column names
    col_names = ["EEP", "log10_isochrone_age_yr", "initial_mass", "star_mass", "star_mdot", "he_core_mass",
                 "c_core_mass", "o_core_mass", "log_L", "log_L_div_Ledd", "log_LH", "log_LHe", "log_LZ",
                 "log_Teff", "log_abs_Lgrav", "log_R", "log_g", "log_surf_z", "surf_avg_omega",
                 "surf_avg_v_rot", "surf_num_c12_div_num_o16", "v_wind_Km_per_s", "surf_avg_omega_crit",
                 "surf_avg_omega_div_omega_crit", "surf_avg_v_crit", "surf_avg_v_div_v_crit",
                 "surf_avg_Lrad_div_Ledd", "v_div_csound_surf", "surface_h1", "surface_he3",
                 "surface_he4", "surface_li7", "surface_be9", "surface_b11", "surface_c12", "surface_c13",
                 "surface_n14", "surface_o16", "surface_f19", "surface_ne20", "surface_na23",
                 "surface_mg24", "surface_si28", "surface_s32", "surface_ca40", "surface_ti48",
                 "surface_fe56", "log_center_T", "log_center_Rho", "center_degeneracy", "center_omega",
                 "center_gamma", "mass_conv_core", "center_h1", "center_he4", "center_c12", "center_n14",
                 "center_o16", "center_ne20", "center_mg24", "center_si28", "pp", "cno", "tri_alfa",
                 "burn_c", "burn_n", "burn_o", "c12_c12", "delta_nu", "delta_Pg", "nu_max",
                 "acoustic_cutoff", "max_conv_vel_div_csound", "max_gradT_div_grada", "gradT_excess_alpha",
                 "min_Pgas_div_P", "max_L_rad_div_Ledd", "e_thermal", "phase"]

    data.columns = col_names[:data.shape[1]]  # Trim column names to the correct number
    data["log_L_spec"] = np.log10 ( ( (10**data["log_Teff"].astype(float))**4.0 / 10.0**data["log_g"].astype(float)) / sol_sp )
    return data

def load_tracks(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    header_end_index = None
    for i, line in enumerate(lines):
        if line.startswith("# EEPs:"):
            header_end_index = i + 1  # Data starts after this line
            break

    if header_end_index is None:
        raise ValueError("Could not find the end of the header.")

    data = pd.read_csv(file_path, skiprows=header_end_index, delim_whitespace=True, comment='#')

    # Adding column names
    col_names = ["star_age", "star_mass", "star_mdot", "he_core_mass", "c_core_mass", "o_core_mass",
                 "log_L", "log_L_div_Ledd", "log_LH", "log_LHe", "log_LZ", "log_Teff", "log_abs_Lgrav",
                 "log_R", "log_g", "log_surf_z", "surf_avg_omega", "surf_avg_v_rot", "surf_num_c12_div_num_o16",
                 "v_wind_Km_per_s", "surf_avg_omega_crit", "surf_avg_omega_div_omega_crit", "surf_avg_v_crit",
                 "surf_avg_v_div_v_crit", "surf_avg_Lrad_div_Ledd", "v_div_csound_surf", "surface_h1",
                 "surface_he3", "surface_he4", "surface_li7", "surface_be9", "surface_b11", "surface_c12",
                 "surface_c13", "surface_n14", "surface_o16", "surface_f19", "surface_ne20", "surface_na23",
                 "surface_mg24", "surface_si28", "surface_s32", "surface_ca40", "surface_ti48", "surface_fe56",
                 "log_center_T", "log_center_Rho", "center_degeneracy", "center_omega", "center_gamma",
                 "mass_conv_core", "center_h1", "center_he4", "center_c12", "center_n14", "center_o16",
                 "center_ne20", "center_mg24", "center_si28", "pp", "cno", "tri_alfa", "burn_c", "burn_n",
                 "burn_o", "c12_c12", "delta_nu", "delta_Pg", "nu_max", "acoustic_cutoff",
                 "max_conv_vel_div_csound", "max_gradT_div_grada", "gradT_excess_alpha", "min_Pgas_div_P",
                 "max_L_rad_div_Ledd", "e_thermal", "phase"]

    data.columns = col_names[:data.shape[1]]  # Trim column names to the correct number
    data["log_L_spec"] = np.log10 ( ( (10**data["log_Teff"].astype(float))**4.0 / 10.0**data["log_g"].astype(float)) / sol_sp )
    return data

class MyApp(QtWidgets.QMainWindow):
    def __init__(self, filepath, stepgen, parent=None):


        self.filepath = filepath
        self.stepgen = stepgen
        self.passband = None
        self.load_config()

        print(len(self.dfs_gens))
        
        # mist isochrones and tracks for [Fe/H] closest to mainconfig['params_init']['metallicity1']
        grid_nodes = np.array([-4.00, -3.50, -3.00, -2.50, -2.00, -1.75, -1.50, -1.25, 
                               -1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50])
        closest_node = grid_nodes[np.argmin(np.abs(grid_nodes - self.mainconfig['params_init']['metallicity1']))]
        formatted_str = f"p{closest_node:.2f}" if closest_node > -0.001 else f"m{abs(closest_node):.2f}"
        
        isochrones_file = f'/Users/nadezhda/Dropbox/MIST/full/MIST_v1.2_vvcrit0.0_full_isos/MIST_v1.2_feh_{formatted_str}_afe_p0.0_vvcrit0.0_full.iso'
        path = f'/Users/nadezhda/Dropbox/MIST/full/MIST_v1.2_feh_{formatted_str}_afe_p0.0_vvcrit0.0_EEPS/'
        tracks_files = [f for f in os.listdir(path) if f.endswith('.eep')]
        tracks_files = sorted(tracks_files)
        print(path)
    
        self.isochrones = load_isochrones(isochrones_file)
        
        self.tracks = {}
        mass_prev = 0
        for track_file in tracks_files:
            mass = float(track_file.split('M.track.eep')[0]) / 100
            
            if mass-mass_prev < 0.1999:
                continue
            mass_prev = mass
                
            self.tracks[mass] = load_tracks(path+track_file)
            

        super(MyApp, self).__init__(parent)
        uic.loadUi('gui.ui', self)

        plt.rcParams['figure.facecolor'] = '#ececec'
        plt.rcParams['axes.facecolor'] = '#ececec'

        # tab1
        self.mainPlotWidget = self.findChild(QtWidgets.QWidget, 'wMain')
        self.scrollArea = self.findChild(QtWidgets.QScrollArea, 'scrlMain')
        self.scrollAreaWidgetContents = self.findChild(QtWidgets.QWidget, 'scrlThumb')
        self.zoomInButton = self.findChild(QtWidgets.QPushButton, 'btnZoomIn')
        self.zoomOutButton = self.findChild(QtWidgets.QPushButton, 'btnZoomOut')


        # Create a QLabel for displaying the image
        self.imageLabel = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.scrollArea.setWidget(self.imageLabel)

        # Initialize buffers for figures and thumbnails
        self.figure_files = []
        self.thumbnail_files = []

        # Initialize zoom factor
        self.zoom_factor = 0.5
        self.current_index = 0

        # Generate figures and thumbnails
        self.generate_figures_and_thumbnails()

        # Create thumbnail panel
        self.create_thumbnail_panel()

        # Connect buttons to their respective slots
        self.zoomInButton.clicked.connect(self.zoom_in)
        self.zoomOutButton.clicked.connect(self.zoom_out)


        self.imageLabel.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.imageLabel.customContextMenuRequested.connect(self.show_context_menu)


        # tab 2

        self.tabSelector = self.findChild(QtWidgets.QWidget, 'tabSelector')
        self.scrollableFrameSelector = self.findChild(QtWidgets.QWidget, 'scrlSelector')
        self.scrollableFrameSelected = self.findChild(QtWidgets.QWidget, 'scrlSelected')

        self.tab2_layout_created = False
        self.tabWidget.currentChanged.connect(self.on_tab_change)

        self.btnFilterLCChi2 = self.findChild(QtWidgets.QPushButton, 'btnFilterLCChi2')
        self.txtFilterLCChi2 = self.findChild(QtWidgets.QLineEdit, 'txtFilterLCChi2')
        self.btnFilterDisChi2 = self.findChild(QtWidgets.QPushButton, 'btnFilterDisChi2')
        self.txtFilterDisChi2 = self.findChild(QtWidgets.QLineEdit, 'txtFilterDisChi2')
        self.btnFilterSp1Chi2 = self.findChild(QtWidgets.QPushButton, 'btnFilterSp1Chi2')
        self.txtFilterSp1Chi2 = self.findChild(QtWidgets.QLineEdit, 'txtFilterSp1Chi2')
        self.btnFilterSp2Chi2 = self.findChild(QtWidgets.QPushButton, 'btnFilterSp2Chi2')
        self.txtFilterSp2Chi2 = self.findChild(QtWidgets.QLineEdit, 'txtFilterSp2Chi2')
        self.btnFilterNgen = self.findChild(QtWidgets.QPushButton, 'btnFilterNgen')
        self.txtFilterNgen = self.findChild(QtWidgets.QLineEdit, 'txtFilterNgen')
        self.btnFilterLong = self.findChild(QtWidgets.QPushButton, 'btnFilterLong')
        self.txtFilterLong = self.findChild(QtWidgets.QLineEdit, 'txtFilterLong')

        self.rbByGen = self.findChild(QtWidgets.QRadioButton, 'rbByGen')
        self.rbByRating = self.findChild(QtWidgets.QRadioButton, 'rbByRating')

        self.radioGroup2 = QtWidgets.QButtonGroup(self)
        self.radioGroup2.addButton(self.rbByGen)
        self.radioGroup2.addButton(self.rbByRating)
        self.rbByGen.setChecked(True)

        self.btnFilterZoom = self.findChild(QtWidgets.QPushButton, 'btnFilterZoom')

        self.btnFilterReset = self.findChild(QtWidgets.QPushButton, 'btnFilterReset')

        self.btnFilterLCChi2.clicked.connect(lambda: self.update_tab2_filter('LCChi2'))
        self.btnFilterDisChi2.clicked.connect(lambda: self.update_tab2_filter('DisChi2'))
        self.btnFilterSp1Chi2.clicked.connect(lambda: self.update_tab2_filter('Sp1Chi2'))
        self.btnFilterSp2Chi2.clicked.connect(lambda: self.update_tab2_filter('Sp2Chi2'))
        self.btnFilterNgen.clicked.connect(lambda: self.update_tab2_filter('Ngen'))
        self.btnFilterLong.clicked.connect(lambda: self.update_tab2_filter('Long'))
        self.btnFilterReset.clicked.connect(lambda: self.update_tab2_filter('Reset'))
        self.btnFilterZoom.clicked.connect(lambda: self.update_tab2_filter('Zoom'))


        self.rbConstrNone = self.findChild(QtWidgets.QRadioButton, 'rbConstrNone')
        self.rbConstrMass = self.findChild(QtWidgets.QRadioButton, 'rbConstrMass')
        self.rbConstrMassSum = self.findChild(QtWidgets.QRadioButton, 'rbConstrMassSum')
        self.rbConstrAge = self.findChild(QtWidgets.QRadioButton, 'rbConstrAge')

        self.radioGroup3 = QtWidgets.QButtonGroup(self)
        self.radioGroup3.addButton(self.rbConstrNone)
        self.radioGroup3.addButton(self.rbConstrMass)
        self.radioGroup3.addButton(self.rbConstrMassSum)
        self.radioGroup3.addButton(self.rbConstrAge)
        self.rbConstrNone.setChecked(True)


        # tab 3

        self.tabDetailed = self.findChild(QtWidgets.QWidget, 'tabDetailed')
        self.scrlConfig = self.findChild(QtWidgets.QScrollArea, 'scrlConfig')
        self.txtConfig = self.scrlConfig.findChild(QtWidgets.QPlainTextEdit, 'txtConfig')
        self.scrlParams = self.findChild(QtWidgets.QScrollArea, 'scrlParams')
        self.txtParams = self.scrlParams.findChild(QtWidgets.QPlainTextEdit, 'txtParams')
        self.scrlTransfParams = self.findChild(QtWidgets.QScrollArea, 'scrlTransfParams')
        self.txtTransfParams = self.scrlTransfParams.findChild(QtWidgets.QPlainTextEdit, 'txtTransfParams')
        self.scrlPlottedSolution = self.findChild(QtWidgets.QScrollArea, 'scrlPlottedSolution')
        self.scrlHRD = self.findChild(QtWidgets.QScrollArea, 'scrlHRD')
        self.computeButton = self.findChild(QtWidgets.QPushButton, 'btnRecompute')
        self.saveButton = self.findChild(QtWidgets.QPushButton, 'btnSaveAll')

        self.tab3_layout_created = False

        self.computeButton.clicked.connect(self.update_tab3_content)
        self.saveButton.clicked.connect(self.save_solution)


        self.rbPhotLum = self.findChild(QtWidgets.QRadioButton, 'rbPhotLum')
        self.rbSpecLum = self.findChild(QtWidgets.QRadioButton, 'rbSpecLum')
        self.radioGroup = QtWidgets.QButtonGroup(self)
        self.radioGroup.addButton(self.rbPhotLum)
        self.radioGroup.addButton(self.rbSpecLum)
        self.rbSpecLum.setChecked(True)



        # tab 4

        self.tabOpt = self.findChild(QtWidgets.QWidget, 'tabOpt')
        self.scrlConfig_opt = self.findChild(QtWidgets.QScrollArea, 'scrlConfig_opt')
        self.txtConfig_opt = self.scrlConfig_opt.findChild(QtWidgets.QPlainTextEdit, 'txtConfig_opt')
        self.scrlParams_opt = self.findChild(QtWidgets.QScrollArea, 'scrlParams_opt')
        self.txtParams_opt = self.scrlParams_opt.findChild(QtWidgets.QPlainTextEdit, 'txtParams_opt')
        self.scrlTransfParams_opt = self.findChild(QtWidgets.QScrollArea, 'scrlTransfParams_opt')
        self.txtTransfParams_opt = self.scrlTransfParams_opt.findChild(QtWidgets.QPlainTextEdit, 'txtTransfParams_opt')
        self.scrlPlottedSolution_opt = self.findChild(QtWidgets.QScrollArea, 'scrlPlottedSolution_opt')
        self.scrlHRD_opt = self.findChild(QtWidgets.QScrollArea, 'scrlHRD_opt')
        self.btnRun_opt = self.findChild(QtWidgets.QPushButton, 'btnRun_opt')

        self.tab4_layout_created = False

        self.btnRun_opt.clicked.connect(self.run_opt_button)


        self.rbPhotLum_opt = self.findChild(QtWidgets.QRadioButton, 'rbPhotLum_opt')
        self.rbSpecLum_opt = self.findChild(QtWidgets.QRadioButton, 'rbSpecLum_opt')
        self.radioGroup_opt = QtWidgets.QButtonGroup(self)
        self.radioGroup_opt.addButton(self.rbPhotLum_opt)
        self.radioGroup_opt.addButton(self.rbSpecLum_opt)
        self.rbSpecLum_opt.setChecked(True)


    def on_tab_change(self, index):

        if self.tabWidget.widget(index) == self.tabSelector and not self.tab2_layout_created:
            self.create_tab2_layout()
            self.tab2_layout_created = True
        if self.tabWidget.widget(index) == self.tabDetailed and not self.tab3_layout_created:
            self.create_tab3_layout()
            self.tab3_layout_created = True
        if self.tabWidget.widget(index) == self.tabOpt and not self.tab4_layout_created:
            self.create_tab4_layout()
            self.tab4_layout_created = True


    def create_tab3_layout(self):

        plt.rcParams['figure.facecolor'] = '#ececec'
        plt.rcParams['axes.facecolor'] = '#ececec'

        config_output = {
                'selected_row': self.selected_row.to_dict(),
                'params_init': self.mainconfig['params_init']
            }
        self.txtConfig.setPlainText(json.dumps(config_output, indent=4))

        if hasattr(self, 'selected_row') and self.selected_row is not None:
            params = self.mainconfig['params_init'].copy()
            for param in self.selected_row.index:
                if param in params:
                    params[param] = self.selected_row[param]


        self.update_tab3_content()

    def create_tab4_layout(self):

        #print(self.weights)

        plt.rcParams['figure.facecolor'] = '#ececec'
        plt.rcParams['axes.facecolor'] = '#ececec'

        #config_output = {
        #        'optimization': {"obj_weights":[self.weights[0],self.weights[1],self.weights[2],self.weights[3]], "labels":['r1','r2'], "bounds":[[0.2,0.3],[0.01,0.1]],
        #                         "method":"LM", "tol":1e-8, "options":{'maxfun': 500, 'disp': True} },
        #        'selected_row': self.selected_row.to_dict(),
        #        'params_init': self.mainconfig['params_init']
        #    }
        config_output = {
                'optimization': {"labels":['r1','r2','q','incl','a'], "obj_weights":[self.weights[0],self.weights[1],self.weights[2],self.weights[3]],
                                 "method":"LM", "tol":1e-8,'epsfcn':1e-6, "options":{'maxfun': 500, 'disp': True} },
                'selected_row': self.selected_row.to_dict(),
                'params_init': self.mainconfig['params_init']
            }
        self.txtConfig_opt.setPlainText(json.dumps(config_output, indent=4))

        if hasattr(self, 'selected_row') and self.selected_row is not None:
            params = self.mainconfig['params_init'].copy()
            for param in self.selected_row.index:
                if param in params:
                    params[param] = self.selected_row[param]


        if (not hasattr(self, 'isochrones')) or self.isochrones is None:
            isochrones_file = '/Users/nadezhda/Dropbox/MIST/full/MIST_v1.2_vvcrit0.0_full_isos/MIST_v1.2_feh_p0.00_afe_p0.0_vvcrit0.0_full.iso'
            path = '/Users/nadezhda/Dropbox/MIST/full/MIST_v1.2_feh_p0.00_afe_p0.0_vvcrit0.0_EEPS/'
            tracks_files = [f for f in os.listdir(path) if f.endswith('.eep')]
            tracks_files = sorted(tracks_files)

            self.isochrones = load_isochrones(isochrones_file)

            self.tracks = {}
            mass_prev = 0
            for track_file in tracks_files:
                mass = float(track_file.split('M.track.eep')[0]) / 100

                if mass-mass_prev < 0.1999:
                    continue
                mass_prev = mass

                self.tracks[mass] = load_tracks(path+track_file)


        self.update_tab4_content()



    def update_plot_in_scroll_area(self, scroll_area, fig, add_toolbar=False):

        canvas = FigureCanvas(fig)

        container = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(container)

        if add_toolbar:
            toolbar = NavigationToolbar(canvas, self)

            toolbar_container = QtWidgets.QWidget()
            toolbar_layout = QtWidgets.QVBoxLayout(toolbar_container)
            toolbar_layout.setContentsMargins(0, 0, 0, 0)
            toolbar_layout.addWidget(toolbar)

            toolbar_container.setFixedHeight(25)

            layout.addWidget(toolbar_container)

            if fig.axes:
                ax = fig.axes[0]

                self.update_labels(ax)
                ax.callbacks.connect('xlim_changed', partial(self.on_axes_limits_change, ax, canvas))
                ax.callbacks.connect('ylim_changed', partial(self.on_axes_limits_change, ax, canvas))

        canvas.draw()

        layout.addWidget(canvas)

        scroll_area.setWidget(container)


    def format_params(self, params):
        formatted_dict = "{\n"
        formatted_dict += ",\n".join([f"    '{key}': {repr(value)}" for key, value in params.items()])
        formatted_dict += "\n}"
        return formatted_dict


    def format_transf_output(self, params, unicode = True):

        formatted_label_unicode = {
            'r1': 'r₁ [relative to a]',
            'r2': 'r₂ [relative to a]',
            'q': 'q',
            'incl': 'i [deg]',
            'a': 'a [R☉]',
            'Teff1': 'Teff₁ [K]',
            'Teff2': 'Teff₂ [K]',
            'Ve1sini': 'vsini₁ [km s⁻¹]',
            'Ve2sini': 'vsini₂ [km s⁻¹]',
            'f_c': 'f_c',
            'f_s': 'f_s',
            'M1': 'M₁ [M☉]',
            'M2': 'M₂ [M☉]',
            'Ve1': 'Vₑ₁ [km s⁻¹]',
            'Ve2': 'Vₑ₂ [km s⁻¹]',
            'Vsync1': 'Vsynch₁ [km s⁻¹]',
            'Vsync2': 'Vsynch₂ [km s⁻¹]',
            'obj1': 'χ²_LC',
            'obj2': 'χ²_Sep',
            'obj3': 'χ²_Syn₁',
            'obj4': 'χ²_Syn₂',
            'ecc': 'e',
            'omega': 'Ω [deg]',
            'logg1': 'logg₁ [dex]',
            'logg2': 'logg₂ [dex]',
            'R1': 'R₁ [R☉]',
            'R2': 'R₂ [R☉]',
            'K1': 'K₁ [km s⁻¹]',
            'K2': 'K₂ [km s⁻¹]',
            'porb': 'Pₒᵣᵦ [d]',
            'gamma': 'Vᵧ [km s⁻¹]',
            'l3': 'Third light',
            't0': 'T₀ [BJD]',
            'metallicity1': '[M/H]₁ [dex]',
            'metallicity2': '[M/H]₂ [dex]',
            'log_L1_Lsun': 'Luminosity log L₁/L☉',
            'log_L2_Lsun': 'Luminosity log L₂/L☉',
            'log_specL1_Lsun': 'Spectroscopic luminosity log ℒ₁/ℒ☉',
            'log_specL2_Lsun': 'Spectroscopic luminosity log ℒ₂/ℒ☉'
        }

        if unicode:
            formatted_label = formatted_label_unicode

        strout = ""
        strout += "Computed parameters: \n"
        
        def CALCULATE_SEMI_AMPLITUDE(a_i, porb, incl, ecc):
            
            a_i *= au.R_sun
            a_i = a_i.to('km').value
            period = porb * au.day
            period = period.to('second').value
            return 2.*np.pi*a_i*np.sin(np.deg2rad(incl)) / np.sqrt(1.-ecc**2) / period

        if self.selected_system is not None:
            
            strout += formatted_label['M1'] + f':\t {self.selected_system.primary.mass:.3f}' + "\n"
            strout += formatted_label['M2'] + f':\t {self.selected_system.secondary.mass:.3f}' + "\n"
            strout += formatted_label['R1'] + f':\t {self.selected_system.primary.radius_sol:.3f}' + "\n"
            strout += formatted_label['R2'] + f':\t {self.selected_system.secondary.radius_sol:.3f}' + "\n"
            strout += formatted_label['logg1'] + f':\t {self.selected_system.primary.logg:.3f}' + "\n"
            strout += formatted_label['logg2'] + f':\t {self.selected_system.secondary.logg:.3f}' + "\n"

            strout += "\n"

            Teff1 = self.selected_system.primary.Teff
            Teff2 = self.selected_system.secondary.Teff
            logg1 = self.selected_system.primary.logg
            logg2 = self.selected_system.secondary.logg

            loglspec1 = np.log10( Teff1**4 / 10**logg1 / sol_sp)
            loglspec2 = np.log10( Teff2**4 / 10**logg2 / sol_sp)

            R1 = self.selected_system.primary.radius_sol
            R2 = self.selected_system.secondary.radius_sol

            loglphot1 = np.log10( (R1)**2 * (Teff1 / sun_teff)**4 )
            loglphot2 = np.log10( (R2)**2 * (Teff2 / sun_teff)**4 )

            strout += formatted_label['log_L1_Lsun'] + f':\t {loglphot1:.3f}' + "\n"
            strout += formatted_label['log_L2_Lsun'] + f':\t {loglphot2:.3f}' + "\n"
            strout += formatted_label['log_specL1_Lsun'] + f':\t {loglspec1:.3f}' + "\n"
            strout += formatted_label['log_specL2_Lsun'] + f':\t {loglspec2:.3f}' + "\n"

            strout += "\n"
            a1 = self.selected_system.orbit.a/(1.+(1/self.selected_system.orbit.q))
            a2 = self.selected_system.orbit.a/(1.+self.selected_system.orbit.q)
            strout += formatted_label['K1'] + f':\t {CALCULATE_SEMI_AMPLITUDE(a1,self.selected_system.orbit.porb,self.selected_system.orbit.inc_deg,self.selected_system.orbit.ecc):.2f}' + "\n"
            strout += formatted_label['K2'] + f':\t {CALCULATE_SEMI_AMPLITUDE(a2,self.selected_system.orbit.porb,self.selected_system.orbit.inc_deg,self.selected_system.orbit.ecc):.2f}' + "\n"

            strout += "\n"

            strout += formatted_label['Ve1'] + f':\t {self.selected_system.primary.v_equatorial:.3f}' + "\n"
            strout += formatted_label['Ve2'] + f':\t {self.selected_system.secondary.v_equatorial:.3f}' + "\n"
            strout += formatted_label['Vsync1'] + f':\t {self.selected_system.primary.ve_sync:.3f}' + "\n"
            strout += formatted_label['Vsync2'] + f':\t {self.selected_system.secondary.ve_sync:.3f}' + "\n"

            strout += formatted_label['ecc'] + f':\t {self.selected_system.orbit.ecc:.3f}' + "\n"
            strout += formatted_label['omega'] + f':\t {180.0/np.pi*self.selected_system.orbit.omega:.3f}' + "\n"


        strout += "\n"

        return strout

    def format_orig_output(self, params, unicode = True):

        formatted_label = {
            'r1' : r'$r_1$ [relative to $a$]',
            'r2' : r'$r_2$ [relative to $a$]',
            'q' : r'$q$',
            'incl' : r'$i$ [deg]',
            'a' : r'$a$ [$R_{\odot}$]',
            'Teff1' : r'$T_{eff1}$ [K]',
            'Teff2' : r'$T_{eff2}$ [K]',
            'Ve1sini' : r'$vsini_1$ [$km$ $s^[-1]$]',
            'Ve2sini' : r'$vsini_2$ [$km$ $s^[-1]$]',

            'f_c' : r'$f_c$',
            'f_s' : r'$f_s$',

            'M1' : r'$M_1$ [$M_{\odot}$]',
            'M2' : r'$M_2$ [$M_{\odot}$]',

            'Ve1' : r'$V_{e1}$ [$km$ $s^[-1]$]',
            'Ve2' : r'$V_{e2}$ [$km$ $s^[-1]$]',
            'Vsync1' : r'$V_{synch1}$ [$km$ $s^[-1]$]',
            'Vsync2' : r'$V_{synch2}$ [$km$ $s^[-1]$]',

            'obj1' : r'$\chi^2_{LC}$',
            'obj2' : r'$\chi^2_{Sep}$',
            'obj3' : r'$\chi^2_{Syn1}$',
            'obj4' : r'$\chi^2_{Syn2}$',

            'ecc' : r'$e$',
            'omega' : r'$\Omega$ [deg]',

            'logg1' : r'$logg_1$ [dex]',
            'logg2' : r'$logg_2$ [dex]',
            'R1' : r'$R_1$ [$R_{\odot}$]',
            'R2' : r'$R_2$ [$R_{\odot}$]',

            'porb' : r'$P_{orb}$ [d]',
            'gamma' : r'$V_{\gamma}$ [$km$ $s^[-1]$]',
            'l3' : r'Third light',
            't0' : r'$T_0$ [BJD]',

            'metallicity1' : r'$[M/H]_1$ [dex]',
            'metallicity2' : r'$[M/H]_2$ [dex]',

            '' : r'$$',
            '' : r'$$',
            '' : r'$$',
            '' : r'$$',
            '' : r'$$',
            '' : r'$$',
            '' : r'$$',


            }

        formatted_label_unicode = {
            'r1': 'r₁ [relative to a]',
            'r2': 'r₂ [relative to a]',
            'q': 'q',
            'incl': 'i [deg]',
            'a': 'a [R☉]',
            'Teff1': 'Teff₁ [K]',
            'Teff2': 'Teff₂ [K]',
            'Ve1sini': 'vsini₁ [km s⁻¹]',
            'Ve2sini': 'vsini₂ [km s⁻¹]',
            'f_c': 'f_c',
            'f_s': 'f_s',
            'M1': 'M₁ [M☉]',
            'M2': 'M₂ [M☉]',
            'Ve1': 'Vₑ₁ [km s⁻¹]',
            'Ve2': 'Vₑ₂ [km s⁻¹]',
            'Vsync1': 'Vsynch₁ [km s⁻¹]',
            'Vsync2': 'Vsynch₂ [km s⁻¹]',
            'obj1': 'χ²_LC',
            'obj2': 'χ²_Sep',
            'obj3': 'χ²_Syn₁',
            'obj4': 'χ²_Syn₂',
            'ecc': 'e',
            'omega': 'Ω [deg]',
            'logg1': 'logg₁ [dex]',
            'logg2': 'logg₂ [dex]',
            'R1': 'R₁ [R☉]',
            'R2': 'R₂ [R☉]',
            'porb': 'Pₒᵣᵦ [d]',
            'gamma': 'Vᵧ [km s⁻¹]',
            'l3': 'Third light',
            't0': 'T₀ [BJD]',
            'metallicity1': '[M/H]₁ [dex]',
            'metallicity2': '[M/H]₂ [dex]',

            'log_L1_Lsun': 'Luminosity log L₁/L☉',
            'log_L2_Lsun': 'Luminosity log L₂/L☉',
            'log_specL1_Lsun': 'Spectroscopic luminosity log ℒ₁/ℒ☉',
            'log_specL2_Lsun': 'Spectroscopic luminosity log ℒ₂/ℒ☉'
        }

        if unicode:
            formatted_label = formatted_label_unicode

        strout = ""
        strout += f"ch2_LC: {self.selected_system.obj_lc:.5f}\n"
        strout += f"ch2_Sep: {self.selected_system.obj_dis:.5f}\n"
        strout += f"ch2_Synth1: {self.selected_system.obj_synth1:.5f}\n"
        strout += f"ch2_Synth2: {self.selected_system.obj_synth2:.5f}\n"

        strout += "\n"
        strout += "Optimized parameters: \n"

        for ip in range(len(self.mainconfig['params_opt'])):
            if self.mainconfig['params_opt'][ip] in formatted_label:
                strout += formatted_label[self.mainconfig['params_opt'][ip]] + f':\t {params[self.mainconfig["params_opt"][ip]]:.4f}' + "\n"
            else:
                strout += self.mainconfig['params_opt'][ip] + f':\t {params[self.mainconfig["params_opt"][ip]]:.4f}' + "\n"
        strout += "\n"
        strout += "Fixed parameters: \n"

        for pp in self.mainconfig['params_init']:
            if pp not in self.mainconfig['params_opt']:
                try:
                    if self.mainconfig['params_opt'][ip] in formatted_label:
                        strout += formatted_label[pp] + f':\t {params[pp]}' + "\n"
                    else:
                        if params[pp] is not None:
                            strout += pp + f':\t {params[pp]}' + "\n"
                        else:
                            strout += pp + r':\t OFF' + "\n"
                except:
                    print('')

        strout += "\n"


        return strout


    def local_opt_scipy(self, init_guess, obj_weights, bounds, labels, ind_id=1, minimum=None, method='BFGS',tol=1e-8, options={'maxfun': 500, 'disp': True}):

        def calculate_weighted_sum(fitness, weights):
            return sum((w * f) for w, f in zip(weights, fitness))

        cconfig = self.mainconfig.copy()
        cconfig['params_opt'] = labels
        cconfig['params_init'] = self.ccparams
        for jp in range(len(labels)):
            if labels[jp] in cconfig['params_bounds']:
                cconfig['params_bounds'][labels[jp]] = bounds[jp]



        self.transformed_params = {
        "sbratio": 1.0,
        }


        orb = orbit(porb=cconfig["params_init"]["porb"], inc_deg=cconfig["params_init"]["incl"], t0=cconfig["params_init"]["t0"],
                    a=cconfig["params_init"]["a"], q=cconfig["params_init"]["q"], f_c=cconfig["params_init"]["f_c"], f_s=cconfig["params_init"]["f_s"])

        system_obs = binary(orbit=orb, primary=star(), secondary=star(), sbratio=self.transformed_params["sbratio"], l3=cconfig["params_init"]["l3"], gamma=cconfig["params_init"]["gamma"],
                            obsLC=self.lcs, obsSpec=self.specs, nbins=cconfig["nbins"],  gdc_option="opt", phase1=cconfig["phase1"], phase2=cconfig["phase2"], oversampling=cconfig["oversampling"], width=cconfig["width"])

        system_obs.UPDATE_PARAMS(cconfig['params_init'])

        synthV_instance = self.synthV_core1
        synthV_imu_instance = self.synthV_imu_core1

        system_obs = repair.update_system(
            cconfig, system_obs, cconfig['params_init'], synthV_instance, synthV_imu_instance, ifplot_main=False, if_save_mods=False, passband=self.passband)

        self.selected_system = system_obs


        def objective_wrapper(x):
            print(x)
            selected_row = cconfig['params_init'].copy()
            for jp in range(len(labels)):
                if labels[jp] in selected_row:
                    selected_row[labels[jp]] = float(x[jp])
            fitness = repair.evaluate(x, mainconfig=cconfig, ind_id=ind_id)
            return calculate_weighted_sum(fitness, obj_weights)

        res = minimize(
            fun=objective_wrapper,
            x0=init_guess,
            tol=tol,
            bounds=bounds,
            method=method,
            options=options
            #options={'disp': True, }
            #options={'maxfun': max_calls, 'disp': True}  # maximum number of function evaluations
        )

        print("Optimized parameters:", res.x)
        print("Objective function value at optimum:", res.fun)


        def calculate_hessian(x, func, epsilon=1e-5):
            """Numerically estimate the Hessian matrix of `func` at point `x`."""
            n = len(x)
            hessian = np.zeros((n, n))
            fx = func(x)
            # Incremental shifts to calculate second derivatives
            for i in range(n):
                x1 = np.array(x, dtype=float)
                x1[i] += epsilon
                f1 = approx_fprime(x1, func, epsilon)
                for j in range(i, n):
                    x2 = np.array(x1, dtype=float)
                    x2[j] += epsilon
                    f2 = approx_fprime(x2, func, epsilon)
                    hessian[i, j] = (f2[j] - f1[j]) / epsilon
                    hessian[j, i] = hessian[i, j]  # Exploit symmetry
            return hessian


        if hasattr(res, 'hess_inv'):  # Check if Hessian inverse is available from the optimization
            # The Hessian inverse approximation can be used as an estimate of the covariance matrix
            hess_inv = res.hess_inv.todense() if hasattr(res.hess_inv, 'todense') else res.hess_inv
            errors = np.sqrt(np.diag(hess_inv))
            print("Approximate standard errors of the parameters:", errors)
        else:
            errors = None
            print("Hessian inverse not available for error estimation.")
            try:
                optimal_x = res.x
                hessian_matrix = calculate_hessian(optimal_x, objective_wrapper)
                hessian_inv = np.linalg.inv(hessian_matrix)
                errors = np.sqrt(np.diag(hessian_inv))
            except Exception as e:
                print(e)
                errors = None


        return res.x, res, errors

    def local_opt_curve_fit(self, init_guess, obj_weights, labels, epsfcn, ind_id=1, tol=1e-8, options={'maxfun': 500}):

        cconfig = self.mainconfig.copy()

        cconfig['params_opt'] = labels
        cconfig['params_init'] = self.ccparams

        self.transformed_params = {
            "sbratio": 1.0,
        }
        max_weight = max(obj_weights)
        threshold = max_weight / 1000.0

        orb = orbit(porb=cconfig["params_init"]["porb"], inc_deg=cconfig["params_init"]["incl"], t0=cconfig["params_init"]["t0"],
                    a=cconfig["params_init"]["a"], q=cconfig["params_init"]["q"], f_c=cconfig["params_init"]["f_c"], f_s=cconfig["params_init"]["f_s"])

        system_obs = binary(orbit=orb, primary=star(), secondary=star(), sbratio=self.transformed_params["sbratio"], l3=cconfig["params_init"]["l3"], gamma=cconfig["params_init"]["gamma"],
                            obsLC=self.lcs, obsSpec=self.specs, nbins=cconfig["nbins"],  gdc_option="opt", phase1=cconfig["phase1"], phase2=cconfig["phase2"], oversampling=cconfig["oversampling"], width=cconfig["width"])

        system_obs.UPDATE_PARAMS(cconfig['params_init'])

        synthV_instance = self.synthV_core1
        synthV_imu_instance = self.synthV_imu_core1

        system_obs = repair.update_system(
            cconfig, system_obs, cconfig['params_init'], synthV_instance, synthV_imu_instance, ifplot_main=False, if_save_mods=False, passband=self.passband)

        self.selected_system = system_obs


        def calculate_residual_for_dataset(system_obs):
            print(system_obs.primary.Teff, system_obs.secondary.Teff)
            try:
                lc_residuals = (system_obs.modLC - system_obs.obsLC.binned_flux) / system_obs.modLC * obj_weights[0] #* 200.0
            except:
                lc_residuals = np.ones_like(system_obs.obsLC.binned_flux) * obj_weights[0]

            try:
                interp_primary = interp1d(system_obs.primary.spectrum[:, 0], system_obs.primary.spectrum[:, 1], fill_value="extrapolate")
                interp_secondary = interp1d(system_obs.secondary.spectrum[:, 0], system_obs.secondary.spectrum[:, 1], fill_value="extrapolate")

                model_spec1_interp = interp_primary(system_obs.obsSpec.wl_eqlog_to_wl) 
                model_spec2_interp = interp_secondary(system_obs.obsSpec.wl_eqlog_to_wl)

                spec1_residuals_ob = (model_spec1_interp - system_obs.obsSpec.mod[0]) / model_spec1_interp * obj_weights[2]
                spec2_residuals_ob = (model_spec2_interp - system_obs.obsSpec.mod[1]) / model_spec2_interp * obj_weights[3]
            except:
                spec1_residuals_ob = np.linspace(1,1.1,1000)
                spec2_residuals_ob = np.linspace(1,1.1,1000)

            try:
                resid_spec_residuals = system_obs.residSpecNoSqNorm.flatten()[1000:-1000]* obj_weights[1]
            except:
                resid_spec_residuals = np.linspace(1,1.1,1000)

            residuals_to_stack = []
            if obj_weights[0] > threshold:
                residuals_to_stack.append(lc_residuals)
            if obj_weights[1] > threshold:
                residuals_to_stack.append(resid_spec_residuals)
            if obj_weights[2] > threshold:
                residuals_to_stack.append(spec1_residuals_ob)
            if obj_weights[3] > threshold:
                residuals_to_stack.append(spec2_residuals_ob)
            #residuals_to_stack = [lc_residuals, spec1_residuals_ob, spec2_residuals_ob, resid_spec_residuals]
            
            all_residuals = np.hstack(residuals_to_stack)

            all_residuals = np.where(np.isnan(all_residuals), 1e3, all_residuals)

            return all_residuals

        def residual_wrapper(xdata, *params, system_obs, params_init, labels):

            updated_params = params_init.copy()
            for label, value in zip(labels, params):
                updated_params[label] = value

            system_obs.UPDATE_PARAMS(updated_params)

            self.selected_row = pd.Series(updated_params)

            synthV_instance = self.synthV_core1
            synthV_imu_instance = self.synthV_imu_core1

            system_obs = repair.update_system(
                cconfig, system_obs, updated_params, synthV_instance, synthV_imu_instance, ifplot_main=False, if_save_mods=False, passband=self.passband)


            ffi = interp1d(
                system_obs.primary.spectrum[:, 0], system_obs.primary.spectrum[:, 1], fill_value="extrapolate")
            mm = ffi(system_obs.obsSpec.wl_eqlog_to_wl)
            ob1 = repair.fetch_spec_new(system_obs.obsSpec.wl_eqlog_to_wl,
                                 system_obs.obsSpec.mod[0] + 0.5, system_obs.obsSpec.wl_eqlog_to_wl, mm)


            ffi = interp1d(
                system_obs.secondary.spectrum[:, 0], system_obs.secondary.spectrum[:, 1], fill_value="extrapolate")
            mm = ffi(system_obs.obsSpec.wl_eqlog_to_wl)
            ob2 = repair.fetch_spec_new(system_obs.obsSpec.wl_eqlog_to_wl,
                                 system_obs.obsSpec.mod[1] + 0.5, system_obs.obsSpec.wl_eqlog_to_wl, mm)

            system_obs.obsSpec.mod[0] = ob1
            system_obs.obsSpec.mod[1] = ob2

            residuals = calculate_residual_for_dataset(system_obs)

            if xdata is not None:
                if len(residuals) > len(xdata):
                    residuals = residuals[:len(xdata)]
                elif len(residuals) < len(xdata):
                    residuals_buff = residuals
                    residuals = np.ones_like(xdata)
                    for hh in range(len(residuals_buff)):
                        residuals[hh] = residuals_buff[hh]

            return residuals

        initial_params = [cconfig['params_init'][label] for label in labels]

        self.selected_system = system_obs

        try:

            initial_residuals = residual_wrapper(
                None,
                *initial_params,
                system_obs=system_obs,
                params_init=cconfig['params_init'],
                labels=labels
                )


            # Create dummy xdata and ydata with the same size as the residuals
            xdata_dummy = np.zeros_like(initial_residuals)
            ydata_dummy = np.zeros_like(initial_residuals)

            optimized_params, pcov, infodict, mesg, ier  = curve_fit(
                f= lambda xdata, *params: residual_wrapper(xdata, *params, system_obs=system_obs, params_init=cconfig['params_init'], labels=labels),
                xdata=xdata_dummy,  # Pass the dummy xdata
                ydata=ydata_dummy,  # ydata is zero since we're minimizing residuals
                p0=initial_params,  # Initial guess for parameters
                xtol=tol,
                maxfev=options.get('maxfun', 500),
                absolute_sigma=True,
                full_output=True,
                epsfcn = epsfcn
                )

        except Exception as e:
            print(e)
            return None,None,None

        print(mesg)
        print(ier)
        print("Optimized parameters:", optimized_params)

        # Calculate standard errors from the covariance matrix
        if pcov is not None:
            errors = np.sqrt(np.diag(pcov))
            print("Approximate standard errors of the parameters:", errors)
        else:
            errors = None
            print("Covariance matrix not available for error estimation.")

        return optimized_params, infodict, errors


    def run_optimization_lm(self, row, obj_weights, labels, epsfcn, tol, options={'maxfun': 500}):
        init_guess = [row[label] for label in labels]
        return self.local_opt_curve_fit(init_guess, obj_weights, labels, epsfcn, ind_id=1, tol=tol, options=options)

    def run_optimization(self, row, obj_weights, labels, bounds, method, tol, options={'maxfun': 500, 'disp': True}):
        init_guess = [row[label] for label in labels]
        return self.local_opt_scipy(init_guess, obj_weights, bounds, labels, ind_id=1, method=method, tol=tol, options=options)



    def update_tab4_content(self, opt=False):

        plt.rcParams['figure.facecolor'] = '#ececec'
        plt.rcParams['axes.facecolor'] = '#ececec'

        try:
            config_text_opt = self.txtConfig_opt.toPlainText()
            config_data_opt = json.loads(config_text_opt)

            if 'params_init' in config_data_opt:
                self.mainconfig['params_init'] = config_data_opt['params_init']

            if 'selected_row' in config_data_opt:
                self.selected_row = pd.Series(config_data_opt['selected_row'])
                try:
                    self.selected_row.drop(columns=['rating'])
                except:
                    print('')
                try:
                    self.selected_row.drop(columns=['Unnamed: 0'])
                except:
                    print('')
                try:
                    self.selected_row = self.selected_row.drop(columns=['gen','obj1',
                                                  'obj2', 'obj3', 'obj4', 'sum_obj', 'wsum_obj', 'M1' ])
                except:
                    print('nothing dropped')
            if 'optimization' in config_data_opt:
                self.optimization = pd.Series(config_data_opt['optimization'])


            toadd_key = ["ldc_1","ldc_2","ld_1","ld_2","gdc_1","gdc_2","heat_1","heat_2","l3","bfac_1","bfac_2"]
            toadd_val = [None,None,"mugrid","mugrid",None,None,None,None,0.0,None,None]
            for kk in range(len(toadd_key)):
                if toadd_key[kk] not in self.selected_row:
                    self.selected_row[toadd_key[kk]] = toadd_val[kk]



        except json.JSONDecodeError as e:
            QtWidgets.QMessageBox.critical(self, "JSON Error", f"Error parsing JSON in config: {str(e)}")
            return


        params = self.mainconfig["params_init"].copy()
        toadd_key = ["ldc_1","ldc_2","ld_1","ld_2","gdc_1","gdc_2","heat_1","heat_2","l3","bfac_1","bfac_2"]
        toadd_val = [None,None,"mugrid","mugrid",None,None,None,None,0.0,None,None]
        for kk in range(len(toadd_key)):
            if toadd_key[kk] not in params:
                params[toadd_key[kk]] = toadd_val[kk]

        for param in self.selected_row.index:
            if param in params:
                params[param] = self.selected_row[param]

        if params['ldc_1'] is not None:
            params['ld_1'] = "lin"
        if params['ldc_2'] is not None:
            params['ld_2'] = "lin"



        self.ccparams = params

        config_output_opt = {
                'optimization': self.optimization.to_dict(),
                'selected_row': self.selected_row.to_dict(),
                'params_init': self.mainconfig['params_init']
            }



        self.txtConfig_opt.setPlainText(json.dumps(config_output_opt, indent=4))

        self.txtParams_opt.setPlainText(self.format_orig_output(params))

        self.txtTransfParams_opt.setPlainText(self.format_transf_output(params))

        self.update_plot_in_scroll_area(self.scrlPlottedSolution_opt, self.plot_solution(params))

        if self.rbSpecLum_opt.isChecked():
            self.update_plot_in_scroll_area(self.scrlHRD_opt, self.plot_hrd(), add_toolbar=True)
        else:
            self.update_plot_in_scroll_area(self.scrlHRD_opt, self.plot_phot_hrd(), add_toolbar=True)

        if opt:


            #if self.optimization['method'] == "LM":
            self.opt_results=self.run_optimization_lm(self.selected_row,self.optimization['obj_weights'],self.optimization['labels'],
                                                 self.optimization['epsfcn'], self.optimization['tol'],
                                                 options=self.optimization['options'])

            #else:
                #self.opt_results=self.run_optimization(self.selected_row,self.optimization['obj_weights'],self.optimization['labels'],
                #                                     self.optimization['bounds'],self.optimization['method'],self.optimization['tol'],
                #                                     options=self.optimization['options'])

            self.optimized = self.opt_results[0]


            params =  self.mainconfig['params_init'].copy()
            for param in self.selected_row.index:
                if param in params:
                    params[param] = self.selected_row[param]
            for jp in range(len(self.optimization['labels'])):
                if self.optimization['labels'][jp] in params:
                    params[self.optimization['labels'][jp]] = self.optimized[jp]


            for jp in range(len(self.optimization['labels'])):
                if self.optimization['labels'][jp] in self.selected_row:
                    self.selected_row[self.optimization['labels'][jp]] = self.optimized[jp]


            config_output_opt = {
                    'optimization': self.optimization.to_dict(),
                    'selected_row': self.selected_row.to_dict(),
                    'params_init': self.mainconfig['params_init']
                }
            self.txtConfig_opt.setPlainText(json.dumps(config_output_opt, indent=4))

            self.txtParams_opt.setPlainText(self.format_orig_output(params))

            self.txtTransfParams_opt.setPlainText(self.format_transf_output(params))

            self.update_plot_in_scroll_area(self.scrlPlottedSolution_opt, self.plot_solution(params))

            if self.rbSpecLum_opt.isChecked():
                self.update_plot_in_scroll_area(self.scrlHRD_opt, self.plot_hrd(), add_toolbar=True)
            else:
                self.update_plot_in_scroll_area(self.scrlHRD_opt, self.plot_phot_hrd(), add_toolbar=True)



    def run_opt_button(self):
        try:
            self.update_tab4_content(opt=True)
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", f"{e}")
            return

    def update_tab3_content(self):


        plt.rcParams['figure.facecolor'] = '#ececec'
        plt.rcParams['axes.facecolor'] = '#ececec'

        try:
            config_text = self.txtConfig.toPlainText()
            config_data = json.loads(config_text)

            if 'params_init' in config_data:
                self.mainconfig['params_init'] = config_data['params_init']

            if 'selected_row' in config_data:
                self.selected_row = pd.Series(config_data['selected_row'])

        except json.JSONDecodeError as e:
            QtWidgets.QMessageBox.critical(self, "JSON Error", f"Error parsing JSON in config: {str(e)}")
            return

        params = self.mainconfig['params_init'].copy()
        for param in self.selected_row.index:
            if param in params:
                params[param] = self.selected_row[param]

        config_output = {
                'selected_row': self.selected_row.to_dict(),
                'params_init': self.mainconfig['params_init']
            }
        self.txtConfig.setPlainText(json.dumps(config_output, indent=4))


        self.txtParams.setPlainText(self.format_orig_output(params))

        self.txtTransfParams.setPlainText(self.format_transf_output(params))

        self.update_plot_in_scroll_area(self.scrlPlottedSolution, self.plot_solution(params))

        if self.rbSpecLum.isChecked():
            self.update_plot_in_scroll_area(self.scrlHRD, self.plot_hrd(), add_toolbar=True)
        else:
            self.update_plot_in_scroll_area(self.scrlHRD, self.plot_phot_hrd(), add_toolbar=True)


    def save_solution(self):


        plt.rcParams['figure.facecolor'] = '#ececec'
        plt.rcParams['axes.facecolor'] = '#ececec'

        # Open a dialog for selecting a folder to save files
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Select Folder to Save Solution")

        if not folder:
            return

        config_file_path = os.path.join(folder, "config.json")
        params_file_path = os.path.join(folder, "params.txt")
        transf_params_file_path = os.path.join(folder, "transformed_params.txt")
        solution_plot_path = os.path.join(folder, "solution_plot.png")
        hrd_plot_path = os.path.join(folder, "hrd_plot.png")
        hrd_plot_path_spec = os.path.join(folder, "hrd_plot_spec.png")
        corner_ellipses_f_path = os.path.join(folder, "corner_ellipses_cutoff.png")
        corner_ellipses_last_gen_path = os.path.join(folder, "corner_ellipses_last_gen.png")
        corner_alpha_path = os.path.join(folder, "corner_alpha.png")

        lc_data_file_path = os.path.join(folder, "model_lc.txt")
        rv_curve_data_file_path = os.path.join(folder, "model_rv.txt")
        disentangled_spectrum_primary_path = os.path.join(folder, "disentangled_spectrum_primary.txt")
        disentangled_spectrum_secondary_path = os.path.join(folder, "disentangled_spectrum_secondary.txt")
        synthetic_spectrum_primary_path = os.path.join(folder, "synthetic_spectrum_primary.txt")
        synthetic_spectrum_secondary_path = os.path.join(folder, "synthetic_spectrum_secondary.txt")


        # Save config as JSON
        config_text = self.txtConfig.toPlainText()
        with open(config_file_path, 'w') as config_file:
            config_file.write(config_text)

        # Save params as plain text
        params_text = self.txtParams.toPlainText()
        with open(params_file_path, 'w') as params_file:
            params_file.write(params_text)

        # Save transformed params as plain text
        transf_params_text = self.txtTransfParams.toPlainText()
        with open(transf_params_file_path, 'w') as transf_params_file:
            transf_params_file.write(transf_params_text)

        # Save solution plot as PNG
        solution_fig = self.plot_solution(self.selected_row)
        solution_fig.savefig(solution_plot_path, format='png', dpi=300)

        # Save HRD plots as PNG

        hrd_fig_spec = self.plot_hrd()
        hrd_fig_spec.savefig(hrd_plot_path_spec, format='png', dpi=300)

        hrd_fig = self.plot_phot_hrd()
        hrd_fig.savefig(hrd_plot_path, format='png', dpi=300)

        # Save corner plots as PNG
        corner_ellipses_f = self.plot_corner_filtered(self.dfs_gens, cutoff=0.1)
        corner_ellipses_f.savefig(corner_ellipses_f_path, format='png', dpi=300)
        corner_ellipses = self.plot_corner_filtered_last(self.dfs_gens)
        corner_ellipses.savefig(corner_ellipses_last_gen_path, format='png', dpi=300)
        corner_alpha = self.plot_corner_alpha(self.dfs_gens)
        corner_alpha.savefig(corner_alpha_path, format='png', dpi=300)


        lc_data = np.column_stack((
                                    self.selected_system.obsLC.binned_phase,
                                    self.selected_system.obsLC.binned_flux,
                                    self.selected_system.modLC
                                    ))
        np.savetxt(lc_data_file_path, lc_data, header="Phase Observed_Flux Model_Flux", comments="")

        # Save RV curve data as text file (if applicable)
        rv_data = np.column_stack((
            self.selected_system.modRV1_full_phase,
            self.selected_system.modRV1_full,
            self.selected_system.modRV2_full
        ))
        np.savetxt(rv_curve_data_file_path, rv_data, header="Phase Model_RV_Primary Model_RV_Secondary", comments="")

        # Save disentangled spectra as text files
        disentangled_spectrum_primary = np.column_stack((
        self.selected_system.obsSpec.wl_eqlog_to_wl,
        self.selected_system.obsSpec.mod[0]
        ))
        disentangled_spectrum_secondary = np.column_stack((
            self.selected_system.obsSpec.wl_eqlog_to_wl,
            self.selected_system.obsSpec.mod[1]
        ))
        np.savetxt(disentangled_spectrum_primary_path, disentangled_spectrum_primary, header="Wavelength Res_Flux", comments="")
        np.savetxt(disentangled_spectrum_secondary_path, disentangled_spectrum_secondary, header="Wavelength Res_Flux", comments="")

        # Save synthetic spectra as text files
        np.savetxt(synthetic_spectrum_primary_path, self.selected_system.primary.spectrum, header="Wavelength Res_Flux", comments="")
        np.savetxt(synthetic_spectrum_secondary_path, self.selected_system.secondary.spectrum, header="Wavelength Res_Flux", comments="")


    def create_tab2_layout(self):

        plt.rcParams['figure.facecolor'] = '#ececec'
        plt.rcParams['axes.facecolor'] = '#ececec'

        extra_params = {
            'axes.titlesize': 6,
            'axes.labelsize': 6,
            'xtick.labelsize': 6,
            'ytick.labelsize': 6,
            'legend.fontsize': 6,
            'axes.linewidth': 0.5,
            'grid.linewidth': 0.5,
            'lines.linewidth': 0.5,
            'xtick.major.width': 0.5,
            'xtick.minor.width': 0.3,
            'ytick.major.width': 0.5,
            'ytick.minor.width': 0.3,
            'xtick.major.size': 3,
            'xtick.minor.size': 1.5,
            'ytick.major.size': 3,
            'ytick.minor.size': 1.5,
            'figure.figsize': (6, 6),
            'figure.dpi': 150,
            'figure.facecolor':'#ececec',
            'axes.facecolor':'#ececec'
        }

        # Create the interactive corner plot
        fig1 = self.plot_corner_alpha(self.dfs_gens, selector=True, extra_params=extra_params)
        self.fig_canvas = FigureCanvas(fig1)
        self.scrollableFrameSelectorLayout = QtWidgets.QVBoxLayout(self.scrollableFrameSelector)
        self.scrollableFrameSelectorLayout.addWidget(self.fig_canvas)

        self.selected_row = self.dfs_gens[-1].iloc[0]

        self.update_selected_plot(self.selected_row)

        aa = pd.concat(self.dfs_gens)
        self.txtFilterLCChi2.setText(f"{aa['obj1'].max():.3f}")
        self.txtFilterDisChi2.setText(f"{aa['obj2'].max():.3f}")
        self.txtFilterSp1Chi2.setText(f"{aa['obj3'].max():.3f}")
        self.txtFilterSp2Chi2.setText(f"{aa['obj4'].max():.3f}")


    def update_tab2_filter(self, filter_type):
        custom_bounds = None

        if not hasattr(self, 'filtered_dfs'):
            self.filtered_dfs = self.dfs_gens

        aa = pd.concat(self.filtered_dfs)

        try:
            

            if filter_type == 'LCChi2':
                filter_value = float(self.txtFilterLCChi2.text())
                if filter_value <= aa['obj1'].min():
                    filter_value = aa['obj1'].min() * 1.01
                    self.txtFilterLCChi2.setText(f"{filter_value:.3f}")
                filtered_data = [df[df['obj1'] < filter_value] for df in self.filtered_dfs]
                self.filtered_dfs = filtered_data

            elif filter_type == 'DisChi2':
                filter_value = float(self.txtFilterDisChi2.text())
                if filter_value <= aa['obj2'].min():
                    filter_value = aa['obj2'].min() * 1.01
                    self.txtFilterDisChi2.setText(f"{filter_value:.3f}")
                filtered_data = [df[df['obj2'] < filter_value] for df in self.filtered_dfs]
                self.filtered_dfs = filtered_data

            elif filter_type == 'Sp1Chi2':
                filter_value = float(self.txtFilterSp1Chi2.text())
                if filter_value <= aa['obj3'].min():
                    filter_value = aa['obj3'].min() * 1.01
                    self.txtFilterSp1Chi2.setText(f"{filter_value:.3f}")
                filtered_data = [df[df['obj3'] < filter_value] for df in self.filtered_dfs]
                self.filtered_dfs = filtered_data

            elif filter_type == 'Sp2Chi2':
                filter_value = float(self.txtFilterSp2Chi2.text())
                if filter_value <= aa['obj4'].min():
                    filter_value = aa['obj4'].min() * 1.01
                    self.txtFilterSp2Chi2.setText(f"{filter_value:.3f}")
                filtered_data = [df[df['obj4'] < filter_value] for df in self.filtered_dfs]
                self.filtered_dfs = filtered_data

            elif filter_type == 'Reset':
                filtered_data = self.dfs_gens
                self.filtered_dfs = filtered_data

            elif filter_type == 'Zoom':
                try:
                    filtered_data = self.filtered_dfs
                except:
                    filtered_data = self.dfs_gens
                custom_bounds = []
                aa = pd.concat(filtered_data)
                for ii in range(len(self.mainconfig['params_opt'])):
                    custom_bounds.append([])
                    bb = aa[self.mainconfig['params_opt'][ii]]
                    custom_bounds[ii] = [bb.min(), bb.max()]
                custom_bounds.append(self.custom_bounds[-1])
                print(self.custom_bounds)
                print(custom_bounds)
            else:
                pass


            self.update_tab2_plot(filtered_data, custom_bounds=custom_bounds)

        except Exception as e:
            print(e)
            QtWidgets.QMessageBox.critical(self, "Filter Error", "Please enter a valid number for the filter.")



    def update_selected_plot(self, selected_row):
        if not hasattr(self, 'scrollableFrameSelectedLayout'):
            self.scrollableFrameSelectedLayout = QtWidgets.QVBoxLayout(self.scrollableFrameSelected)
        else:
            # Clear existing widgets from layout
            for i in reversed(range(self.scrollableFrameSelectedLayout.count())):
                widget_to_remove = self.scrollableFrameSelectedLayout.itemAt(i).widget()
                self.scrollableFrameSelectedLayout.removeWidget(widget_to_remove)
                widget_to_remove.setParent(None)

        fig2 = self.plot_solution(selected_row)
        self.fig_canvas2 = FigureCanvas(fig2)
        self.scrollableFrameSelectedLayout.addWidget(self.fig_canvas2)

    def clear_layout(self, layout):
        if layout is not None:
            while layout.count():
                child = layout.takeAt(0)
                if child.widget() is not None:
                    child.widget().deleteLater()

    def update_tab2_plot(self, filtered_data=None, custom_bounds=None):
        extra_params = {
            'axes.titlesize': 6,
            'axes.labelsize': 6,
            'xtick.labelsize': 6,
            'ytick.labelsize': 6,
            'legend.fontsize': 6,
            'axes.linewidth': 0.5,
            'grid.linewidth': 0.5,
            'lines.linewidth': 0.5,
            'xtick.major.width': 0.5,
            'xtick.minor.width': 0.3,
            'ytick.major.width': 0.5,
            'ytick.minor.width': 0.3,
            'xtick.major.size': 3,
            'xtick.minor.size': 1.5,
            'ytick.major.size': 3,
            'ytick.minor.size': 1.5,
            'figure.figsize': (6, 6),
            'figure.dpi': 150,
            'figure.facecolor':'#ececec',
            'axes.facecolor':'#ececec'
        }

        data_to_plot = filtered_data if filtered_data is not None else self.dfs_gens

        self.clear_layout(self.scrollableFrameSelectorLayout)

        self.fig1 = self.plot_corner_alpha(data_to_plot, selector=True, extra_params=extra_params, custom_bounds=custom_bounds)
        self.fig_canvas = FigureCanvas(self.fig1)
        self.scrollableFrameSelectorLayout.addWidget(self.fig_canvas)
        self.fig_canvas.draw()


    def show_context_menu(self, pos):
        context_menu = QtWidgets.QMenu(self)
        save_action = context_menu.addAction("Save...")
        save_action.triggered.connect(self.save_plot)
        context_menu.exec_(self.imageLabel.mapToGlobal(pos))

    def save_plot(self):
        options = QtWidgets.QFileDialog.Options()
        file_path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save Image", "", "PNG Files (*.png);;All Files (*)", options=options)
        if file_path:
            fig_file = self.figure_files[self.current_index]
            img = Image.open(fig_file)
            img.save(file_path)

    def zoom_in(self):
        self.zoom_factor *= 1.1
        self.update_main_plot(self.current_index)

    def zoom_out(self):
        self.zoom_factor /= 1.1
        self.update_main_plot(self.current_index)

    def generate_figures_and_thumbnails(self):
        objective_type = self.mainconfig.get('objective_type', 'chi2') # chi2 = default value if key is not specified in config file

        if objective_type == 'chi2':
            plot_functions = [
                lambda: self.plot_columns_against_gen(self.dfs_gens),
                lambda: self.plot_objectives_convergence(self.dfs_gens),
                lambda: self.plot_pareto_front(self.dfs_gens),
                lambda: self.plot_corner_alpha(self.dfs_gens),
                #lambda: self.plot_corner_alpha_colour(self.dfs_gens, colour='sum_obj'),
                #lambda: self.plot_corner_alpha_colour(self.dfs_gens, colour='wsum_obj'),
                #lambda: self.plot_corner_alpha_colour(self.dfs_gens, colour='obj1'),
                lambda: self.plot_corner_filtered_last(self.dfs_gens),
                lambda: self.plot_corner_filtered(self.dfs_gens, cutoff=0.1),
            ]
        else:
            plot_functions = [
                lambda: self.plot_columns_against_gen(self.dfs_gens),
                lambda: self.plot_objectives_convergence(self.dfs_gens),
                lambda: self.plot_pareto_front(self.dfs_gens),
                lambda: self.plot_corner_alpha(self.dfs_gens),
            ]

        for i, plot_func in enumerate(plot_functions):
            fig_file, thumb_file = self.create_plot_and_thumbnail(plot_func, i)
            self.figure_files.append(fig_file)
            self.thumbnail_files.append(thumb_file)


    def create_plot_and_thumbnail(self, plot_func, index):
        
        print(index)
        
        fig = plot_func()

        # Save the figure as PNG
        fig_file = f'fig_{index}.png'
        fig.savefig(fig_file, format='png', bbox_inches='tight', dpi=150)
        plt.close(fig)

        # Create a thumbnail
        img = Image.open(fig_file)
        img.thumbnail((50, 50))
        thumb_file = f'thumb_{index}.png'
        img.save(thumb_file)

        return fig_file, thumb_file


    def create_thumbnail_panel(self):
        layout = QtWidgets.QVBoxLayout(self.scrollAreaWidgetContents)

        for i, thumb_file in enumerate(self.thumbnail_files):
            img = Image.open(thumb_file)
            qimage = self.pil_to_pixmap(img)
            pixmap = QtGui.QPixmap.fromImage(qimage)

            label = QtWidgets.QLabel(self.scrollAreaWidgetContents)
            label.setPixmap(pixmap)
            layout.addWidget(label)

            label.mousePressEvent = lambda event, idx=i: self.update_main_plot(idx)


    def pil_to_pixmap(self, im):
        """ Convert PIL Image to QPixmap """
        if im.mode == "RGB":
            r, g, b = im.split()
            im = Image.merge("RGB", (b, g, r))
        elif im.mode == "RGBA":
            r, g, b, a = im.split()
            im = Image.merge("RGBA", (b, g, r, a))
        elif im.mode == "L":
            im = im.convert("RGBA")
        im2 = im.convert("RGBA")
        data = im2.tobytes("raw", "RGBA")
        qim = QtGui.QImage(data, im2.size[0], im2.size[1], QtGui.QImage.Format_RGBA8888)
        return qim

    def update_main_plot(self, idx):
        self.current_index = idx
        fig_file = self.figure_files[idx]
        img = Image.open(fig_file)
        img = img.resize((int(img.width * self.zoom_factor), int(img.height * self.zoom_factor)))

        qimage = self.pil_to_pixmap(img)
        if qimage:
            pixmap = QtGui.QPixmap.fromImage(qimage)
            self.imageLabel.setPixmap(pixmap)
            self.imageLabel.adjustSize()

    def plot_hrd(self):

        M1 = self.selected_system.primary.mass
        M2 = self.selected_system.secondary.mass

        Teff1 = self.selected_system.primary.Teff
        Teff2 = self.selected_system.secondary.Teff

        logg1 = self.selected_system.primary.logg
        logg2 = self.selected_system.secondary.logg

        lspec1 = np.log10( Teff1**4 / 10**logg1 / sol_sp)
        lspec2 = np.log10( Teff2**4 / 10**logg2 / sol_sp)

        fig, _ = plt.subplots(figsize=(5, 5))

        main_sequence_phases = [0,1,2]

        # Filter isochrones by phase
        if 'phase' in self.isochrones.columns:
            isochrones_main_sequence = self.isochrones[self.isochrones['phase'].isin(main_sequence_phases)]
        else:
            raise ValueError("Phase column not found in isochrones data.")

        unique_ages = sorted(isochrones_main_sequence['log10_isochrone_age_yr'].unique())
        print(unique_ages)
        filtered_ages = unique_ages[::3]
        print(filtered_ages)

        # Plot filtered isochrones
        for age in filtered_ages:
            subset = isochrones_main_sequence[isochrones_main_sequence['log10_isochrone_age_yr'] == age]
            plt.plot(10**subset['log_Teff'], subset['log_L_spec'], '--k',lw=0.5,alpha=0.35, label=f'Age: {10**age:.0f} years')

        lum_min = min(lspec1, lspec2)
        lum_max = max(lspec1, lspec2)

        xlim = [min(Teff1, Teff2) - 0.1 * min(Teff1, Teff2), max(Teff1, Teff2) + 0.1 * max(Teff1, Teff2)]
        ylim = [lum_min - 0.05 * lum_min, lum_max + 0.03 * lum_max]

        if self.tracks:
            # Filter and plot tracks that touch the area between lspec1 and lspec2
            for mass, track in self.tracks.items():
                if 'phase' in track.columns:
                    track_main_sequence = track[track['phase'].isin(main_sequence_phases)]
                    if not track_main_sequence.empty:
                        track_lum_min = track_main_sequence['log_L_spec'].min()
                        track_lum_max = track_main_sequence['log_L_spec'].max()

                        # Check if the track touches the area between lspec1 and lspec2
                        if (track_lum_min <= lum_max and track_lum_max >= lum_min):
                            # Plot the track
                            plt.plot(10**track_main_sequence['log_Teff'], track_main_sequence['log_L_spec'], "-k", lw=0.5,
                                     label=f'Mass: {mass:.2f} $M_\\odot$')

                            # Update the xlim[0] only if necessary
                            start_idx = track_main_sequence.index[0]
                            start_teff = 10**(track_main_sequence.at[start_idx, 'log_Teff'])

                            xlim[1] = max(xlim[1], start_teff)

                            # Add text labels for the track start points
                            if xlim[0] * 1.02 <= start_teff <= xlim[1] * 0.98:
                                if mass < 8:
                                    label = f'{mass:.1f} $M_\\odot$'
                                else:
                                    label = f'{mass:.0f} $M_\\odot$'
                                plt.text(start_teff, track_main_sequence.at[start_idx, 'log_L_spec'] - 0.02, label, fontsize=6, ha='right', va='bottom')


        plt.xlim(xlim)
        plt.ylim(ylim)


        l1, = plt.plot(Teff1,lspec1,"*b",markersize=12, alpha=0.8,mew=0)
        l1.set_markerfacecolor((0.404, 0.514, 1, 0.878))
        label = f'{M1:.2f} $M_\\odot$'
        plt.text(Teff1*0.9,lspec1*1.003, label, fontsize=6, ha='right', va='bottom', c="royalblue")

        l2, = plt.plot(Teff2,lspec2,"*b",markersize=12, alpha=0.8,mew=0)
        l2.set_markerfacecolor((0.404, 0.514, 1, 0.878))
        label = f'{M2:.2f} $M_\\odot$'
        plt.text(Teff2*0.9,lspec2*1.003, label, fontsize=6, ha='right', va='bottom', c="royalblue")

        plt.xlabel(r'$T_{eff}$ [K]')
        plt.ylabel(r'log $\mathcal{L}$/$\mathcal{L}_{\odot}$')
        plt.xscale('log')
        #plt.yscale('log')
        plt.gca().invert_xaxis()  # Hotter stars on the left

        return fig


    def update_labels(self, ax):

        for text in ax.texts:
            text.remove()

        x_limits = ax.get_xlim()
        y_limits = ax.get_ylim()

        M1 = self.selected_system.primary.mass
        M2 = self.selected_system.secondary.mass

        Teff1 = self.selected_system.primary.Teff
        Teff2 = self.selected_system.secondary.Teff

        logg1 = self.selected_system.primary.logg
        logg2 = self.selected_system.secondary.logg

        R1 = self.selected_system.primary.radius_sol
        R2 = self.selected_system.secondary.radius_sol

        if self.rbSpecLum.isChecked():
            lspec1 = np.log10( Teff1**4 / 10**logg1 / sol_sp)
            lspec2 = np.log10( Teff2**4 / 10**logg2 / sol_sp)
        else:
            lspec1 = np.log10( (R1)**2 * (Teff1 / sun_teff)**4 )
            lspec2 = np.log10( (R2)**2 * (Teff2 / sun_teff)**4 )


        if self.tracks:

            main_sequence_phases = [0,1,2]
            for mass, self.track in self.tracks.items():
                if 'phase' in self.track.columns:
                    track_main_sequence = self.track[self.track['phase'].isin(main_sequence_phases)]
                else:
                    raise ValueError("Phase column not found in track data.")

                start_idx = track_main_sequence.index[0]
                start_teff = 10**(track_main_sequence.at[start_idx, 'log_Teff']+0.005)

                if self.rbSpecLum.isChecked():
                    start_lum = track_main_sequence.at[start_idx, 'log_L_spec']-0.02
                else:
                    start_lum = track_main_sequence.at[start_idx, 'log_L']-0.02

                x_limits = plt.gca().get_xlim()
                y_limits = plt.gca().get_ylim()

                if x_limits[0]*0.98 >= start_teff >= x_limits[1]*1.02 and y_limits[0]*1.03 <= start_lum <= y_limits[1]*0.97:
                    if mass < 8:
                        label = f'{mass:.1f} $M_\\odot$'
                    else:
                        label = f'{mass:.0f} $M_\\odot$'
                    ax.text(start_teff, start_lum, label, fontsize=6, ha='right', va='bottom')


        label = f'{M1:.2f} $M_\\odot$'
        ax.text(Teff1*0.95,lspec1*1.003, label, fontsize=6, ha='right', va='bottom', c="royalblue")

        label = f'{M2:.2f} $M_\\odot$'
        ax.text(Teff2*0.95,lspec2*1.003, label, fontsize=6, ha='right', va='bottom', c="royalblue")


    def on_axes_limits_change(self, ax, canvas, event_ax):
        self.update_labels(ax)
        canvas.draw()


    def plot_phot_hrd(self):

        M1 = self.selected_system.primary.mass
        M2 = self.selected_system.secondary.mass

        R1 = self.selected_system.primary.radius_sol
        R2 = self.selected_system.secondary.radius_sol

        Teff1 = self.selected_system.primary.Teff
        Teff2 = self.selected_system.secondary.Teff

        logg1 = self.selected_system.primary.logg
        logg2 = self.selected_system.secondary.logg

        lspec1 = np.log10( (R1)**2 * (Teff1 / sun_teff)**4 )
        lspec2 = np.log10( (R2)**2 * (Teff2 / sun_teff)**4 )

        fig, _ = plt.subplots(figsize=(5, 5))

        main_sequence_phases = [0,1,2]

        # Filter isochrones by phase
        if 'phase' in self.isochrones.columns:
            isochrones_main_sequence = self.isochrones[self.isochrones['phase'].isin(main_sequence_phases)]
        else:
            raise ValueError("Phase column not found in isochrones data.")

        unique_ages = sorted(isochrones_main_sequence['log10_isochrone_age_yr'].unique())
        print(unique_ages)
        filtered_ages = unique_ages[::3]
        print(filtered_ages)

        # Plot filtered isochrones
        for age in filtered_ages:
            subset = isochrones_main_sequence[isochrones_main_sequence['log10_isochrone_age_yr'] == age]
            plt.plot(10**subset['log_Teff'], subset['log_L'], '--k',lw=0.5,alpha=0.35, label=f'Age: {10**age:.0f} years')

        lum_min = min(lspec1, lspec2)
        lum_max = max(lspec1, lspec2)

        xlim = [min(Teff1, Teff2) - 0.1 * min(Teff1, Teff2), max(Teff1, Teff2) + 0.1 * max(Teff1, Teff2)]
        ylim = [lum_min - 0.05 * lum_min, lum_max + 0.03 * lum_max]

        if self.tracks:
            for mass, track in self.tracks.items():
                if 'phase' in track.columns:
                    track_main_sequence = track[track['phase'].isin(main_sequence_phases)]
                    if not track_main_sequence.empty:
                        track_lum_min = track_main_sequence['log_L'].min()
                        track_lum_max = track_main_sequence['log_L'].max()

                        if (track_lum_min <= lum_max and track_lum_max >= lum_min):
                            plt.plot(10**track_main_sequence['log_Teff'], track_main_sequence['log_L'], "-k", lw=0.5,
                                     label=f'Mass: {mass:.2f} $M_\\odot$')

                            start_idx = track_main_sequence.index[0]
                            start_teff = 10**(track_main_sequence.at[start_idx, 'log_Teff'])

                            xlim[1] = max(xlim[1], start_teff)

                            if xlim[0] * 1.02 <= start_teff <= xlim[1] * 0.98:
                                if mass < 8:
                                    label = f'{mass:.1f} $M_\\odot$'
                                else:
                                    label = f'{mass:.0f} $M_\\odot$'
                                plt.text(start_teff, track_main_sequence.at[start_idx, 'log_L'] - 0.02, label, fontsize=6, ha='right', va='bottom')


        plt.xlim(xlim)
        plt.ylim(ylim)


        l1, = plt.plot(Teff1,lspec1,"*b",markersize=12, alpha=0.8,mew=0)
        l1.set_markerfacecolor((0.404, 0.514, 1, 0.878))
        label = f'{M1:.2f} $M_\\odot$'
        plt.text(Teff1*0.9,lspec1*1.003, label, fontsize=6, ha='right', va='bottom', c="royalblue")

        l2, = plt.plot(Teff2,lspec2,"*b",markersize=12, alpha=0.8,mew=0)
        l2.set_markerfacecolor((0.404, 0.514, 1, 0.878))
        label = f'{M2:.2f} $M_\\odot$'
        plt.text(Teff2*0.9,lspec2*1.003, label, fontsize=6, ha='right', va='bottom', c="royalblue")

        plt.xlabel(r'$T_{eff}$ [K]')
        plt.ylabel(r'log $L$/$L_{\odot}$')
        plt.xscale('log')
        #plt.yscale('log')
        plt.gca().invert_xaxis()  # Hotter stars on the left

        return fig


    def plot_columns_against_gen(self, dfs, custom_values_dict=None, custom_bounds=None, top_minima_df=None, color=False):

        labels = self.mainconfig['params_opt'].copy()
        labels.append('M1')

        fig, axs = plt.subplots(len(labels), 1, sharex=True, figsize=(6, 10), dpi=150)
        cmap = cm.seismic

        if color:
            all_wsum_obj = np.concatenate([df['wsum_obj'].values for df in dfs])
            norm = LogNorm(vmin=all_wsum_obj.min(), vmax=all_wsum_obj.max())


        for i, label in enumerate(labels):
            # Plot each column against 'gen' in a subplot
            for j in range(len(dfs)):
                df = dfs[j]
                if color:
                    colors = cmap(norm(df['wsum_obj'].values))
                    sc = axs[i].scatter(df['gen'], df[label], c=colors, s=0.5, alpha=0.5)
                else:
                    axs[i].scatter(df['gen'], df[label], 0.5, "k", alpha=0.5)

                if custom_values_dict is not None:
                    axs[i].plot([0, max(df['gen'])], [custom_values_dict[label], custom_values_dict[label]], "-r", lw=1)

            if top_minima_df is not None:
                for _, row in top_minima_df.iterrows():
                    axs[i].plot([0, max(df['gen'])], [row[label], row[label]], "--r", lw=1)

            axs[i].set_ylabel(label)
            if label == "m2":
                axs[i].set_ylabel("q")

            if custom_bounds is not None:
                axs[i].set_ylim(custom_bounds[i])

        axs[-1].set_xlabel("Generation")
        plt.subplots_adjust(hspace=0.1)


        return fig

    def plot_solution(self, params_buff):

        input_params = dict(params_buff)
        if isinstance(input_params, dict):
            params = self.mainconfig['params_init'].copy()

            for param in input_params:
                if param in params:
                    params[param] = input_params[param]

        self.transformed_params = {
        "sbratio": 1.0,
        }



        orb = orbit(porb=self.mainconfig["params_init"]["porb"], inc_deg=self.mainconfig["params_init"]["incl"], t0=self.mainconfig["params_init"]["t0"],
                    a=self.mainconfig["params_init"]["a"], q=self.mainconfig["params_init"]["q"], f_c=self.mainconfig["params_init"]["f_c"], f_s=self.mainconfig["params_init"]["f_s"])

        system_obs = binary(orbit=orb, primary=star(), secondary=star(), sbratio=self.transformed_params["sbratio"], l3=self.mainconfig["params_init"]["l3"], gamma=self.mainconfig["params_init"]["gamma"],
                            obsLC=self.lcs, obsSpec=self.specs, nbins=self.mainconfig["nbins"],  gdc_option="opt", phase1=self.mainconfig["phase1"], phase2=self.mainconfig["phase2"], oversampling=self.mainconfig["oversampling"], width=self.mainconfig["width"])

        system_obs.UPDATE_PARAMS(params)



        synthV_instance = self.synthV_core1
        synthV_imu_instance = self.synthV_imu_core1

        system_obs = repair.update_system(
            self.mainconfig, system_obs, params, synthV_instance, synthV_imu_instance, ifplot_main=False, if_save_mods=False, passband=self.passband)

        self.selected_system = system_obs


        ffi = interp1d(
            system_obs.primary.spectrum[:, 0], system_obs.primary.spectrum[:, 1], fill_value="extrapolate")
        mm = ffi(system_obs.obsSpec.wl_eqlog_to_wl)
        ob1 = repair.fetch_spec_new(system_obs.obsSpec.wl_eqlog_to_wl,
                             system_obs.obsSpec.mod[0] + 0.5, system_obs.obsSpec.wl_eqlog_to_wl, mm)


        ffi = interp1d(
            system_obs.secondary.spectrum[:, 0], system_obs.secondary.spectrum[:, 1], fill_value="extrapolate")
        mm = ffi(system_obs.obsSpec.wl_eqlog_to_wl)
        ob2 = repair.fetch_spec_new(system_obs.obsSpec.wl_eqlog_to_wl,
                             system_obs.obsSpec.mod[1] + 0.5, system_obs.obsSpec.wl_eqlog_to_wl, mm)

        system_obs.obsSpec.mod[0] = ob1
        system_obs.obsSpec.mod[1] = ob2


        def plot_spectrum(ax, wl_range, ob, system_spectrum, system_wl_eqlog):
            indices = np.where((system_wl_eqlog >= wl_range[0]) & (
                system_wl_eqlog <= wl_range[1]))
            ax.plot(system_wl_eqlog[indices],
                    ob[indices], lw=0.5, c='cornflowerblue')

            system_spectrum = np.array(system_spectrum)
            spectrum_indices = np.where((system_spectrum[:, 0] >= wl_range[0]) & (
                system_spectrum[:, 0] <= wl_range[1]))
            ax.plot(system_spectrum[spectrum_indices, 0].flatten(
            ), system_spectrum[spectrum_indices, 1].flatten(), lw=0.5, c="crimson")

            ax.set_xlim(wl_range)

        fig = plt.figure(figsize=(4, 6), dpi=150)  # Adjusted figure size and DPI
        plt.subplots_adjust(hspace=0.5, wspace=0.3)  # Adjusted spacing between subplots

        extra_params = {
            'axes.titlesize': 6,
            'axes.labelsize': 6,
            'xtick.labelsize': 6,
            'ytick.labelsize': 6,
            'legend.fontsize': 6,
            'axes.linewidth': 0.5,
            'grid.linewidth': 0.5,
            'lines.linewidth': 0.5,
            'xtick.major.width': 0.5,
            'xtick.minor.width': 0.3,
            'ytick.major.width': 0.5,
            'ytick.minor.width': 0.3,
            'xtick.major.size': 3,
            'xtick.minor.size': 1.5,
            'ytick.major.size': 3,
            'ytick.minor.size': 1.5,
            'figure.subplot.top': 0.99,
            'figure.subplot.bottom': 0.05,
            'figure.subplot.left': 0.1,
            'figure.subplot.right': 0.95,
            'figure.facecolor':'#ececec',
            'axes.facecolor':'#ececec'
        }
        plt.rcParams.update(extra_params)



        ax_lc = plt.subplot2grid(shape=(5, 2), loc=(0, 0), colspan=2)

        ax_dis = plt.subplot2grid(shape=(5, 2), loc=(1, 0), colspan=2)

        ax1 = plt.subplot2grid(shape=(5, 2), loc=(2, 0), colspan=1)
        ax2 = plt.subplot2grid(shape=(5, 2), loc=(2, 1), colspan=1)
        ax3 = plt.subplot2grid(shape=(5, 2), loc=(3, 0), colspan=1)
        ax4 = plt.subplot2grid(shape=(5, 2), loc=(3, 1), colspan=1)
        ax5 = plt.subplot2grid(shape=(5, 2), loc=(4, 0), colspan=1)
        ax6 = plt.subplot2grid(shape=(5, 2), loc=(4, 1), colspan=1)

        ax_lc.plot(system_obs.obsLC.binned_phase-1, system_obs.obsLC.binned_flux,
                   marker='.', markersize=1, c='cornflowerblue', alpha=0.5)
        ax_lc.plot(system_obs.obsLC.binned_phase-1,
                   system_obs.modLC, lw=0.5, c="crimson")
        ax_lc.plot(system_obs.obsLC.binned_phase, system_obs.obsLC.binned_flux,
                   marker='.', markersize=1, c='cornflowerblue', alpha=0.5)
        ax_lc.plot(system_obs.obsLC.binned_phase,
                   system_obs.modLC, lw=0.5, c="crimson")
        ax_lc.plot(system_obs.obsLC.binned_phase+1, system_obs.obsLC.binned_flux,
                   marker='.', markersize=1, c='cornflowerblue', alpha=0.5)
        ax_lc.plot(system_obs.obsLC.binned_phase+1,
                   system_obs.modLC, lw=0.5, c="crimson")
        ax_lc.set_xlim([-0.4, 1.4])


        shift = (np.average(
            np.abs(system_obs.residSpec[:, 0]-system_obs.residSpec[:, -1])))*1.5
        for i in range(len(system_obs.residSpec[0])):

            xid = int((len(system_obs.obsSpec.wl_eqlog_to_wl)-len(system_obs.residSpec[:, i]))/2.0)
            #print(len(system_obs.obsSpec.wl_eqlog_to_wl), len(system_obs.residSpec[:, i]), xid)
            c = "cornflowerblue"
            if i % 2 == 0:
                c = "crimson"
            ax_dis.plot(system_obs.obsSpec.wl_eqlog_to_wl[xid:-xid], system_obs.residSpec[:, i] -
                        shift*i, lw=0.5, c=c, alpha=0.5)

        ax_dis.set_ylim([-(len(system_obs.residSpec[0])+1)*shift, shift])
        ax_dis.set_xlim([system_obs.obsSpec.wl_eqlog_to_wl[3000], system_obs.obsSpec.wl_eqlog_to_wl[-xid]])
        #ax_dis.set_xticks([])

        plot_spectrum(ax1, [self.mainconfig["wl11"], self.mainconfig["wl21"]],
                      ob1, system_obs.primary.spectrum, system_obs.obsSpec.wl_eqlog_to_wl)
        plot_spectrum(ax3, [4300, 4380], ob1,
                      system_obs.primary.spectrum, system_obs.obsSpec.wl_eqlog_to_wl)
        plot_spectrum(ax5, [4450, 4490], ob1,
                      system_obs.primary.spectrum, system_obs.obsSpec.wl_eqlog_to_wl)

        plot_spectrum(ax2, [self.mainconfig["wl12"], self.mainconfig["wl22"]],
                      ob2, system_obs.secondary.spectrum, system_obs.obsSpec.wl_eqlog_to_wl)
        plot_spectrum(ax4, [4300, 4380], ob2,
                      system_obs.secondary.spectrum, system_obs.obsSpec.wl_eqlog_to_wl)
        plot_spectrum(ax6, [4450, 4490], ob2,
                      system_obs.secondary.spectrum, system_obs.obsSpec.wl_eqlog_to_wl)

        self.selected_system = system_obs

        return fig


    def plot_objectives_convergence(self,dfs):
        
        mainconfig = self.mainconfig
        objective_type = mainconfig.get('objective_type', 'chi2') # chi2 = default value if key is not specified in config file

        fig, axs = plt.subplots(ncols=3,nrows=4,sharex=True,figsize=(6,5),dpi=150)
        for j in range(len(dfs)):

            df = dfs[j]

            ngen = dfs[j]['gen'].values[0]

            df['sum_obj'] = df['obj1']+df['obj2']+df['obj3']+df['obj4']
            min_sum_row = df.loc[df['sum_obj'].idxmin()]
            max_sum_row = df.loc[df['sum_obj'].idxmax()]

            if objective_type == 'logP':
                axs[0][0].set_title('Max separate')
                axs[0][0].scatter(ngen, max(df['obj1']), 1, "k")
                axs[1][0].scatter(ngen, max(df['obj2']), 1, "k")
                axs[2][0].scatter(ngen, max(df['obj3']), 1, "k")
                axs[3][0].scatter(ngen, max(df['obj4']), 1, "k")

                axs[0][1].set_title('Max coupled')
                axs[0][1].scatter(ngen, max_sum_row['obj1'], 1, "k")
                axs[1][1].scatter(ngen, max_sum_row['obj2'], 1, "k")
                axs[2][1].scatter(ngen, max_sum_row['obj3'], 1, "k")
                axs[3][1].scatter(ngen, max_sum_row['obj4'], 1, "k")

                axs[0][2].set_title('Median')
                axs[0][2].scatter(ngen, np.median(df['obj1']), 1, "k")
                axs[1][2].scatter(ngen, np.median(df['obj2']), 1, "k")
                axs[2][2].scatter(ngen, np.median(df['obj3']), 1, "k")
                axs[3][2].scatter(ngen, np.median(df['obj4']), 1, "k")
            else:
                axs[0][0].set_title('Min separate')
                axs[0][0].scatter(ngen, min(df['obj1']), 1, "k")
                axs[1][0].scatter(ngen, min(df['obj2']), 1, "k")
                axs[2][0].scatter(ngen, min(df['obj3']), 1, "k")
                axs[3][0].scatter(ngen, min(df['obj4']), 1, "k")
    
                axs[0][1].set_title('Min coupled')
                axs[0][1].scatter(ngen, min_sum_row['obj1'], 1, "k")
                axs[1][1].scatter(ngen, min_sum_row['obj2'], 1, "k")
                axs[2][1].scatter(ngen, min_sum_row['obj3'], 1, "k")
                axs[3][1].scatter(ngen, min_sum_row['obj4'], 1, "k")
    
                axs[0][2].set_title('Median')
                axs[0][2].scatter(ngen, np.median(df['obj1']), 1, "k")
                axs[1][2].scatter(ngen, np.median(df['obj2']), 1, "k")
                axs[2][2].scatter(ngen, np.median(df['obj3']), 1, "k")
                axs[3][2].scatter(ngen, np.median(df['obj4']), 1, "k")

        return fig


    def read_all_gens(self, ngen, stepgen):

        listdf = []
        for ngen in range(0,ngen,stepgen):
            mainconfig = self.mainconfig
            
            
            objective_type = mainconfig.get('objective_type', 'chi2') # chi2 = default value if key is not specified in config file
            
            try:
                params_loc = self.mainconfig['params_init']
                data = pd.read_csv(f"{mainconfig['saveto']}/gen{ngen}_initPorb{params_loc['porb']}_initA{params_loc['a']}_ngen{mainconfig['ngen']}_popsize{mainconfig['popsize']}.dat",sep="\t")

                data = data.apply(pd.to_numeric, errors='coerce')

                if objective_type == 'logP':
                    data = data[data["obj1"]> -1e10]
                    data = data[data["obj2"]> -1e10]
                    data = data[data["obj3"]> -1e10]
                    data = data[data["obj4"]> -1e10]
                else:
                    data = data[data["obj1"]<1000.0]
                    data = data[data["obj2"]<1000.0]
                    data = data[data["obj3"]<1000.0]
                    data = data[data["obj4"]<1000.0]


                if "q" not in data.columns:
                    data = data.rename(columns={"m2": "q"})

                    try:
                        data = data.rename(columns={"Ve1": "Ve1sini"})
                        data = data.rename(columns={"Ve2": "Ve2sini"})
                    except:
                        print("")

                data["sum_obj"] = data["obj1"] +  data["obj2"] + data["obj3"] + data["obj4"]

                data['r1'] = data['r1'].apply(lambda x: np.abs(x))
                data['r2'] = data['r2'].apply(lambda x: np.abs(x))
                data['q'] = data['q'].apply(lambda x: np.abs(x))
                data['a'] = data['a'].apply(lambda x: np.abs(x))


                data = data[data['r1'] > 0]
                data = data[data['r1'] < 1]
                data = data[data['r2'] > 0]
                data = data[data['r2'] < 1]
                data = data[data['incl'] <= 92.0]
                data = data[data['q'] < 3.0]

                porb = self.mainconfig['params_init']['porb']


                M1_all = []
                for i,row in data.iterrows():

                    ecc = 0.0
                    if 'f_c' in data.columns:
                        ecc, omega = SOLVE_ECC_OMEGA(row['f_c'],row['f_s'])


                    a1 = row["a"]/(1.+(1/row["q"]))
                    a2 = row["a"]/(1.+row["q"])
                    incl = row["incl"]
                    K1 = CALCULATE_SEMI_AMPLITUDE(a1, porb, incl, ecc)
                    K2 = CALCULATE_SEMI_AMPLITUDE(a2, porb, incl, ecc)
                    M1,M2 = CALCULATE_MASSES(K1, K2, porb, incl, ecc)
                    M1_all.append(M1)

                data["M1"] = M1_all

                data['gen'] = ngen


                population = df_to_deap_population(data, labels=self.mainconfig['params_opt'], sort_by_fitness=True)
                data = deap_population_to_df(population, labels=self.mainconfig['params_opt'])

                #print(data)

                listdf.append(data)

            except Exception as e:
                print(e)
        return listdf


    def plot_corner_filtered(self, df_full, custom_bounds = None, custom_min = None, ngen=None,fpath=None, cutoff=1):
        
        
        try:
            if custom_min is None:
                custom_min = self.custom_min
            if custom_bounds is None:
                custom_bounds = self.custom_bounds
    
            df_full = pd.concat(df_full[-10:])
            
            cut_obj1 = min(df_full['obj1']) + cutoff*(max(df_full['obj1'])-min(df_full['obj1']))
            cut_obj2 = min(df_full['obj2']) + cutoff*(max(df_full['obj2'])-min(df_full['obj2']))
            cut_obj3 = min(df_full['obj3']) + cutoff*(max(df_full['obj3'])-min(df_full['obj3']))
            cut_obj4 = min(df_full['obj4']) + cutoff*(max(df_full['obj4'])-min(df_full['obj4']))
            
            df_full = df_full[df_full['obj1']<= cut_obj1]
            df_full = df_full[df_full['obj2']<= cut_obj2]
            df_full = df_full[df_full['obj3']<= cut_obj3]
            df_full = df_full[df_full['obj4']<= cut_obj4]
            
            
            try:
                df_full=df_full.drop(columns=['rating'])
            except:
                print("")
            try:
                df_full=df_full.drop(columns=['Unnamed: 0'])
            except:
                print('')
            try:
                df_full=df_full.drop(columns=['gen','obj1',
                                              'obj2', 'obj3', 'obj4', 'sum_obj', 'wsum_obj' ])
            except:
                df_full=df_full.drop(columns=['gen','obj1',
                                              'obj2', 'obj3', 'obj4', 'sum_obj', ])
            n_cols = len(df_full.columns)
    
    
            fig, axs = plt.subplots(n_cols, n_cols, figsize=(2*n_cols, 2*n_cols), dpi=150)
    
            sigmas = {}
    
    
            for i in range(n_cols):
                for j in range(n_cols):
    
                    df = df_full
                    if custom_bounds is not None:
                        df = df[df.iloc[:, j] > custom_bounds[j][0]]
                        df = df[df.iloc[:, j] < custom_bounds[j][1]]
                        df = df[df.iloc[:, i] > custom_bounds[i][0]]
                        df = df[df.iloc[:, i] < custom_bounds[i][1]]
    
                    ax = axs[i, j]
    
                    if i == j:  # Histogram on the diagonal
                        ax.hist(df.iloc[:, i], bins=20, color='gray', alpha=0.7)
                        ax.set_yticklabels([])
                        ax.set_yticks([])
                        if custom_bounds is not None:
                            ax.set_xlim(custom_bounds[j])
    
                    elif i < j:  # Empty plots above the diagonal
                        ax.axis('off')
    
                    elif i > j:  # Scatter plot and confidence ellipses below the diagonal
                        x = df.iloc[:, j]
                        y = df.iloc[:, i]
    
    
                        ax.scatter(x, y, s=1, color='k', alpha=0.5)
    
                        try:
                            _,pearson,sigma_x,sigma_y = confidence_ellipse(x, y, ax, n_std=1, edgecolor='red', facecolor='none', label='1σ')
                            sigmas[(j, i)] = (sigma_x, sigma_y)
                        except:
                            print("")
                        try:
                            confidence_ellipse(x, y, ax, n_std=2, edgecolor='orange', linestyle='--', label='2σ')
                        except:
                            print("")
                        try:
                            confidence_ellipse(x, y, ax, n_std=3, edgecolor='yellow', linestyle=':', label='3σ')
                        except:
                            print("")
    
    
                        if custom_bounds is not None:
                            ax.set_xlim(custom_bounds[j])
                            ax.set_ylim(custom_bounds[i])
    
    
                            ax.axvline(custom_min[j], color="g",lw=0.5)
                            ax.axhline(custom_min[i], color="g",lw=0.5)
                            ax.plot(custom_min[j], custom_min[i], 3, "sg")
    
    
                    if i == n_cols - 1:
                        ax.set_xlabel(df_full.columns[j])
                    else:
                        ax.set_xticklabels([])
                        ax.set_xticks([])
                    if j == 0:
                        ax.set_ylabel(df_full.columns[i])
                    else:
                        ax.set_yticklabels([])
                        ax.set_yticks([])
    
    
            if custom_min is not None:
                for (j, i), (sigma_x, sigma_y) in sigmas.items():
                    if i != j:  # Avoid the diagonal
                        truex = custom_min[j]
                        truey = custom_min[i]
                        ax = axs[j, i]  # Place text in the cell at (i, j) in the grid
                        sigma_text = f'1σ {df_full.columns[j]} = {100*sigma_x/truex:.1f}%\n1σ {df_full.columns[i]} = {100*sigma_y/truey:.1f}%'
                        ax.text(0.5, 0.5, sigma_text, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=14)
                        ax.axis('off')
    
    
            plt.tight_layout()
            return fig
        except Exception as e:
            print(f"!!! Error filtering cutoff {cutoff} !!!")
            print(e)
            return

    def plot_pareto_front(self, dfs, labels=None):

        if labels is None:
            labels = self.mainconfig['params_opt']
        
        objective_type = self.mainconfig.get('objective_type', 'chi2') # chi2 = default value if key is not specified in config file


        fig, ax = plt.subplots(3,4,figsize=(12,8), dpi=150)

        extra_params = {
            'axes.labelsize': 6,
            'xtick.labelsize': 6,
            'ytick.labelsize': 6,
            'legend.fontsize': 6,
            'axes.linewidth': 0.5,
            'grid.linewidth': 0.5,
            'lines.linewidth': 0.5,
            'xtick.major.width': 0.5,
            'xtick.minor.width': 0.3,
            'ytick.major.width': 0.5,
            'ytick.minor.width': 0.3,
            'xtick.major.size': 3,
            'xtick.minor.size': 1.5,
            'ytick.major.size': 3,
            'ytick.minor.size': 1.5,
            'figure.facecolor':'#ececec',
            'axes.facecolor':'#ececec'
        }
        plt.rcParams.update(extra_params)
        plt.subplots_adjust(hspace=0.8, wspace=0.8)

        for j in range(len(dfs)):

            population = df_to_deap_population(dfs[j], labels)
            pareto_front = get_first_pareto_front(population)

            ngen=dfs[j]['gen'].values[0]

            objectives = np.array([ind.fitness.values for ind in pareto_front])


            ax[0][0].plot(objectives[:,1],objectives[:,0],'ok',markersize=1, alpha = 1.0/(self.mainconfig['ngen']) * ngen)

            ax[0][1].plot(objectives[:,1],objectives[:,2],'ok',markersize=1, alpha = 1.0/(self.mainconfig['ngen']) * ngen)

            ax[0][2].plot(objectives[:,1],objectives[:,3],'ok',markersize=1, alpha = 1.0/(self.mainconfig['ngen']) * ngen)

            ax[0][3].plot(objectives[:,0],objectives[:,2],'ok',markersize=1, alpha = 1.0/(self.mainconfig['ngen']) * ngen)

            if len(dfs)-10 < j < len(dfs):


                ax[1][0].plot(objectives[:,1],objectives[:,0],'ok',markersize=1, alpha = 1.0/(self.mainconfig['ngen']) * ngen)

                ax[1][1].plot(objectives[:,1],objectives[:,2],'ok',markersize=1, alpha = 1.0/(self.mainconfig['ngen']) * ngen)

                ax[1][2].plot(objectives[:,1],objectives[:,3],'ok',markersize=1, alpha = 1.0/(self.mainconfig['ngen']) * ngen)

                ax[1][3].plot(objectives[:,0],objectives[:,2],'ok',markersize=1, alpha = 1.0/(self.mainconfig['ngen']) * ngen)



                ax[2][0].plot(objectives[:,1],objectives[:,0],'ok',markersize=1, alpha = 1.0/(self.mainconfig['ngen']) * ngen)

                ax[2][1].plot(objectives[:,1],objectives[:,2],'ok',markersize=1, alpha = 1.0/(self.mainconfig['ngen']) * ngen)

                ax[2][2].plot(objectives[:,1],objectives[:,3],'ok',markersize=1, alpha = 1.0/(self.mainconfig['ngen']) * ngen)

                ax[2][3].plot(objectives[:,0],objectives[:,2],'ok',markersize=1, alpha = 1.0/(self.mainconfig['ngen']) * ngen)


        ax[0][0].set_xlabel(f'$\chi^2$ sp.sep.')
        ax[0][0].set_ylabel(f'$\chi^2$ LC')

        ax[0][1].set_xlabel(f'$\chi^2$ sp.sep.')
        ax[0][1].set_ylabel(f'$\chi^2$ syn1')

        ax[0][2].set_xlabel(f'$\chi^2$ sp.sep.')
        ax[0][2].set_ylabel(f'$\chi^2$ syn2')

        ax[0][3].set_xlabel(f'$\chi^2$ LC')
        ax[0][3].set_ylabel(f'$\chi^2$ syn1')

        ax[1][0].set_xlabel(f'$\chi^2$ sp.sep.')
        ax[1][0].set_ylabel(f'$\chi^2$ LC')

        ax[1][1].set_xlabel(f'$\chi^2$ sp.sep.')
        ax[1][1].set_ylabel(f'$\chi^2$ syn1')

        ax[1][2].set_xlabel(f'$\chi^2$ sp.sep.')
        ax[1][2].set_ylabel(f'$\chi^2$ syn2')

        ax[1][3].set_xlabel(f'$\chi^2$ LC')
        ax[1][3].set_ylabel(f'$\chi^2$ syn1')

        if objective_type == 'chi2':
            ax[2][0].set_xscale('log')
            ax[2][0].set_yscale('log')
    
            ax[2][1].set_xscale('log')
            ax[2][1].set_yscale('log')
    
            ax[2][2].set_xscale('log')
            ax[2][2].set_yscale('log')
    
            ax[2][3].set_xscale('log')
            ax[2][3].set_yscale('log')


        ax[2][0].set_xlabel(f'$\chi^2$ sp.sep.')
        ax[2][0].set_ylabel(f'$\chi^2$ LC')

        ax[2][1].set_xlabel(f'$\chi^2$ sp.sep.')
        ax[2][1].set_ylabel(f'$\chi^2$ syn1')

        ax[2][2].set_xlabel(f'$\chi^2$ sp.sep.')
        ax[2][2].set_ylabel(f'$\chi^2$ syn2')

        ax[2][3].set_xlabel(f'$\chi^2$ LC')
        ax[2][3].set_ylabel(f'$\chi^2$ syn1')


        return fig


    def load_config(self):

        filepath = self.filepath



        if filepath:
            with open(filepath, 'r') as file:
                self.mainconfig = json.load(file)


            passband = np.loadtxt("tess_transmission.txt", skiprows=6,delimiter=",")
            passband[:,0] *= 10.0
            if "custom_passband" in self.mainconfig:
                passband = np.loadtxt(self.mainconfig["custom_passband"], delimiter="\t")
            self.passband = passband

            self.dfs_gens = self.read_all_gens(self.mainconfig["ngen"], self.stepgen)

            if "q" not in self.mainconfig['params_init']:

                self.mainconfig['params_init'] = dict(('q', v) if k == 'm2' else (k, v) for k, v in self.mainconfig['params_init'].items())
                self.mainconfig['params_bounds'] = dict(('q', v) if k == 'm2' else (k, v) for k, v in self.mainconfig['params_bounds'].items())

                for jj in range(len(self.mainconfig['params_opt'])):
                    if self.mainconfig['params_opt'][jj] == "m2":
                        self.mainconfig['params_opt'][jj] = "q"

            if "Ve1sini" not in self.mainconfig['params_init']:

                self.mainconfig['params_init'] = dict(('Ve1sini', v*np.sin(self.mainconfig['params_init']["incl"]/180.0*np.pi)) if k == 'Ve1' else (k, v) for k, v in self.mainconfig['params_init'].items())
                self.mainconfig['params_init'] = dict(('Ve2sini', v*np.sin(self.mainconfig['params_init']["incl"]/180.0*np.pi)) if k == 'Ve2' else (k, v) for k, v in self.mainconfig['params_init'].items())
                self.mainconfig['params_bounds'] = dict(('Ve1sini', v) if k == 'Ve1' else (k, v) for k, v in self.mainconfig['params_bounds'].items())
                self.mainconfig['params_bounds'] = dict(('Ve2sini', v) if k == 'Ve2' else (k, v) for k, v in self.mainconfig['params_bounds'].items())

                for jj in range(len(self.mainconfig['params_opt'])):
                    if self.mainconfig['params_opt'][jj] == "Ve1":
                        self.mainconfig['params_opt'][jj] = "Ve1sini"
                    if self.mainconfig['params_opt'][jj] == "Ve2":
                        self.mainconfig['params_opt'][jj] = "Ve2sini"

            print(len(self.dfs_gens))
            population = df_to_deap_population(self.dfs_gens[-1], labels=self.mainconfig['params_opt'])
            weights = calculate_weights_from_pareto_front(population)
            #weights=[1000,1,0.1,0.01]
            print(weights)
            self.weights=weights

            for j in range(len(self.dfs_gens)):
                self.dfs_gens[j]['wsum_obj'] = self.dfs_gens[j]['obj1']*weights[0]+self.dfs_gens[j]['obj2']*weights[1]+self.dfs_gens[j]['obj3']*weights[2]+self.dfs_gens[j]['obj4']*weights[3]


            self.custom_min = [ self.mainconfig['params_init'][ip] for ip in self.mainconfig['params_opt']]
            self.custom_bounds = [ self.mainconfig['params_bounds'][ip] for ip in self.mainconfig['params_opt']]

            ecc, omega = SOLVE_ECC_OMEGA(self.mainconfig['params_init']['f_c'],self.mainconfig['params_init']['f_s'])



            a1 = self.mainconfig['params_init']["a"]/(1.+(1/self.mainconfig['params_init']["q"]))
            a2 = self.mainconfig['params_init']["a"]/(1.+self.mainconfig['params_init']["q"])
            incl = self.mainconfig['params_init']["incl"]
            porb = self.mainconfig['params_init']['porb']
            K1 = CALCULATE_SEMI_AMPLITUDE(a1, porb, incl, ecc)
            K2 = CALCULATE_SEMI_AMPLITUDE(a2, porb, incl, ecc)
            M1,M2 = CALCULATE_MASSES(K1, K2, porb, incl, ecc)

            self.custom_min.append(M1)
            self.custom_bounds.append([M1*0.6,M1*1.4])

            print(self.custom_min , self.custom_bounds )
            
            specs = specSeries()
            if mainconfig["obsmode"] == "test":
                specs.LOAD_TESTSUIT(mainconfig["path_specs"], period = mainconfig["params_init"]["porb"], t0 = mainconfig["params_init"]["t0"], wl1 = mainconfig["wl11"], wl2 = mainconfig["wl21"])
            elif mainconfig["obsmode"] == "obs-list":
                specs.LOAD_LIST(mainconfig["path_specs"], wl1 = mainconfig["wl11"], wl2 = mainconfig["wl21"])
            else:
                specs.LOAD_HERMES(mainconfig["path_specs"], wl1 = mainconfig["wl11"], wl2 = mainconfig["wl21"])
            
            specs.TRIM(mainconfig["wl11"]+10, mainconfig["wl21"]-10)
            self.specs = specs


            lcf = np.loadtxt(self.mainconfig["path_lc"])
            if mainconfig["obsmode"] == "test":
                lcs = lightcurve(lcf[:,0], lcf[:,1], np.zeros_like(lcf[:,0]))
            else:
                lcs = lightcurve(lcf[:,0], lcf[:,1], lcf[:,2])
            self.lcs = lcs



            self.synthV_core1 = synthetic(wd_synthv="synthV_Imu_core1/", db_atmos_models=self.mainconfig["atmmodels"], if_pyastr=True,
                                     if_convolve=False, if_imu=False, atmos_mode="lin_interp", abund_tables=self.mainconfig["abunds"], pref=os.getcwd())
            self.synthV_imu_core1 = synthetic(wd_synthv="synthV_Imu_core1/", db_atmos_models=self.mainconfig["atmmodels"], if_convolve=False,
                                         if_imu=True, if_lines=False, atmos_mode="lin_interp", abund_tables=self.mainconfig["abunds"],  pref=os.getcwd())

            
    def plot_corner_filtered_last(self, df_full, custom_bounds = None, custom_min = None, ngen=None,fpath=None):


        if custom_min is None:
            custom_min = self.custom_min
        if custom_bounds is None:
            custom_bounds = self.custom_bounds

        df_full = pd.concat(df_full[-10:])
        try:
            df_full=df_full.drop(columns=['rating'])
        except:
            print("")
        try:
            df_full=df_full.drop(columns=['Unnamed: 0'])
        except:
            print('')
        try:
            df_full=df_full.drop(columns=['gen','obj1',
                                          'obj2', 'obj3', 'obj4', 'sum_obj', 'wsum_obj' ])
        except:
            df_full=df_full.drop(columns=['gen','obj1',
                                          'obj2', 'obj3', 'obj4', 'sum_obj', ])
        n_cols = len(df_full.columns)


        fig, axs = plt.subplots(n_cols, n_cols, figsize=(2*n_cols, 2*n_cols), dpi=150)

        sigmas = {}


        for i in range(n_cols):
            for j in range(n_cols):

                df = df_full
                if custom_bounds is not None:
                    df = df[df.iloc[:, j] > custom_bounds[j][0]]
                    df = df[df.iloc[:, j] < custom_bounds[j][1]]
                    df = df[df.iloc[:, i] > custom_bounds[i][0]]
                    df = df[df.iloc[:, i] < custom_bounds[i][1]]

                ax = axs[i, j]

                if i == j:  # Histogram on the diagonal
                    ax.hist(df.iloc[:, i], bins=20, color='gray', alpha=0.7)
                    ax.set_yticklabels([])
                    ax.set_yticks([])
                    if custom_bounds is not None:
                        ax.set_xlim(custom_bounds[j])

                elif i < j:  # Empty plots above the diagonal
                    ax.axis('off')

                elif i > j:  # Scatter plot and confidence ellipses below the diagonal
                    x = df.iloc[:, j]
                    y = df.iloc[:, i]


                    ax.scatter(x, y, s=1, color='k', alpha=0.5)

                    try:
                        _,pearson,sigma_x,sigma_y = confidence_ellipse(x, y, ax, n_std=1, edgecolor='red', facecolor='none', label='1σ')
                        sigmas[(j, i)] = (sigma_x, sigma_y)
                    except:
                        print("")
                    try:
                        confidence_ellipse(x, y, ax, n_std=2, edgecolor='orange', linestyle='--', label='2σ')
                    except:
                        print("")
                    try:
                        confidence_ellipse(x, y, ax, n_std=3, edgecolor='yellow', linestyle=':', label='3σ')
                    except:
                        print("")


                    if custom_bounds is not None:
                        ax.set_xlim(custom_bounds[j])
                        ax.set_ylim(custom_bounds[i])


                        ax.axvline(custom_min[j], color="g",lw=0.5)
                        ax.axhline(custom_min[i], color="g",lw=0.5)
                        ax.plot(custom_min[j], custom_min[i], 3, "sg")


                if i == n_cols - 1:
                    ax.set_xlabel(df_full.columns[j])
                else:
                    ax.set_xticklabels([])
                    ax.set_xticks([])
                if j == 0:
                    ax.set_ylabel(df_full.columns[i])
                else:
                    ax.set_yticklabels([])
                    ax.set_yticks([])


        if custom_min is not None:
            for (j, i), (sigma_x, sigma_y) in sigmas.items():
                if i != j:  # Avoid the diagonal
                    truex = custom_min[j]
                    truey = custom_min[i]
                    ax = axs[j, i]  # Place text in the cell at (i, j) in the grid
                    sigma_text = f'1σ {df_full.columns[j]} = {100*sigma_x/truex:.1f}%\n1σ {df_full.columns[i]} = {100*sigma_y/truey:.1f}%'
                    ax.text(0.5, 0.5, sigma_text, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=14)
                    ax.axis('off')


        plt.tight_layout()
        return fig

    def interpolate_mass_age_2d(self, Teff_array, L_array, if_Lphot = False):
        print(len(Teff_array),len(L_array))
        points = []
        masses = []
        ages = []

        for mm, track in self.tracks.items():
            if 'phase' in track.columns:
                track_main_sequence = track[track['phase'].isin([0, 1, 2])]
                if not track_main_sequence.empty:
                    teff_values = 10**track_main_sequence['log_Teff'].values
                    if if_Lphot:
                        lum_values = track_main_sequence['log_L'].values
                    else:
                        lum_values = track_main_sequence['log_L_spec'].values

                    age_values = track_main_sequence['star_age'].values
                    mass_values = track_main_sequence['star_mass'].values

                    for i in range(len(teff_values)):
                        points.append([teff_values[i], lum_values[i]])
                        masses.append(mass_values[i])
                        ages.append(age_values[i])

        points = np.array(points)
        masses = np.array(masses)
        ages = np.array(ages)


        Teff_array = np.asarray(Teff_array)
        L_array = np.asarray(L_array)


        age_interpolator = LinearNDInterpolator(points, ages)
        interpolated_age = age_interpolator(Teff_array, L_array)

        mass_interpolator = LinearNDInterpolator(points, masses)
        interpolated_mass = mass_interpolator(Teff_array, L_array)

        if np.any(np.isnan(interpolated_mass)):
            print("Warning: Some interpolations returned NaN. Some requested points might be out of the tracks' range.")

        return interpolated_mass, interpolated_age


    def find_evolutionary_pars(self, params_array, if_Lphot = False):

        array_Teffs_primary = []
        array_Lspec_primary = []
        array_Lphot_primary = []
        array_dyn_mass_primary = []
        array_evo_mass_primary = []
        array_evo_age_primary = []

        array_Teffs_secondary = []
        array_Lspec_secondary = []
        array_Lphot_secondary = []
        array_dyn_mass_secondary = []
        array_evo_mass_secondary = []
        array_evo_age_secondary = []

        orb = orbit(porb=self.mainconfig["params_init"]["porb"], inc_deg=self.mainconfig["params_init"]["incl"], t0=self.mainconfig["params_init"]["t0"],
                    a=self.mainconfig["params_init"]["a"], q=self.mainconfig["params_init"]["q"], f_c=self.mainconfig["params_init"]["f_c"], f_s=self.mainconfig["params_init"]["f_s"])

        system_obs = binary(orbit=orb, primary=star(), secondary=star(), sbratio=1.0, l3=self.mainconfig["params_init"]["l3"], gamma=self.mainconfig["params_init"]["gamma"],
                            obsLC=self.lcs, obsSpec=self.specs, nbins=self.mainconfig["nbins"],  gdc_option="opt", phase1=self.mainconfig["phase1"], phase2=self.mainconfig["phase2"], oversampling=self.mainconfig["oversampling"], width=self.mainconfig["width"])


        try:
            for jj in range(len(params_array)):

                input_params = dict(params_array[jj])
                if isinstance(input_params, dict):
                    params = self.mainconfig['params_init'].copy()

                    for param in input_params:
                        if param in params:
                            params[param] = input_params[param]


                system_obs.UPDATE_PARAMS(params)

                Teff1 = system_obs.primary.Teff
                Teff2 = system_obs.secondary.Teff
                logg1 = system_obs.primary.logg
                logg2 = system_obs.secondary.logg

                loglspec1 = np.log10( Teff1**4 / 10**logg1 / sol_sp)
                loglspec2 = np.log10( Teff2**4 / 10**logg2 / sol_sp)

                R1 = system_obs.primary.radius_sol
                R2 = system_obs.secondary.radius_sol

                loglphot1 = np.log10( (R1)**2 * (Teff1 / sun_teff)**4 )
                loglphot2 = np.log10( (R2)**2 * (Teff2 / sun_teff)**4 )

                M1_dyn = system_obs.primary.mass
                M2_dyn = system_obs.secondary.mass

                array_Teffs_primary.append(Teff1)
                array_Lspec_primary.append(loglspec1)
                array_Lphot_primary.append(loglphot1)
                array_dyn_mass_primary.append(M1_dyn)

                array_Teffs_secondary.append(Teff2)
                array_Lspec_secondary.append(loglspec2)
                array_Lphot_secondary.append(loglphot2)
                array_dyn_mass_secondary.append(M2_dyn)


            if if_Lphot:
                array_evo_mass_primary, array_evo_age_primary = self.interpolate_mass_age_2d(array_Teffs_primary, array_Lphot_primary, if_Lphot=True)
                array_evo_mass_secondary, array_evo_age_secondary  = self.interpolate_mass_age_2d(array_Teffs_secondary, array_Lphot_secondary, if_Lphot=True)
            else:
                array_evo_mass_primary, array_evo_age_primary = self.interpolate_mass_age_2d(array_Teffs_primary, array_Lphot_primary, if_Lphot=False)
                array_evo_mass_secondary, array_evo_age_secondary = self.interpolate_mass_age_2d(array_Teffs_secondary, array_Lphot_secondary, if_Lphot=False)

            result = {}

            result['primary_Teffs'] = array_Teffs_primary
            result['psecondary_Teffs'] = array_Teffs_secondary

            if if_Lphot:
                result['primary_L'] = array_Lphot_primary
                result['secondary_L'] = array_Lphot_secondary
            else:
                result['primary_L'] = array_Lspec_primary
                result['secondary_L'] = array_Lspec_secondary

            result['primary_mass_evo'] = array_evo_mass_primary
            result['secondary_mass_evo'] = array_evo_mass_secondary
            result['primary_age_evo'] = array_evo_age_primary
            result['secondary_age_evo'] = array_evo_age_secondary
            result['primary_mass_dyn'] = array_dyn_mass_primary
            result['secondary_mass_dyn'] = array_dyn_mass_secondary

            return result

        except Exception as e:
            print(e)
            return None


    def plot_corner_alpha(self, df_full, custom_bounds=None, custom_min=None, ngen=None, fpath=None, noaxes=False, extra_params = None, selector = False):


        if selector:
            if len(df_full) > 30:
                df_full = pd.concat(df_full[::10])
            else:
                df_full = pd.concat(df_full)
            if custom_bounds is None:
                custom_bounds = self.custom_bounds
        else:
            df_full = pd.concat(df_full)

        df_full = df_full.drop_duplicates(subset=self.mainconfig['params_opt'])

        full_population = df_to_deap_population(df_full, labels=self.mainconfig['params_opt'], sort_by_fitness=True)
        df_full = deap_population_to_df(full_population, labels=self.mainconfig['params_opt'], if_add_columns=True, weights=self.weights)

        if self.rbByRating.isChecked():
            gens = df_full['rating']
        else:
            gens = df_full['gen']

        top_5_percent_count = int(0.05 * len(df_full))
        top_5_percent_df = df_full.head(top_5_percent_count)
        top_5_percent_indices = top_5_percent_df.index

        weights = gens / np.sum(gens)

        try:
            df_full=df_full.drop(columns=['rating'])
        except:
            print("")
        try:
            df_full=df_full.drop(columns=['Unnamed: 0'])
        except:
            print('')
        try:
            df_full = df_full.drop(columns=['gen', 'obj1', 'obj2', 'obj3', 'obj4', 'sum_obj', 'wsum_obj'])
        except:
            df_full = df_full.drop(columns=['gen', 'obj1', 'obj2', 'obj3', 'obj4', 'sum_obj'])
        n_cols = len(df_full.columns)

        if extra_params is not None:
            plt.rcParams.update(extra_params)
        else:
            plt.rcParams.update({'figure.figsize': (2*n_cols, 2*n_cols)})

        fig, axs = plt.subplots(n_cols, n_cols)
        fig.tight_layout()
        plt.subplots_adjust(hspace=0.0, wspace=0.0)  # Adjust the spacing between subplots


        xlims = [None] * n_cols
        ylims = [None] * n_cols

        sigma_values = {i: None for i in range(n_cols)}

        for i in range(n_cols):
            for j in range(n_cols):
                x = df_full.iloc[:, j]
                y = df_full.iloc[:, i]

                mean_x = np.average(x, weights=weights)
                mean_y = np.average(y, weights=weights)
                center = (mean_x, mean_y)
                cov = np.cov(x, y, aweights=weights)
                sigma_x = np.sqrt(cov[0, 0])
                sigma_y = np.sqrt(cov[1, 1])

                if not selector:
                    xlims[j] = [center[0] - 1.3 * sigma_x, center[0] + 1.3 * sigma_x]
                    ylims[i] = [center[1] - 1.3 * sigma_y, center[1] + 1.3 * sigma_y]
                else:
                    xlims[j] = [custom_bounds[j][0], custom_bounds[j][1]]
                    ylims[i] = [custom_bounds[i][0], custom_bounds[i][1]]
                sigma_values[i] = sigma_y

        scatter_plots = []

        if self.rbConstrMass.isChecked() or self.rbConstrAge.isChecked() or self.rbConstrMassSum.isChecked():
            params_result = self.find_evolutionary_pars(df_full.to_dict(orient='records'), if_Lphot=True)
            if self.rbConstrMass.isChecked() :
                diff_evo = np.abs(np.array(params_result['primary_mass_evo']) - np.array(params_result['primary_mass_dyn']))
            elif self.rbConstrMassSum.isChecked() :
                diff_evo = np.abs(np.array(params_result['primary_mass_evo']) - np.array(params_result['primary_mass_dyn']))
                diff_evo += np.abs(np.array(params_result['secondary_mass_evo']) - np.array(params_result['secondary_mass_dyn']))
            else:
                diff_evo = np.abs(np.array(params_result['primary_age_evo']) - np.array(params_result['secondary_age_evo']))


            median = np.median(diff_evo)
            spread = median - np.min(diff_evo)
            threshold = median + spread
            filtered_data = diff_evo[diff_evo < threshold]
            vmax_filtered = np.max(filtered_data)

            norm = Normalize(vmin=np.min(diff_evo), vmax=vmax_filtered)
            cmap = plt.cm.rainbow

        for i in range(n_cols):
            for j in range(n_cols):
                df = df_full
                ax = axs[i, j]

                if i == j:  # Histogram on the diagonal
                    ax.hist(df.iloc[:, i], bins=25, color='gray', alpha=0.7, weights=weights, range=xlims[i])
                    ax.set_xlim(xlims[i])
                    ax.set_yticklabels([])
                    ax.set_yticks([])
                    if noaxes:
                        ax.set_xticks([])

                elif i < j:  # Empty plots above the diagonal
                    ax.axis('off')

                elif i > j:  # Scatter plot and confidence ellipses below the diagonal
                    x = df.iloc[:, j]
                    y = df.iloc[:, i]

                    if self.rbConstrMass.isChecked() or self.rbConstrAge.isChecked() or self.rbConstrMassSum.isChecked():
                        sc = ax.scatter(x, y, c=diff_evo, cmap=cmap, norm=norm, s=1, alpha=0.5)
                        scatter_plots.append(sc)
                    else:
                        sc = ax.scatter(x, y, s=1, color='k', alpha=1.0 / max(gens) * gens, marker='.', edgecolors='none')
                        scatter_plots.append(sc)


                    mean_x = np.average(x, weights=weights)
                    mean_y = np.average(y, weights=weights)
                    center = (mean_x, mean_y)
                    cov = np.cov(x, y, aweights=weights)
                    sigma_x = np.sqrt(cov[0, 0])
                    sigma_y = np.sqrt(cov[1, 1])

                    ax.set_xlim(xlims[j])
                    ax.set_ylim(ylims[i])


                    if custom_min is not None:
                        ax.axvline(custom_min[j], color="grey", lw=0.5)
                        ax.axhline(custom_min[i], color="grey", lw=0.5)

                    if not noaxes:
                        ax.text(0.95, 0.95, f"{center[0]:.3f} ± {sigma_x:.3f}\n"
                                            f"{center[1]:.3f} ± {sigma_y:.3f}",
                                verticalalignment='top', horizontalalignment='right',
                                transform=ax.transAxes, color='k', fontsize=2)

                    if i == n_cols - 1:
                        ax.set_xlabel(df_full.columns[j])
                    else:
                        ax.set_xticklabels([])
                        ax.set_xticks([])
                    if j == 0:
                        ax.set_ylabel(df_full.columns[i])
                    else:
                        ax.set_yticklabels([])
                        ax.set_yticks([])

                    if noaxes:
                        ax.set_xticks([])
                        ax.set_yticks([])



        if selector:
            cursor = mplcursors.cursor(scatter_plots, hover=False)

            @cursor.connect("add")
            def on_add(sel):
                index = sel.target.index
                selected_row = df_full.iloc[index]
                self.selected_row = selected_row
                config_output = {
                        'selected_row': self.selected_row.to_dict(),
                        'params_init': self.mainconfig['params_init']
                    }
                self.txtConfig.setPlainText(json.dumps(config_output, indent=4))

                formatted_text = "Selected point:\n" + "\n".join([f"{k}: {v:.3f}" for k, v in selected_row.items()])

                for text in fig.texts:
                    text.set_visible(False)

                text = fig.text(0.8, 0.5, formatted_text, ha='center', fontsize=8, wrap=True)

                for sc in scatter_plots:
                    offsets = sc.get_offsets()

                    facecolors = sc.get_facecolors()

                    # Update facecolors
                    new_facecolors = []
                    for i in range(len(offsets)):
                        if i == index:
                            new_facecolors.append((0.8627, 0.0784, 0.2353, 0.75))  # Highlight selected point in crimson
                        else:
                            if self.rbConstrMass.isChecked() or self.rbConstrAge.isChecked() or self.rbConstrMassSum.isChecked():
                                original_color = facecolors[i]
                                new_facecolors.append(original_color)
                            else:
                                if df_full.index[i] in top_5_percent_indices:
                                    new_facecolors.append((0, 0.596, 1, 1))  # Keep top 5% in blue
                                else:
                                    new_facecolors.append((0, 0, 0, 0.5))  # Others in black

                    sc.set_facecolor(new_facecolors)

                    #if index < len(offsets):
                    #    sc.set_facecolor([(0, 0, 1, 1) if i == index else (0, 0, 0, 0.5) for i in range(len(offsets))])
                    if self.rbConstrMass.isChecked() or self.rbConstrAge.isChecked() or self.rbConstrMassSum.isChecked():
                        sc.set_sizes([5 if i == index else 1 for i in range(len(offsets))])
                    else:
                        sc.set_sizes([100 if i == index else 10 for i in range(len(offsets))])
                    sc.figure.canvas.draw_idle()


                params = self.mainconfig["params_init"].copy()
                toadd_key = ["ldc_1","ldc_2","ld_1","ld_2","gdc_1","gdc_2","heat_1","heat_2","l3","bfac_1","bfac_2"]
                toadd_val = [None,None,"mugrid","mugrid",None,None,None,None,0.0,None,None]
                for kk in range(len(toadd_key)):
                    if toadd_key[kk] not in params:
                        params[toadd_key[kk]] = toadd_val[kk]

                for param in selected_row.index:
                    if param in params:
                        params[param] = selected_row[param]

                if params['ldc_1'] is not None:
                    params['ld_1'] = "lin"
                if params['ldc_2'] is not None:
                    params['ld_2'] = "lin"



                self.update_selected_plot(params)



        return fig

    def plot_corner_alpha_colour(self, df_full, custom_bounds=None, custom_min=None, ngen=None, fpath=None, noaxes=False, extra_params = None, selector = False, colour='obj1', cmap='seismic'):

        if selector:
            df_full = pd.concat(df_full[::10])
            if custom_bounds is None:
                custom_bounds = self.custom_bounds
        else:
            df_full = pd.concat(df_full)

        df_full = df_full.drop_duplicates(subset=self.mainconfig['params_opt'])

        full_population = df_to_deap_population(df_full, labels=self.mainconfig['params_opt'], sort_by_fitness=True)
        df_full = deap_population_to_df(full_population, labels=self.mainconfig['params_opt'], if_add_columns=True, weights=self.weights)

        if self.rbByRating.isChecked():
            gens = df_full['rating']
        else:
            gens = df_full['gen']

        weights = gens / np.sum(gens)

        try:
            df_full=df_full.drop(columns=['rating'])
        except:
            print("")

        try:
            df_full=df_full.drop(columns=['Unnamed: 0'])
        except:
            print('')

        sum_obj = df_full[colour].copy()
        norm = LogNorm(vmin=sum_obj.min(), vmax=sum_obj.max())

        try:
            df_full = df_full.drop(columns=['gen',  'obj1', 'obj2', 'obj3', 'obj4', 'sum_obj', 'wsum_obj'])
        except:
            df_full = df_full.drop(columns=['gen',  'obj1', 'obj2', 'obj3', 'obj4', 'sum_obj'])
        n_cols = len(df_full.columns)

        if extra_params is not None:
            plt.rcParams.update(extra_params)
        else:
            plt.rcParams.update({'figure.figsize': (2*n_cols, 2*n_cols)})

        fig, axs = plt.subplots(n_cols, n_cols)
        fig.tight_layout()
        plt.subplots_adjust(hspace=0.0, wspace=0.0)

        xlims = [None] * n_cols
        ylims = [None] * n_cols

        sigma_values = {i: None for i in range(n_cols)}

        for i in range(n_cols):
            for j in range(n_cols):
                x = df_full.iloc[:, j]
                y = df_full.iloc[:, i]

                mean_x = np.average(x, weights=weights)
                mean_y = np.average(y, weights=weights)
                center = (mean_x, mean_y)
                cov = np.cov(x, y, aweights=weights)
                sigma_x = np.sqrt(cov[0, 0])
                sigma_y = np.sqrt(cov[1, 1])

                if not selector:
                    xlims[j] = [center[0] - 1.3 * sigma_x, center[0] + 1.3 * sigma_x]
                    ylims[i] = [center[1] - 1.3 * sigma_y, center[1] + 1.3 * sigma_y]
                else:
                    xlims[j] = [custom_bounds[j][0], custom_bounds[j][1]]
                    ylims[i] = [custom_bounds[i][0], custom_bounds[i][1]]
                sigma_values[i] = sigma_y

        scatter_plots = []

        for i in range(n_cols):
            for j in range(n_cols):
                ax = axs[i, j]

                if i == j:  # Histogram on the diagonal
                    ax.hist(df_full.iloc[:, i], bins=25, color='gray', alpha=0.7, weights=weights, range=xlims[i])
                    ax.set_xlim(xlims[i])
                    ax.set_yticklabels([])
                    ax.set_yticks([])
                    if noaxes:
                        ax.set_xticks([])

                elif i < j:  # Empty plots above the diagonal
                    ax.axis('off')

                elif i > j:  # Scatter plot and confidence ellipses below the diagonal
                    x = df_full.iloc[:, j]
                    y = df_full.iloc[:, i]
                    sc = ax.scatter(x, y, s=3, c=sum_obj, cmap=cmap, alpha=1.0/(max(gens)) * gens, marker='.', edgecolors='none', norm=norm)
                    scatter_plots.append(sc)

                    ax.set_xlim(xlims[j])
                    ax.set_ylim(ylims[i])

                    if custom_min is not None:
                        ax.axvline(custom_min[j], color="grey", lw=0.5)
                        ax.axhline(custom_min[i], color="grey", lw=0.5)

                    if not noaxes:
                        ax.text(0.95, 0.95, f"{center[0]:.3f} ± {sigma_x:.3f}\n"
                                            f"{center[1]:.3f} ± {sigma_y:.3f}",
                                verticalalignment='top', horizontalalignment='right',
                                transform=ax.transAxes, color='k', fontsize=2)

                    if i == n_cols - 1:
                        ax.set_xlabel(df_full.columns[j])
                    else:
                        ax.set_xticklabels([])
                        ax.set_xticks([])
                    if j == 0:
                        ax.set_ylabel(df_full.columns[i])
                    else:
                        ax.set_yticklabels([])
                        ax.set_yticks([])

                    if noaxes:
                        ax.set_xticks([])
                        ax.set_yticks([])

        return fig






if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='PyQt5 Application')
    parser.add_argument('filepath', type=str, help='Path to the data file (e.g., JSON file)')
    parser.add_argument('--stepgen', type=int, default=5, help='Step in generations while loading full evolution (default: 5)')

    args = parser.parse_args()

    if args.filepath:
        with open(args.filepath, 'r') as file:
            mainconfig = json.load(file)
            objective_type = mainconfig.get('objective_type', 'chi2') # chi2 = default value if key is not specified in config file
    else:
        objective_type = "chi2"

    if objective_type == "logP":
        if not hasattr(creator, "FitnessMulti"):
            creator.create("FitnessMulti", base.Fitness, weights=(1.0, 1.0, 1.0, 1.0))
    else:
        if not hasattr(creator, "FitnessMulti"):
            creator.create("FitnessMulti", base.Fitness, weights=(-1.0, -1.0, -1.0, -1.0))

    if not hasattr(creator, "Individual"):
        creator.create("Individual", list, fitness=creator.FitnessMulti)

    app = QtWidgets.QApplication(sys.argv)
    mainWindow = MyApp(args.filepath, stepgen=args.stepgen)
    mainWindow.show()
    sys.exit(app.exec_())
