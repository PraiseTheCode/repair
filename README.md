

<img src="https://github.com/PraiseTheCode/repair/assets/32466853/f319ba55-fd3e-4970-a00b-2e97c19fbf28" width="200" height="68">

***re:pair [ Revising Eclipsing binaries analysis : Photometry And spectroscopy Infused Recipe ]***

`re:pair` is a code for simultaneous modelling of photometric and spectroscopic time series of double-lined eclipsing binaries. The code is described by Serebriakova et al. 2025, aa53605-24. The code creates a single self-consistent model of a binary that produces both its lightcurve and spectra of components, and performs a search of optimal orbital and atmospheric parameters through Multi-Objective Optimisation, using [`DEAP`](https://github.com/DEAP/) library.  


## 🛠 Installation

This code depends on two external tools:  
- [`ellc`](https://github.com/pmaxted/ellc) for light curve modeling ([P. Maxted, 2016, A&A 591, A111](https://www.aanda.org/articles/aa/pdf/2016/07/aa28579-16.pdf))
- `SynthV` for spectrum synthesis (included as source + executables) ([V. Tsymbal, 1996](https://articles.adsabs.harvard.edu/pdf/1996ASPC..108..198T), [V. Tsymbal, 2019](https://articles.adsabs.harvard.edu/pdf/2019ASPC..518..247T))

### 1. Create environment and install `ellc` and required Python packages

Make sure `ellc` is installed **first**, as it sets important dependencies.  
Then install additional packages:

```bash
# create conda environment (use any convenient name instead of repair-env)
conda create -n repair-env python=3.10
conda activate repair-env

# install ellc first
conda install -c conda-forge ellc

# install libraries required by re:pair
conda install -c conda-forge numpy scipy matplotlib pandas seaborn astropy deap

```

### 2. Set up `SynthV`

#### Option A: Use included binaries

Precompiled executables are provided in `synthV/synthV_executable/` for both Linux and MacOS (select which one you need and remove the part "_"system"" part after "SynthV")
To test if it works:

```bash
cd synthV/synthV_executable
./SynthV
```

If it prints version info and produces spectrum output files (using `SynthV.config`), it works.

#### Option B: Compile from source

You need a Fortran compiler (`ifort` recommended):

```bash
cd synthV/synthV_source
ifort -o ../synthV_executable/SynthV *.f* -zero
```


### 3. Download atmospheric models grid (required by SynthV)

Download, for example, `LLModels.tar.gz` from  
[https://fys.kuleuven.be/ster/meetings/binary-2015/gssp-software-package](https://fys.kuleuven.be/ster/meetings/binary-2015/gssp-software-package)  
and extract it somewhere. You will later set the path to these models in the `re:pair` config file.
In principle, any grids in Kurucz's format are supported, but naming convention may need to be addressed for re:pair to correctly identify grid nods.

### 4. (Optional) Set up GUI environment for interactive post-processing

The optional graphical interface requires a separate environment with PyQt5 and related packages.

```bash
conda create -n repair-gui python=3.10
conda activate repair-gui

# Base packages
conda install -c conda-forge numpy scipy matplotlib pandas seaborn astropy deap pillow mplcursors

# GUI library
pip install pyqt5
```

> ⚠️ Make sure your current working directory includes the main `re:pair` repo,  
> or adjust your `PYTHONPATH` to access the local modules (like `repair`, `binary`, `orbit`, `star`, etc.).



## 🚀 Usage

After activating your environment and navigating to the main directory:

```bash
conda activate repair-env
cd repair/main/
python repair.py
```

This will launch the GA optimization routine using the settings from your provided `config.json`.

An example configuration is included in `tests/B5/config.json`, together with the "observed" data for this test (it is the simulated noised data of 'B5' benchmark system from Serebriakova et al. 2025)

---


### ⚙️ Configuration file: `config.json`

All settings needed for running `re:pair` are defined in a single JSON file. Key sections include:

#### 🔧 Paths
| Key            | Description |
|----------------|-------------|
| `saveto`       | Where to save outputs |
| `synthVpath`   | Path to SynthV executable |
| `convolvepath` | Path to line broadening tool | - not needed by default (default behaviour is broadening via PyAstronomy modules)
| `atmmodels`    | Folder with LLModels or other atmosphere models grid |
| `abunds`       | Folder with abundance tables |
| `path_specs`   | Observed spectra directory | - fits files for epochs info + normalised spectra 
| `path_lc`      | Light curve file path |

#### 🧬 Optimization control
| Key               | Description |
|------------------|-------------|
| `nthreads`       | Number of threads to use |
| `popsize`        | Size of population for GA (genetic algorithm) optimization |
| `ngen`           | Number of generations for GA |
| `checkpoint_interval` | Save frequency for chekpoints files - needed to run from certain generation if previous run was interrupted |

#### 🪐 Initial parameters
The parameters used to fix non-optimized parameters. The optimized parameters will be overwritten later,  but it is still useful to set some realistic parameters here because the code, before running optimization, saves a plot with model with these initial parameters, which is useful to see if the config is working as intended.
Set in `params_init` — includes stellar radii (`r1`, `r2`), temperatures (`Teff1`, `Teff2`), mass ratio `q`, inclination, semi-major axis `a`, orbital period `porb`, etc. Most parameters names repeat those of [`ellc`](https://github.com/pmaxted/ellc); most of the parameters supported by `ellc` can be passed here. Added parameters needed for spectroscopy, such as Teff, metallicity, etc.

#### 🧠 Optimization targets
| Key            | Description |
|----------------|-------------|
| `params_opt`   | List of parameters to optimize |
| `params_bounds`| Allowed ranges for all parameters | - the code will generate a Sobol grid in these ranges to use as Generation 0 of GA optimisation.

---

### 📈 Output

The results (figures, logs, population parameters etc.) will be saved in the folder defined under `"saveto"`. You can monitor progress with intermediate output, interrupt at any stage, and continue from checkpoints. Logs may contain info and warnings omitted or lost in console output.

---

## 🖼 Optional: Visualize results with `beyond_repair` GUI

If you have installed the optional GUI environment (see Step 4), you can use the PyQt-based interface to explore precomputed model results interactively.

To run it:

```bash
conda activate repair-gui
cd repair/main/
python beyond_repair.py config_B5st.json --stepgen 20
```

This will open an interactive GUI useful for parameter space exploration, estimating errors, evolution tracking, etc.
