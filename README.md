

<img src="https://github.com/PraiseTheCode/repair/assets/32466853/f319ba55-fd3e-4970-a00b-2e97c19fbf28" width="200" height="68">

re:pair [ Revising Eclipsing binaries analysis : Photometry And spectroscopy Infused Recipe ]

re:pair is a code for simultaneous modelling of photometric and spectroscopic time series of double-lined eclipsing binaries. The code is described by Serebriakova et al. 2025, aa53605-24. The code creates a single self-consistent model of a binary that produces both its lightcurve and spectra of components, and performs a search of optimal orbital and atmospheric parameters through Multi-Objective Optimisation, using [DEAP](https://github.com/DEAP/) library.  


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

Precompiled executables are provided in `synthV/synthV_executable/` for both Linux and MacOS (select which one you need and remove the part "_*system*" after "SynthV"
To test if they work:

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
and extract it somewhere. You will later set the path to these models in the `repair` config file.
In principle, any grids in Kurucz's format are supported, but naming convention may need to be addressed for re:pair to correctly identify grid nods.

