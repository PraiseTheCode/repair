

<img src="https://github.com/PraiseTheCode/repair/assets/32466853/f319ba55-fd3e-4970-a00b-2e97c19fbf28" width="200" height="68">

re:pair [ Revising Eclipsing binaries analysis : Photometry And spectroscopy Infused Recipe ]


## 🛠 Installation

This code depends on two external tools:  
- [`ellc`](https://github.com/pmaxted/ellc) for light curve modeling ([P. Maxted, 2016, A&A 591, A111](https://www.aanda.org/articles/aa/pdf/2016/07/aa28579-16.pdf)
- `SynthV` for spectrum synthesis (included as source + binaries) ([V. Tsymbal, 1996](https://articles.adsabs.harvard.edu/pdf/1996ASPC..108..198T), [V. Tsymbal, 2019](https://articles.adsabs.harvard.edu/pdf/2019ASPC..518..247T) )

### 1. Create environment and install `ellc` and required Python packages

Make sure `ellc` is installed **first**, as it sets important environment variables.  
Then install additional packages needed for spectrum and light curve analysis:

```bash
conda create -n repair-env python=3.10
conda activate repair-env

# install ellc first
conda install -c conda-forge ellc

# install core analysis libraries
pip install numpy scipy matplotlib pandas seaborn astropy deap

### 2. Set up `SynthV`

#### Option A: Use included binaries (Linux/macOS only)

Precompiled executables are provided in `synthV/synthV_executable/` and `synthV/convolve_executable/`.  
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

Similarly, for the convolution tool:

```bash
cd ../convolve_source
make
mv convolve ../convolve_executable/
```

### 3. Download atmospheric models (required by SynthV)

Download `LLModels.tar.gz` from  
[https://fys.kuleuven.be/ster/meetings/binary-2015/gssp-software-package](https://fys.kuleuven.be/ster/meetings/binary-2015/gssp-software-package)  
and extract it somewhere. You will later set the path to these models in the `repair` config file.
