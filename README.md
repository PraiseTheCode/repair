

<img src="https://github.com/PraiseTheCode/repair/assets/32466853/f319ba55-fd3e-4970-a00b-2e97c19fbf28" width="200" height="68">

re:pair [ Revising Eclipsing binaries analysis : Photometry And spectroscopy Infused Recipe ]


## 🛠 Installation

This code depends on two external tools:  
- [`ellc`](https://github.com/pmaxted/ellc) for light curve modeling  
- `SynthV` for spectrum synthesis (included as source + binaries)

### 1. Create environment and install `ellc`
Make sure `ellc` is installed **first**, as it sets key environment variables.

```bash
conda create -n repair-env python=3.10
conda activate repair-env
conda install -c conda-forge ellc
```

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
