# PrideGTS

![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-orange)
![Platform](https://img.shields.io/badge/platform-Windows%2010%2F11-blue)
![License](https://img.shields.io/badge/license-GPL--3.0-green)

PrideGTS (Pride GNSS Time Series Analysis Software) is open-source MATLAB software for GNSS coordinate time series analysis. In one graphical interface it covers the processing from multi-source coordinate time series data to environmental loading correction, noise modeling, and time series analysis.

## Features

- Reads PBO (`.pos`), NGL (`.tenv3`) and CMONOC (`.txt`) coordinate files and converts them to Hector-compatible `.mom` files
- Preprocessing: outlier removal, known-offset import, three-component offset detection (MGLR) and missing data reconstruction (KKF, PCA iteration, RegEM, VBPCA)
- Environmental loading correction (NTAL, NTOL, HYDL) and common-mode error filtering (PCA, VBPCA, VBRPCA, FastICA, VBICA)
- Noise-model selection and parameter estimation with Hector
- Batch processing in every module

## Requirements

- MATLAB R2022b or later; base MATLAB only, no additional toolboxes
- Windows 10/11 (the Hector 2.1 executables shipped in `Hector/exe` are Windows programs)

## Installation

1. Download or clone this repository to a local directory with write permission.
2. In MATLAB, add the PrideGTS folder with all subfolders to the path (**Home → Set Path → Add with Subfolders**).
3. Double-click `Install/PrideGTS.mlappinstall`.
4. Launch PrideGTS from the **APPS** tab.

## Documentation

See the [user manual](doc/PrideGTS_Software_User_Manual_EN.pdf). It includes a step-by-step example with the sample data in `sample_data/`.

## Third-party components

PrideGTS incorporates the following third-party source code. The repository keeps the original copyright notices in the corresponding source files whenever available. Components marked as "See original distribution" are included from their original upstream distributions; their full license terms should be checked in the original upstream project before redistribution or commercial reuse.

| Component | Reference | License | Notes |
|---|---|---|---|
| Hector 2.1 | Bos et al. (2013) | GPL-3.0-or-later | Windows executables; source at https://teromovigo.com/product/hector |
| FastICA 2.5 | Hyvärinen (1999) | GPL-2.0-or-later | Copyright notice kept in source |
| VBRPCA | Babacan et al. (2012) | GPL-3.0-or-later | GPL notice kept in source header |
| VBPCA | Ilin and Raiko (2010) | No license stated | Copyright held by the original authors |
| vbICA 1.0 and Netlab | Choudrey (2003); Nabney (2002) | Netlab: BSD-3-Clause | --- |
| RegEM | Schneider (2001) | No license stated | Default parameters modified |
| GMIS (KKF) | Liu et al. (2018) | No license stated | https://doi.org/10.1007/s10291-017-0689-3 |


## License

[PrideGTS is distributed under the GNU General Public License v3.0](LICENSE)

## Contact

- Xianpeng Huang, Guilin University of Technology (gnsshxp20@glut.edu.cn or 2640355140@qq.com)
- Chaoming Liao, Wuhan University (liaochaoming@nnnu.edu.cn or 164981263@qq.com)
