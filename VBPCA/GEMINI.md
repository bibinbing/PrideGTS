# VBPCA: A MATLAB Toolbox for Variational Bayesian Principal Component Analysis

## Project Overview

This directory contains a MATLAB toolbox for performing Principal Component Analysis (PCA) using a Variational Bayesian approach, specifically designed to handle missing data. The project is authored by Alexander Ilin and Tapani Raiko.

The core of the toolbox is the `vbpca_full.m` function, which implements a PCA model with unrestricted Gaussian posterior distributions for the model parameters. The model is defined as `X(:,j) = Mu + A*S(:,j) + Noise`, where `X` is the data matrix, `Mu` is the mean, `A` is the matrix of principal components, `S` is the matrix of component scores (weights), and `Noise` is isotropic Gaussian noise.

The toolbox can handle missing data in two ways:
1.  By representing missing values as `NaN`s in a full data matrix.
2.  By using sparse matrices that only contain the observed values.

A high-level script, `RunVBPCA.m`, is provided to demonstrate a typical workflow for running the analysis and processing the results.

## Key Files

*   **`RunVBPCA.m`**: A wrapper script that demonstrates how to use the `vbpca_full` function. It takes a data matrix, the number of components to extract, and other parameters, and returns a structure containing the results, including the filtered data.
*   **`vbpca_full.m`**: The main function that implements the Variational Bayesian PCA algorithm. It offers several options for initialization, termination criteria, and algorithm choice ('vb', 'map', 'ppca').
*   **`contents.m`**: Provides a summary of all the functions included in the toolbox, categorized into PCA algorithms, cost functions, helpers, and visualization tools.
*   **`readme.txt`**: The original readme file, containing a description of the toolbox, installation instructions, and warranty information.
*   **`demopcamv.m`**: A demonstration script.

## How to Use

This is a MATLAB project. To use it, you need to have MATLAB installed.

### Basic Usage (with NaNs for missing data)

1.  Open MATLAB.
2.  Add the `VBPCA` directory to your MATLAB path.
3.  Prepare your data matrix `X` where rows represent observations (e.g., time steps) and columns represent variables (e.g., measurement stations). Missing data should be marked as `NaN`.
4.  Call the `RunVBPCA` function with your data and desired parameters.

**Example:**
```matlab
% X is your n x p data matrix with missing values as NaN
NumComponents = 5;      % Number of principal components to estimate
NumReconComp = 1;       % Which component to use for reconstruction/filtering
HypUpData = 100;        % Iterations before updating the prior
MaxIters = 1000;        % Maximum number of iterations
Display = true;         % Show progress

VBPCAResult = RunVBPCA(X, NumComponents, NumReconComp, HypUpData, MaxIters, Display);

% The data with the first principal component removed is in:
% VBPCAResult.FilteredData
```

### Advanced Usage (with sparse matrices)

For performance-critical applications with large amounts of missing data, the toolbox can use sparse matrices. This requires compiling MEX files.

1.  Ensure you have a C++ compiler installed and configured for MATLAB (e.g., MinGW on Windows).
2.  From the `VBPCA` directory in MATLAB, run the `install` script:
    ```matlab
    install
    ```
3.  For multi-threaded support, you may need to install pthreads and run `install_threads`. See `readme.txt` for details.
4.  Your data matrix `X` should be a sparse matrix containing only the observed values.

## Development Conventions

*   The code is written in MATLAB script files (`.m`).
*   Functions are well-documented with comments explaining their purpose, parameters, and return values.
*   Helper functions are provided for tasks like removing empty rows/columns (`rmempty.m`) and checking arguments (`argschk.m`).
*   The core algorithm in `vbpca_full.m` allows for flexible configuration through optional parameter/value pairs.
