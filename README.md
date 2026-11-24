# fci-phase-correction-parameter-study
MATLAB code and results for a 3,600-run simulation study of iterative phase-error correction in field-cycling MRI. Evaluates 20 threshold configurations across six phase-error amplitudes with 30 repeats per condition. Includes run-level and iteration-level results, dependencies, and documentation.

Below is the revised README text, with the updated output filenames and an explanation of the `20x6x30` label. The driver script remains named `MATLAB_30Repeat_Stochastic_Study.m`.

````markdown
## Execution instructions

### 1. Purpose

This repository contains the MATLAB driver and numerical outputs for a
3,600-run repeated simulation study of iterative field-cycling MRI
(FCI) phase-error correction.

There are two ways to use these materials:

- Inspect the archived results without rerunning the simulations.
- Regenerate the simulations using MATLAB, the required reference
  image, and the correction functions.

The output label `20x6x30` represents:

- 20 threshold configurations: four background thresholds × five
  phase-convergence thresholds.
- Six nominal phase-error amplitudes.
- Thirty stochastic repeats per condition.

Total: 20 × 6 × 30 = 3,600 simulation runs.

### 2. Software requirements

The original experiments were performed using:

- MATLAB: [MATLAB R2025b Update 5 Version: 25.2.0.3177638]
- Required toolboxes: [Optimization Toolbox: Version 25.2 (R2025b); Image Processing Toolbox: Version 25.2 (R2025b); Statistics and Machine Learning Toolbox:	Version 25.2 (R2025b)].
- Operating system: [Microsoft Windows 11 Pro Version 10.0 (Build 26300)].

The workflow uses functions including `fmincon`, `immse`, `psnr`, `ssim`, and `corr`.

### 3. Download the repository

Download and extract the complete repository, or clone it using Git.

Keep the `primary_repeated_study` directory structure intact.
Do not download only the main MATLAB script: the dependency functions and reference image are also required to regenerate
the simulations.

The driver script is:

```text
MATLAB_30Repeat_Stochastic_Study.m
```

The numerical output filenames are:

```text
ParameterSweep_20x6x30_RunResults.csv
ParameterSweep_20x6x30_IterationResults.csv
ParameterSweep_20x6x30.mat
```

### 4. Preserve the archived results

Before running the simulation, make a working copy of the complete
`primary_repeated_study` folder.

For example, name the copy:

```text
primary_repeated_study_local_rerun
```

Run the code in that copy, not in the preserved archive.

The driver saves results using fixed filenames. Running it again
can overwrite numerical outputs and figures with the same names
in its output location.

Retain the original numerical files and document any renaming or
code changes in `documentation/provenance.md`.

### 5. Provide the reference input

The driver loads:

```text
PhantomRef.mat
```

Place this input in the working copy of `primary_repeated_study`,
alongside `MATLAB_30Repeat_Stochastic_Study.m`.

The file must contain the variable `phantomRef` in the structure
expected by the driver.

If `PhantomRef.mat` is unavailable, the existing CSV results can
still be inspected, but the simulations cannot be regenerated.

### 6. Set the MATLAB working folder and function path

Open MATLAB and use the Current Folder browser to navigate to
the working copy of `primary_repeated_study`.

Then enter:

```matlab
addpath(fullfile(pwd, 'dependencies'));
```

This makes functions directly inside `dependencies` available for
the current MATLAB session.

If `dependencies` contains required subfolders, use this instead:

```matlab
addpath(genpath(fullfile(pwd, 'dependencies')));
```

The `documentation` folder does not need to be added to the MATLAB
function path.

### 7. Check required files and functions

Before starting the experiment, enter:

```matlab
assert(isfile('MATLAB_30Repeat_Stochastic_Study.m'), ...
    'The main MATLAB driver was not found.');

assert(isfolder('dependencies'), ...
    'The dependencies folder was not found.');

assert(isfile('PhantomRef.mat'), ...
    'PhantomRef.mat must be supplied before running the simulation.');

requiredFunctions = {
    'iterative_images_correction'
    'linear_correction'
    'optim_image'
    'fft2c'
    'ifft2c'
    'fmincon'
    'immse'
    'psnr'
    'ssim'
    'corr'
};

for k = 1:numel(requiredFunctions)
    assert(~isempty(which(requiredFunctions{k})), ...
        'Required function not found: %s', requiredFunctions{k});
end
```

These checks confirm that the named functions can be located.
They do not replace checking the full dependency list.

To inspect which correction functions MATLAB will use, enter:

```matlab
which iterative_images_correction -all
which linear_correction -all
which optim_image -all
which fft2c -all
which ifft2c -all
```

Confirm that the intended project functions are selected rather than similarly named functions from another installation.

### 8. Run the simulation

From the working copy of `primary_repeated_study`, enter:

```matlab
MATLAB_30Repeat_Stochastic_Study
```

Alternatively, open the script in the MATLAB Editor and select Run.

The configured design comprises:

- Background thresholds: 0.05, 0.10, 0.15, and 0.20.
- Phase-convergence thresholds: 0.025, 0.05, 0.10, 0.15, and 0.20.
- Nominal phase-error amplitudes: 2°, 10°, 20°, 40°, 60°, and 80°.
- Thirty stochastic repeats per condition.

The four background thresholds and five phase thresholds give
20 threshold configurations. Combining these with six amplitudes
gives 120 experimental conditions.

Repeating each condition 30 times produces 3,600 simulation runs.

Do not change these settings if the purpose is to reproduce the
reported experimental design.

### 9. Capture a new execution log

To record a verification run, start a diary before executing
the driver:

```matlab
diary('verification_execution_log.txt');
diary on

version
ver
computer

MATLAB_30Repeat_Stochastic_Study

diary off
```

Use a fresh log filename for each verification run.

This is a log of the new execution, not the historical log of
the original experiments. If execution stops with an error,
enter `diary off` manually after recording the error.

### 10. Expected outputs

The revised driver saves:

```text
ParameterSweep_20x6x30_RunResults.csv
ParameterSweep_20x6x30_IterationResults.csv
ParameterSweep_20x6x30.mat
```

It also generates and saves figures.

The filename label `20x6x30` describes the simulation design,
not the number of rows in the iteration-level table:

- 20 threshold configurations.
- Six nominal phase-error amplitudes.
- Thirty repeats per condition.

The archived run-level CSV contains 3,600 rows and 41 variables.

The archived iteration-level CSV contains 8,095 rows and
19 variables. Each row represents an iteration within a simulation
run; the number of iterations can differ between runs.

A new execution should be checked against the archived results.
Execution time may differ. Any differences in numerical outcomes
or iteration counts should be investigated and documented.

### 11. Inspect the results without rerunning

If the archived numerical outputs have been saved or renamed
using the updated filenames, open them in MATLAB with:

```matlab
R = readtable('ParameterSweep_20x6x30_RunResults.csv');
I = readtable('ParameterSweep_20x6x30_IterationResults.csv');

size(R)
size(I)

head(R)
head(I)
```

To inspect the MATLAB-native archive:

```matlab
S = load('ParameterSweep_20x6x30.mat');

fieldnames(S)
```

The CSV files can also be opened in spreadsheet software.

If the archived files still use the earlier filenames, either
use those filenames when loading them or rename distribution
copies and document the change. Rerunning the simulations is
not necessary solely to change filenames.

### 12. Interpretation notes

Each row in the run-level table represents one simulation run.

Each row in the iteration-level table represents one iteration
within a run. Iteration rows must not be counted as independent
simulation experiments.

Refer to `documentation/data_dictionary.csv` for variable
definitions and `documentation/provenance.md` for the origin
and processing history of the files.

Phase-recovery measures remain exploratory pending verification
of phase-sign, wrapping, and affine-component conventions.

Numerical convergence indicates that the implemented stopping
criterion was satisfied; it does not necessarily indicate
improved image quality.
````
