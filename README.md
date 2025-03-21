# Enhanced LEEPO: Liposome Encapsulation Efficiency Predictor & Optimizer

## Overview
Enhanced LEEPO is an advanced computational tool for predicting and optimizing the encapsulation efficiency of hydrophilic drugs in liposomes. Based on the mathematical model developed by Xu et al. (2012), this software integrates multiple factors that affect encapsulation including media effects, lipid properties, preparation methods, and geometric considerations.

## Key Features
- **Media Effects Modeling**: Incorporates buffer type and ionic strength effects on lipid molecular area
- **Elliptical Geometry**: Supports both spherical and elliptical (prolate/oblate) liposome geometries
- **Lipid Database**: Built-in properties for common lipids with temperature-dependent calculations
- **Preparation Method Comparison**: Analyzes different preparation techniques for optimal encapsulation
- **Drug-Lipid Interaction Detection**: Identifies potential interactions that may affect encapsulation
- **Advanced Visualization**: Size distribution, pie charts, and comparative plots
- **Defect Optimization**: Iteratively adjusts defect percentage to match experimental values

## Installation

### Prerequisites
```bash
pip install numpy pandas matplotlib scipy
```

### Setup
1. Download the `enhanced_leepo.py` file
2. Ensure you have Python 3.6+ installed
3. Install the required dependencies

## Usage

### Basic Command
```bash
python enhanced_leepo.py mu sigma size_range defect thickness molecular_area concentration volume [expected_ee] [aspect_ratio] [lipid_composition] [buffer] [prep_method] [temperature] [drug_name]
```

### Examples

#### Basic Run (Spherical Liposomes)
```bash
python enhanced_leepo.py 190.0 38.0 1000 0.23 5.2 0.45 140 5
```

#### With Optimization Target
```bash
python enhanced_leepo.py 190.0 38.0 1000 0.23 5.2 0.45 140 5 60.04
```

#### Using Lipid Database (set thickness and mol_area to 0)
```bash
python enhanced_leepo.py 190.0 38.0 1000 0.23 0 0 140 5 60.04 1.0 "DSPC:CHOLESTEROL:DPTAP=6:3:2" "HEPES_10mM"
```

#### With Elliptical Geometry (aspect_ratio = 1.5)
```bash
python enhanced_leepo.py 190.0 38.0 1000 0.23 5.2 0.45 140 5 60.04 1.5
```

#### Full Example with All Parameters
```bash
python enhanced_leepo.py 190.0 38.0 1000 0.23 0 0 140 5 60.04 1.0 "DSPC:CHOLESTEROL:DPTAP=6:3:2" "HEPES_10mM" "THIN_FILM_HYDRATION_EXTRUSION" 25.0 "TENOFOVIR"
```

## Input Parameters

| Parameter | Description | Unit | Example |
|-----------|-------------|------|---------|
| mu | Mean particle size | nm | 190.0 |
| sigma | Standard deviation | nm | 38.0 |
| size_range | Maximum size for calculation | nm | 1000 |
| defect | Anticipated defect percentage | % | 0.23 |
| thickness | Bilayer thickness (0 to use from database) | nm | 5.2 |
| molecular_area | Lipid molecular area (0 to use from database) | nm² | 0.45 |
| concentration | Lipid molar concentration | mMol/L | 140 |
| volume | Total sample volume | ml | 5 |
| expected_ee | Expected/target encapsulation efficiency [optional] | % | 60.04 |
| aspect_ratio | Ratio for elliptical liposomes (>1: prolate, <1: oblate) [optional] | - | 1.5 |
| lipid_composition | Lipid ratio in format "LIPID1:LIPID2=RATIO1:RATIO2" [optional] | - | "DSPC:CHOLESTEROL:DPTAP=6:3:2" |
| buffer | Buffer/media name from database [optional] | - | "HEPES_10mM" |
| prep_method | Preparation method from database [optional] | - | "THIN_FILM_HYDRATION_EXTRUSION" |
| temperature | Experimental temperature [optional] | °C | 25.0 |
| drug_name | Drug from database for interaction detection [optional] | - | "TENOFOVIR" |

## Output
The program generates several outputs:
1. Terminal output with calculated parameters and encapsulation efficiency
2. Multiple plots saved to a `LEEPO-plots` folder:
   - Encapsulation efficiency pie charts
   - Size distribution
   - Comparison plots (if optimization is performed)
   - Preparation method comparison

## Databases

### Lipid Database
Includes properties for common lipids like DSPC, DPPC, DMPC, DSPG, DPPG, DMPG, DPPS, DPTAP, SOPC, and cholesterol.

### Buffer Database
Includes common buffers and their effects on lipid molecular areas:
- DI_WATER
- HEPES_10mM, HEPES_20mM
- PHOSPHATE_10mM, PHOSPHATE_50mM, PHOSPHATE_100mM
- NaCl_5mM, NaCl_20mM

### Preparation Methods
Supported methods with their characteristic parameters:
- THIN_FILM_HYDRATION
- THIN_FILM_HYDRATION_EXTRUSION
- REVERSE_PHASE_EVAPORATION
- SONICATION
- ETHANOL_INJECTION
- FREEZE_THAW

### Drug Database
Information on common drugs and their expected interactions:
- TENOFOVIR
- CALCEIN
- SUCROSE
- HEMOGLOBIN
- 5,6-CF
- IODINE
- SOD

## References
1. Xu, X., Khan, M. A., & Burgess, D. J. (2012). Predicting hydrophilic drug encapsulation inside unilamellar liposomes. International journal of pharmaceutics, 423(2), 410-418. https://doi.org/10.1016/j.ijpharm.2011.12.019
