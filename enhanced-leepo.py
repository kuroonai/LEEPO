#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Enhanced LEEPO: Liposome Encapsulation Efficiency Predictor & Optimizer
Based on research by Xu et al. (2012) "Predicting hydrophilic drug encapsulation inside unilamellar liposomes"
International Journal of Pharmaceutics 423.2: 410-418.
DOI: https://doi.org/10.1016/j.ijpharm.2011.12.019

Enhancements:
- Media effects on lipid molecular area
- Elliptical geometry support
- Lipid database with bilayer thickness data
- Temperature effects
- Preparation method comparison
- Drug-lipid interaction detection
'''

import sys
import os
import time
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import constants as cons
from scipy.stats import lognorm

# Configuration
MAX_ITERATIONS = 1000
CONVERGENCE_THRESHOLD = 0.00001
REFERENCE_THRESHOLD = 0.01
DEFAULT_MAX_EE_VALUE = 74.05
DEFAULT_ASPECT_RATIO = 1.0  # Default is sphere
DEFAULT_TEMPERATURE = 25.0  # Default temperature in °C

# Database of lipids with their properties
LIPID_DATABASE = {
    # Format: 'lipid_name': {'thickness': bilayer thickness in nm at 25°C, 
    #                        'transition_temp': transition temperature in °C,
    #                        'molecular_area': area in nm² at 25°C in water}
    'DSPC': {'thickness': 5.1, 'transition_temp': 55, 'molecular_area': 0.500},
    'DPPC': {'thickness': 4.6, 'transition_temp': 41, 'molecular_area': 0.480},
    'DMPC': {'thickness': 3.8, 'transition_temp': 23, 'molecular_area': 0.460},
    'DSPG': {'thickness': 5.28, 'transition_temp': 53.6, 'molecular_area': 0.520},
    'DPPG': {'thickness': 5.02, 'transition_temp': 40.2, 'molecular_area': 0.500},
    'DMPG': {'thickness': 4.84, 'transition_temp': 22.4, 'molecular_area': 0.480},
    'DPPS': {'thickness': 5.2, 'transition_temp': None, 'molecular_area': 0.550},
    'DPTAP': {'thickness': 4.0, 'transition_temp': 45.4, 'molecular_area': 0.567},
    'SOPC': {'thickness': 4.9, 'transition_temp': -14, 'molecular_area': 0.520},
    'CHOLESTEROL': {'thickness': 2.0, 'transition_temp': None, 'molecular_area': 0.376}
}

# Buffer/media database with their effect on lipid molecular area
BUFFER_DATABASE = {
    # Format: 'buffer_name': {'pH': pH value, 'ionic_strength': in mM,
    #                         'area_modifier': multiplier for molecular area (depends on lipid charge)}
    'DI_WATER': {'pH': 7.0, 'ionic_strength': 0, 
                'area_modifier': {'neutral': 1.00, 'cationic': 1.00, 'anionic': 1.00}},
    'HEPES_10mM': {'pH': 7.4, 'ionic_strength': 10, 
                  'area_modifier': {'neutral': 1.03, 'cationic': 1.01, 'anionic': 0.97}},
    'HEPES_20mM': {'pH': 7.4, 'ionic_strength': 20, 
                  'area_modifier': {'neutral': 1.02, 'cationic': 1.00, 'anionic': 0.96}},
    'PHOSPHATE_10mM': {'pH': 7.4, 'ionic_strength': 10, 
                      'area_modifier': {'neutral': 1.00, 'cationic': 0.94, 'anionic': 0.98}},
    'PHOSPHATE_50mM': {'pH': 7.4, 'ionic_strength': 50, 
                      'area_modifier': {'neutral': 0.99, 'cationic': 0.93, 'anionic': 0.97}},
    'PHOSPHATE_100mM': {'pH': 7.4, 'ionic_strength': 100, 
                       'area_modifier': {'neutral': 0.99, 'cationic': 0.93, 'anionic': 0.97}},
    'NaCl_5mM': {'pH': 7.0, 'ionic_strength': 5, 
                'area_modifier': {'neutral': 1.01, 'cationic': 0.98, 'anionic': 0.99}},
    'NaCl_20mM': {'pH': 7.0, 'ionic_strength': 20, 
                 'area_modifier': {'neutral': 1.00, 'cationic': 0.97, 'anionic': 0.98}}
}

# Preparation method database
PREP_METHOD_DATABASE = {
    # Format: 'method_name': {'lamellarity': uni/multi, 'size_range': [min, max] in nm,
    #                         'distribution_factor': affects width of size distribution,
    #                         'defect_factor': introduces defect % in the bilayer}
    'THIN_FILM_HYDRATION': {'lamellarity': 'multi', 'size_range': [400, 5000], 
                           'distribution_factor': 0.3, 'defect_factor': 0.25},
    'THIN_FILM_HYDRATION_EXTRUSION': {'lamellarity': 'uni', 'size_range': [80, 200], 
                                     'distribution_factor': 0.2, 'defect_factor': 0.20},
    'REVERSE_PHASE_EVAPORATION': {'lamellarity': 'uni', 'size_range': [150, 450], 
                                 'distribution_factor': 0.25, 'defect_factor': 0.23},
    'SONICATION': {'lamellarity': 'uni', 'size_range': [50, 150], 
                  'distribution_factor': 0.2, 'defect_factor': 0.30},
    'ETHANOL_INJECTION': {'lamellarity': 'uni', 'size_range': [100, 300], 
                         'distribution_factor': 0.22, 'defect_factor': 0.22},
    'FREEZE_THAW': {'lamellarity': 'uni', 'size_range': [150, 500], 
                   'distribution_factor': 0.28, 'defect_factor': 0.18}
}

# Drug interaction database - used to detect potential drug-lipid interactions
DRUG_INTERACTION_DATABASE = {
    # Format: 'drug_name': {'charge': +1, 0, -1, 'logP': partition coefficient, 
    #                      'expected_deviation': expected deviation from model}
    'TENOFOVIR': {'charge': -1, 'logP': -2.3, 'expected_deviation': 0.05},
    'CALCEIN': {'charge': -4, 'logP': -5.02, 'expected_deviation': 0.1},
    'SUCROSE': {'charge': 0, 'logP': -3.7, 'expected_deviation': 0.02},
    'HEMOGLOBIN': {'charge': 0, 'logP': None, 'expected_deviation': 0.25},
    '5,6-CF': {'charge': -2, 'logP': -1.5, 'expected_deviation': 0.12},
    'IODINE': {'charge': 0, 'logP': 2.5, 'expected_deviation': 0.15},
    'SOD': {'charge': 0, 'logP': None, 'expected_deviation': 0.20}
}

def create_output_folder():
    """Create output folder for plots if it doesn't exist"""
    cwd = os.getcwd()
    
    if os.name == 'nt':
        folder = f'{cwd}\\LEEPO-plots'
    else:
        folder = f'{cwd}/LEEPO-plots'
        
    if not os.path.isdir(folder):
        os.mkdir('LEEPO-plots')
        
    os.chdir(folder)
    return folder

def validate_inputs(args):
    """Validate command line arguments"""
    usage_msg = '''
    Syntax error encountered:

    Usage:

    Enter the system parameters in the order of:

    1. Mean (μ - nm)
    2. Standard deviation (σ - nm)
    3. Size range (eg: for 1-1000 nm enter "1000")
    4. Anticipated defect (%)
    5. Bilayer thickness (d - nm) or leave as 0 to use value from lipid database
    6. Molecular area of lipid (a - nm^2) or leave as 0 to use value from database
    7. Lipid molar concentration (c - mMol/L)
    8. Total volume (V - ml)
    9. Expected Encapsulation efficiency (%) [optional]
    10. Aspect ratio (a/b) for elliptical liposomes [optional, default=1.0]
    11. Lipid composition (e.g., "DSPC:CHOLESTEROL:DPTAP=6:3:2") [optional]
    12. Buffer/media name (e.g., "HEPES_10mM") [optional, default="DI_WATER"]
    13. Preparation method (e.g., "THIN_FILM_HYDRATION_EXTRUSION") [optional]
    14. Temperature (°C) [optional, default=25.0]
    15. Drug name [optional]

    '''
    
    if len(args) < 9:
        print(usage_msg)
        raise ValueError("Insufficient number of arguments provided")
        
    try:
        # Parse basic parameters
        params = [float(arg) for arg in args[1:9]]
        
        # Convert size range to integer
        params[2] = int(params[2])
        
        # Check for negative values where inappropriate
        if any(p < 0 for p in params):
            raise ValueError("All parameters must be positive values")
        
        # Add optional parameters with defaults
        if len(args) > 9:
            params.append(float(args[9]))  # Expected EE
        else:
            params.append(None)  # No expected EE
            
        if len(args) > 10:
            params.append(float(args[10]))  # Aspect ratio
        else:
            params.append(DEFAULT_ASPECT_RATIO)  # Default aspect ratio
            
        if len(args) > 11:
            params.append(args[11])  # Lipid composition
        else:
            params.append(None)  # No specific lipid composition
            
        if len(args) > 12:
            params.append(args[12])  # Buffer/media
        else:
            params.append("DI_WATER")  # Default buffer
            
        if len(args) > 13:
            params.append(args[13])  # Preparation method
        else:
            params.append(None)  # No specific preparation method
            
        if len(args) > 14:
            params.append(float(args[14]))  # Temperature
        else:
            params.append(DEFAULT_TEMPERATURE)  # Default temperature
            
        if len(args) > 15:
            params.append(args[15])  # Drug name
        else:
            params.append(None)  # No specific drug
            
        return params
        
    except ValueError as e:
        print(usage_msg)
        raise ValueError(f"Invalid parameter value: {e}")

def get_lipid_properties(lipid_composition, temperature=25.0, buffer="DI_WATER"):
    """Calculate bilayer thickness and molecular area from lipid composition, temperature, and buffer"""
    if lipid_composition is None:
        return None, None
    
    # Parse lipid composition, e.g., "DSPC:CHOLESTEROL:DPTAP=6:3:2"
    parts = lipid_composition.split('=')
    if len(parts) != 2:
        raise ValueError(f"Invalid lipid composition format: {lipid_composition}. Use format 'LIPID1:LIPID2=RATIO1:RATIO2'")
    
    lipids = parts[0].split(':')
    ratios_str = parts[1].split(':')
    ratios = [float(r) for r in ratios_str]
    
    # Verify all lipids are in database
    for lipid in lipids:
        if lipid not in LIPID_DATABASE:
            raise ValueError(f"Lipid {lipid} not found in database")
    
    # Calculate total ratio
    total_ratio = sum(ratios)
    
    # Calculate weighted average thickness
    thickness = 0
    for lipid, ratio in zip(lipids, ratios):
        lipid_thickness = LIPID_DATABASE[lipid]['thickness']
        thickness += lipid_thickness * (ratio / total_ratio)
    
    # Calculate molecular area with buffer and temperature effects
    area = 0
    for lipid, ratio in zip(lipids, ratios):
        # Get base molecular area
        base_area = LIPID_DATABASE[lipid]['molecular_area']
        
        # Apply temperature effect
        lipid_tm = LIPID_DATABASE[lipid]['transition_temp']
        if lipid_tm is not None:
            if temperature > lipid_tm:
                # Above transition temperature, membrane becomes more fluid
                temp_factor = 1.0 + 0.01 * min(10, temperature - lipid_tm) / 10.0
            else:
                # Below transition temperature, membrane is more rigid
                temp_factor = 1.0 - 0.005 * min(10, lipid_tm - temperature) / 10.0
        else:
            temp_factor = 1.0
        
        # Determine lipid charge type for buffer effect
        if lipid == 'DPTAP':
            charge_type = 'cationic'
        elif lipid in ['DSPG', 'DPPG', 'DMPG', 'DPPS']:
            charge_type = 'anionic'
        else:
            charge_type = 'neutral'
        
        # Apply buffer effect
        buffer_factor = BUFFER_DATABASE[buffer]['area_modifier'][charge_type]
        
        # Apply weighted factors
        modified_area = base_area * temp_factor * buffer_factor
        area += modified_area * (ratio / total_ratio)
    
    return thickness, area

def get_preparation_method_params(method_name, defect):
    """Get parameters specific to preparation method"""
    if method_name is None or method_name not in PREP_METHOD_DATABASE:
        return None, defect
    
    method = PREP_METHOD_DATABASE[method_name]
    lamellarity = method['lamellarity']
    
    # Adjust defect percentage based on preparation method
    adjusted_defect = defect * method['defect_factor']
    
    return method, adjusted_defect

def calculate_probability(x, mu, sigma):
    """Calculate probability with log-normal distribution"""
    sigbymu = sigma / mu
    
    # Avoid division by zero or negative logs
    if x <= 0:
        return 0
        
    try:
        return (1 / (math.sqrt(2 * math.pi * (sigbymu)**2 * x**2))) * \
               math.exp((math.log(x) - math.log(mu))**2 / (-2 * (sigbymu)**2))
    except (ValueError, ZeroDivisionError):
        return 0

def calculate_theoretical_lognormal_distribution(mu, sigma, size_range):
    """Calculate theoretical log-normal distribution based on parameters"""
    x = np.arange(1, size_range + 1)
    
    # Convert to log-normal parameters
    s = np.sqrt(np.log(1 + (sigma/mu)**2))
    scale = mu / np.sqrt(1 + (sigma/mu)**2)
    
    # Calculate PDF
    pdf = lognorm.pdf(x, s=s, scale=scale)
    
    return x, pdf

def plot_encapsulation(entrapvolume, ee, V, expected_ee=None, is_optimized=False):
    """Create pie chart for encapsulation efficiency"""
    labels = ['Encapsulated', 'Free']
    sizes = [sum(entrapvolume), V - sum(entrapvolume)]
    colors = ['Green', 'Red']
    explode = (0.0, 0.1)   # explode 1st slice
 
    plt.figure()
    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.2f%%', shadow=True, startangle=100)
 
    plt.axis('equal')
    
    if expected_ee is not None and abs(ee - expected_ee) <= REFERENCE_THRESHOLD:
        plt.title("LEEPO-Optimized EE")
        plt.savefig("LEEPO-Optimised.pdf", bbox_inches='tight')
    else:
        plt.title("LEEPO-Predicted EE")
        plt.savefig("LEEPO-predicted.pdf", bbox_inches='tight')
    
    plt.close()

def plot_size_distribution(mu, sigma, size_range):
    """Plot log-normal size distribution"""
    x, pdf = calculate_theoretical_lognormal_distribution(mu, sigma, size_range)
    
    plt.figure()
    plt.plot(x, pdf)
    plt.title("Liposome Size Distribution")
    plt.xlabel("Vesicle size (nm)")
    plt.ylabel("Probability")
    plt.savefig("LEEPO-size-distribution.pdf", bbox_inches='tight')
    plt.close()

def plot_comparison(data1, data2, x_values, linestyle1, linestyle2, color1, color2, 
                    ylabel, title, filename):
    """Create comparison plot between two datasets"""
    plt.figure()
    plt.plot(x_values, data1, linestyle=linestyle1, color=color1, lw=1.0, label='Predicted')
    plt.plot(x_values, data2, linestyle=linestyle2, color=color2, lw=1.0, label='Optimized')
    plt.title(title)
    plt.xlabel('Vesicle size (nm)')
    plt.ylabel(ylabel)
    plt.legend(loc='best')
    plt.savefig(filename, bbox_inches='tight')
    plt.close()

def calculate_ellipsoid_volume(a, b, c):
    """Calculate volume of an ellipsoid"""
    return (4/3) * math.pi * a * b * c

def calculate_ellipsoid_surface_area(a, b, c):
    """Calculate approximate surface area of an ellipsoid using Knud Thomsen's formula"""
    # Using p=1.6075 gives a good approximation
    p = 1.6075
    return 4 * math.pi * ((a**p * b**p + a**p * c**p + b**p * c**p) / 3)**(1/p)

def calculate_vesicle_parameters(size_range, mu, sigma, defect, d, molarea, aspect_ratio):
    """Calculate vesicle parameters based on size distribution"""
    p = np.zeros(size_range)
    x = np.arange(1, size_range + 1)  # Particle sizes (nm)
    
    # Calculate probability for each size
    for i in range(size_range):
        p[i] = calculate_probability(x[i], mu, sigma)
    
    # Prevent division by zero
    p_max = max(p) if max(p) > 0 else 1
    
    # Calculate defect for each size
    dx = np.array([defect * (p_val / p_max) for p_val in p])
    
    # Calculate vesicle parameters
    v_in = np.zeros(size_range)
    area_out = np.zeros(size_range)
    area_in = np.zeros(size_range)
    lnum_out = np.zeros(size_range)
    lnum_in = np.zeros(size_range)
    k = np.zeros(size_range)
    
    for i in range(size_range):
        radius = x[i] / 2
        
        # For ellipsoid, define axes based on aspect ratio
        # For prolate: a > b = c (rugby ball shape)
        # For oblate: a = b > c (disk shape)
        if aspect_ratio > 1.0:  # Prolate
            a = radius * aspect_ratio
            b = radius
            c = radius
        elif aspect_ratio < 1.0:  # Oblate
            a = radius
            b = radius
            c = radius * aspect_ratio
        else:  # Sphere
            a = radius
            b = radius
            c = radius
        
        # Internal volume (micro L)
        # For intact liposomes (1-dx[i])
        if a > d and b > d and c > d:
            # Calculate inner ellipsoid dimensions
            a_in = max(0, a - d)
            b_in = max(0, b - d)
            c_in = max(0, c - d)
            vol_intact = calculate_ellipsoid_volume(a_in, b_in, c_in)
        else:
            vol_intact = 0
            
        # For defect liposomes (dx[i])
        if a > d/2 and b > d/2 and c > d/2:
            # Calculate inner ellipsoid dimensions for defect
            a_defect = max(0, a - d/2)
            b_defect = max(0, b - d/2)
            c_defect = max(0, c - d/2)
            vol_defect = calculate_ellipsoid_volume(a_defect, b_defect, c_defect)
        else:
            vol_defect = 0
            
        # Combine intact and defect volumes
        v_in[i] = (vol_intact * (1-dx[i]) + vol_defect * dx[i]) * 10**-18
        
        # Outer surface area
        area_out[i] = calculate_ellipsoid_surface_area(a, b, c)
        
        # Inner surface area
        if a > d and b > d and c > d:
            # For intact liposomes
            a_in = max(0, a - d)
            b_in = max(0, b - d)
            c_in = max(0, c - d)
            area_intact = calculate_ellipsoid_surface_area(a_in, b_in, c_in)
        else:
            area_intact = 0
            
        if a > d/2 and b > d/2 and c > d/2:
            # For defect liposomes
            a_defect = max(0, a - d/2)
            b_defect = max(0, b - d/2) 
            c_defect = max(0, c - d/2)
            area_defect = calculate_ellipsoid_surface_area(a_defect, b_defect, c_defect)
        else:
            area_defect = 0
            
        area_in[i] = area_intact * (1-dx[i]) + area_defect * dx[i]
        
        # Lipid numbers
        lnum_out[i] = area_out[i] / molarea
        lnum_in[i] = area_in[i] / molarea
        
        # Total lipid units weighted by probability
        k[i] = p[i] * (lnum_out[i] + lnum_in[i])
    
    return p, dx, v_in, area_out, area_in, lnum_out, lnum_in, k, x

def calculate_encapsulation(mu, sigma, size_range, defect, d, molarea, c, V, aspect_ratio):
    """Calculate encapsulation efficiency"""
    p, dx, v_in, area_out, area_in, lnum_out, lnum_in, k, x = calculate_vesicle_parameters(
        size_range, mu, sigma, defect, d, molarea, aspect_ratio)
    
    # Avoid division by zero
    k_sum = sum(k) if sum(k) > 0 else 1
    
    # Calculate vesicle number
    m = c * V * cons.N_A / 10**6 / k_sum
    
    # Calculate entrapment volumes
    vesicle_number = [m * p_val for p_val in p]
    entrap_volume = [vesicle_number[i] * v_in[i] / 1000 for i in range(size_range)]
    
    # Calculate encapsulation efficiency
    ee = sum(entrap_volume) / V * 100
    
    return ee, entrap_volume, vesicle_number, m, p, dx, v_in, area_out, area_in, lnum_out, lnum_in, k, x

def optimize_encapsulation(mu, sigma, size_range, defect_initial, d, molarea, c, V, expected_ee, aspect_ratio):
    """Optimize defect to match expected encapsulation efficiency"""
    defect = defect_initial
    ee_previous, _, _, _, _, _, _, _, _, _, _, _, _ = calculate_encapsulation(
        mu, sigma, size_range, defect, d, molarea, c, V, aspect_ratio)
    
    print('\nOptimization step begins!!\n')
    
    iteration = 0
    while iteration < MAX_ITERATIONS:
        # Adjust defect based on previous result
        if ee_previous < expected_ee:
            defect = defect + (abs(ee_previous - expected_ee) / 10)
        elif ee_previous > expected_ee:
            defect = defect - (abs(ee_previous - expected_ee) / 10)
        
        # Calculate new encapsulation efficiency
        ee, entrap_volume, vesicle_number, m, p, dx, v_in, area_out, area_in, lnum_out, lnum_in, k, x = calculate_encapsulation(
            mu, sigma, size_range, defect, d, molarea, c, V, aspect_ratio)
        
        print(f"Iteration No: {iteration} ##############################")
        print(f"\n\t\tThe encapsulation efficiency at iteration {iteration} is {ee:.4f}%, " 
              f"The defect is {abs(defect):.4f}%, Δ efficiency = {abs(ee - expected_ee):.4f}%\n")
        print(f"\t\tThe entrapment volume is {sum(entrap_volume)/(c*V)*1000:.4f} μl/μmol\n")
        
        # Check if we've reached the desired accuracy
        if abs(ee - expected_ee) <= CONVERGENCE_THRESHOLD:
            break
            
        ee_previous = ee
        iteration += 1
        
        # Avoid infinite loops if not converging
        if iteration == MAX_ITERATIONS:
            print("Warning: Maximum iteration limit reached without convergence")
    
    print('Optimization ends!!\n\n')
    print(f'The optimized values are:\n\n'
          f'Expected EE \t\t- {expected_ee:.4f}%\n\n'
          f'EE predicted by model \t- {ee_previous:.4f}%\n\n'
          f'EE Optimized by model \t- {ee:.4f}%\n\n'
          f'Defect Percentage \t- {abs(defect):.4f}%\n\n'
          f'Total Vesicles \t\t- {m:.4E}\n\n'
          f'Total internal volume \t- {sum(entrap_volume):.4f} ml\n\n'
          f'Entrapment volume \t- {sum(entrap_volume)/(c*V)*1000:.4f} μl/μmol\n\n'
          f'Number of optimization steps - {iteration} steps\n\n')
    
    return ee, entrap_volume, defect, iteration

def determine_max_ee(mu, sigma, size_range, defect_initial, d, molarea, c, V, ee_initial, aspect_ratio, max_ee_target=DEFAULT_MAX_EE_VALUE):
    """Determine maximum encapsulation efficiency range"""
    defect = defect_initial
    ee_new = ee_initial
    
    print('\nOptimal range determination step begins !!\n')
    
    iteration = 0
    while abs(ee_new - max_ee_target) >= REFERENCE_THRESHOLD and iteration < MAX_ITERATIONS:
        # Adjust defect based on previous result
        if ee_new < max_ee_target:
            defect = defect + (abs(ee_new - max_ee_target) / 10)
        elif ee_new > max_ee_target:
            defect = defect - (abs(ee_new - max_ee_target) / 10)
        
        # Calculate new encapsulation efficiency
        ee_new, entrap_volume, _, m, _, _, _, _, _, _, _, _, _ = calculate_encapsulation(
            mu, sigma, size_range, defect, d, molarea, c, V, aspect_ratio)
        
        print(f"Iteration No: {iteration} ##############################")
        print(f"\n\t\tThe encapsulation efficiency at iteration {iteration} is {ee_new:.4f}%, "
              f"The defect is {abs(defect):.4f}%, Δ efficiency = {abs(ee_new - max_ee_target):.4f}%\n")
        print(f"\t\tThe entrapment volume is {sum(entrap_volume)/(c*V)*1000:.4f} μl/μmol\n")
        
        iteration += 1
        
        # Avoid infinite loops if not converging
        if iteration == MAX_ITERATIONS:
            print("Warning: Maximum iteration limit reached without convergence")
            break
    
    print(f"\nThe optimal Encapsulation efficiency range we can expect before liposome degradation is {ee_initial:.4f} - {ee_new:.4f} %")
    print(f"\n\nThe Optimal defect/void % is between {abs(defect_initial):.4f} - {defect:.4f} %\n")
    
    return ee_new, defect

def generate_comparison_plots(mu, sigma, size_range, initial_defect, optimized_defect, d, molarea, c, V, aspect_ratio):
    """Generate comparison plots between initial and optimized parameters"""
    # Calculate parameters with initial defect
    _, _, v_in1, area_out1, area_in1, lnum_out1, lnum_in1, k1, x = calculate_vesicle_parameters(
        size_range, mu, sigma, initial_defect, d, molarea, aspect_ratio)
    
    # Calculate parameters with optimized defect
    _, dx2, v_in2, area_out2, area_in2, lnum_out2, lnum_in2, k2, _ = calculate_vesicle_parameters(
        size_range, mu, sigma, optimized_defect, d, molarea, aspect_ratio)
    
    # Generate comparison plots
    plots_data = [
        (dx2, '-', '-', 'r', 'g', 'Defect %', 'Defect Percentage', 'LEEPO-defect.pdf'),
        (v_in2, '-', '-', 'r', 'g', 'Internal volume (μl)', 'Internal volume', 'LEEPO-V_internal.pdf'),
        (area_out2, '-', '-', 'r', 'g', 'a_out (nm^2)', 'Outer Surface area', 'LEEPO-area_outer.pdf'),
        (area_in2, '-', '-', 'r', 'g', 'a_in (nm^2)', 'Inner Surface area', 'LEEPO-area_inner.pdf'),
        (lnum_out2, '-', '-', 'r', 'g', '# units', 'Outer Lipid units', 'LEEPO-lin.pdf'),
        (lnum_in2, '-', '-', 'r', 'g', '# units', 'Inner Lipid units', 'LEEPO-lout.pdf'),
        (k2, '-', '-', 'r', 'g', '# units', 'Lipid units', 'LEEPO-lipid units.pdf')
    ]
    
    # Calculate vesicle numbers and entrapment volumes for both scenarios
    k1_sum = sum(k1) if sum(k1) > 0 else 1
    k2_sum = sum(k2) if sum(k2) > 0 else 1
    
    m1 = c * V * cons.N_A / 10**6 / k1_sum
    m2 = c * V * cons.N_A / 10**6 / k2_sum
    
    # Calculate with initial defect
    p1 = [calculate_probability(x_val, mu, sigma) for x_val in x]
    vesicle_number1 = [m1 * p_val for p_val in p1]
    entrap_volume1 = [vesicle_number1[i] * v_in1[i] / 1000 for i in range(size_range)]
    
    # Calculate with optimized defect
    p2 = [calculate_probability(x_val, mu, sigma) for x_val in x]
    vesicle_number2 = [m2 * p_val for p_val in p2]
    entrap_volume2 = [vesicle_number2[i] * v_in2[i] / 1000 for i in range(size_range)]
    
    # Add vesicle number and entrapment volume plots
    plots_data.extend([
        (vesicle_number2, '-', '-', 'r', 'g', '# units', 'Vesicle number', 'LEEPO-Vesiclenumber.pdf'),
        (entrap_volume2, '-', '-', 'r', 'g', 'Entrapped volume (ml)', 'Entrapped volume', 'LEEPO-entrapvol.pdf')
    ])
    
    # Generate all comparison plots
    for i, (data2, ls1, ls2, c1, c2, ylabel, title, filename) in enumerate(plots_data):
        if i == 0:  # Defect plot
            data1 = [initial_defect * (p_val / max(p1)) if max(p1) > 0 else 0 for p_val in p1]
        elif i == 7:  # Vesicle number plot
            data1 = vesicle_number1
        elif i == 8:  # Entrapment volume plot
            data1 = entrap_volume1
        else:
            # Use corresponding data from the initial calculations
            data1 = [v_in1, area_out1, area_in1, lnum_out1, lnum_in1, k1][i-1]
        
        plot_comparison(data1, data2, x, ls1, ls2, c1, c2, ylabel, title, filename)
    
    print(f"Plots completed and can be accessed at '{os.getcwd()}'\n")

def check_drug_lipid_interaction(drug_name, predicted_ee, actual_ee):
    """Check for potential drug-lipid interactions based on deviation from model"""
    if drug_name is None or drug_name not in DRUG_INTERACTION_DATABASE:
        return False, 0
    
    drug_info = DRUG_INTERACTION_DATABASE[drug_name]
    expected_deviation = drug_info['expected_deviation']
    
    # Calculate actual deviation
    actual_deviation = abs(predicted_ee - actual_ee) / predicted_ee
    
    # Compare deviations
    if actual_deviation > expected_deviation * 1.5:
        significant_interaction = True
        print(f"\nWarning: Significant drug-lipid interaction detected for {drug_name}!")
        print(f"Predicted EE: {predicted_ee:.2f}%, Actual EE: {actual_ee:.2f}%")
        print(f"Deviation: {actual_deviation*100:.2f}% (Expected: {expected_deviation*100:.2f}%)")
        
        if drug_info['charge'] != 0:
            print(f"This may be due to charge interaction with the lipid bilayer (drug charge: {drug_info['charge']})")
        if drug_info['logP'] is not None and drug_info['logP'] > 0:
            print(f"This may be due to hydrophobic interactions (drug logP: {drug_info['logP']})")
    else:
        significant_interaction = False
    
    return significant_interaction, actual_deviation

def plot_preparation_method_comparison(mu, sigma, size_range, defect, d, molarea, c, V, aspect_ratio):
    """Compare different preparation methods for the given formulation"""
    # List of methods to compare
    methods = list(PREP_METHOD_DATABASE.keys())
    
    # Calculate EE for each method
    ee_values = []
    for method_name in methods:
        method_info, adjusted_defect = get_preparation_method_params(method_name, defect)
        
        # Adjust size parameters based on method
        method_mu = np.mean(method_info['size_range'])
        method_sigma = method_mu * method_info['distribution_factor']
        
        # Calculate EE
        ee, _, _, _, _, _, _, _, _, _, _, _, _ = calculate_encapsulation(
            method_mu, method_sigma, size_range, adjusted_defect, d, molarea, c, V, aspect_ratio)
        
        ee_values.append(ee)
    
    # Create bar chart
    plt.figure(figsize=(10, 6))
    plt.bar(methods, ee_values, color='skyblue')
    plt.title('Encapsulation Efficiency by Preparation Method')
    plt.xlabel('Preparation Method')
    plt.ylabel('Encapsulation Efficiency (%)')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('LEEPO-preparation-methods.pdf', bbox_inches='tight')
    plt.close()
    
    # Create table with results
    results = pd.DataFrame({
        'Method': methods,
        'EE (%)': [f"{ee:.2f}" for ee in ee_values]
    })
    
    print("\nEncapsulation Efficiency by Preparation Method:")
    print(results.to_string(index=False))
    
    return methods, ee_values

def main():
    try:
        # Start timing
        timer = time.time()
        
        # Parse and validate inputs
        params = validate_inputs(sys.argv)
        mu, sigma, size_range, defect, d_input, molarea_input, c, V = params[:8]
        expected_ee, aspect_ratio, lipid_composition, buffer, prep_method, temperature, drug_name = params[8:]
        
        # Create output folder
        folder = create_output_folder()
        
        # Get lipid properties from composition if provided
        calc_thickness, calc_area = None, None
        if lipid_composition is not None:
            calc_thickness, calc_area = get_lipid_properties(lipid_composition, temperature, buffer)
        
        # Use calculated or input values
        d = calc_thickness if d_input == 0 and calc_thickness is not None else d_input
        molarea = calc_area if molarea_input == 0 and calc_area is not None else molarea_input
        
        # Adjust defect based on preparation method
        method_info = None
        if prep_method is not None:
            method_info, defect = get_preparation_method_params(prep_method, defect)
        
        # Print input parameters
        print("\nInput values are as follows:\n")    
        print(f"Mean (μ) - {mu:.4f} nm\n")
        print(f"\tStandard deviation (σ) - {sigma:.4f} nm\n")
        print(f"\t\tSize range - {size_range} nm\n")
        print(f"\t\t\tDefect  - {defect:.4f} %\n")
        print(f"\t\t\t\tLipid bilayer thickness (d) - {d:.4f} nm\n")
        print(f"\t\t\t\t\tLipid molecular area (a) - {molarea:.4f} nm^2\n")
        print(f"\t\t\t\t\t\tLiposome molar concentration (c) - {c:.4f} mMol/L\n")
        print(f"\t\t\t\t\t\t\tTotal volume (V) - {V:.4f} ml\n")
        
        if expected_ee is not None:
            print(f"\t\t\t\t\t\t\t\tExpected/Experimental EE - {expected_ee:.4f} %\n")
            
        print(f"\t\t\t\t\t\t\t\t\tAspect ratio - {aspect_ratio:.4f}\n")
        if aspect_ratio > 1.0:
            print(f"\t\t\t\t\t\t\t\t\tLiposome type - Prolate ellipsoid (rugby ball shape)\n")
        elif aspect_ratio < 1.0:
            print(f"\t\t\t\t\t\t\t\t\tLiposome type - Oblate ellipsoid (disk shape)\n")
        else:
            print(f"\t\t\t\t\t\t\t\t\tLiposome type - Sphere\n")
            
        if lipid_composition is not None:
            print(f"\t\t\t\t\t\t\t\t\tLipid composition - {lipid_composition}\n")
        
        print(f"\t\t\t\t\t\t\t\t\tBuffer/media - {buffer}\n")
        
        if prep_method is not None:
            print(f"\t\t\t\t\t\t\t\t\tPreparation method - {prep_method}\n")
            
        print(f"\t\t\t\t\t\t\t\t\tTemperature - {temperature:.1f}°C\n")
        
        if drug_name is not None:
            print(f"\t\t\t\t\t\t\t\t\tDrug - {drug_name}\n")
        
        # Plot size distribution
        plot_size_distribution(mu, sigma, size_range)
        
        # Calculate encapsulation efficiency with defect
        ee, entrap_volume, _, m, _, _, _, _, _, _, _, _, _ = calculate_encapsulation(
            mu, sigma, size_range, defect, d, molarea, c, V, aspect_ratio)
        
        print(f"\n\t\tThe encapsulation efficiency for the system is {ee:.4f} %\n\n")
        print(f"\n\t\tThe entrapment volume is {sum(entrap_volume)/(c*V)*1000:.4f} μl/μmol\n")
        print(f"\n\t\tDefect Percentage \t- {defect:.4f}%\n\n"
              f"\t\tTotal Vesicles \t\t- {m:.4E}\n\n"
              f"\t\tTotal internal volume \t- {sum(entrap_volume):.4f} ml\n\n"
              f"\t\tEntrapment volume \t- {sum(entrap_volume)/(c*V)*1000:.4f} μl/μmol\n\n")
        
        # Check for potential drug-lipid interactions
        if drug_name is not None and expected_ee is not None:
            interaction_detected, deviation = check_drug_lipid_interaction(drug_name, ee, expected_ee)
        
        # Calculate encapsulation efficiency with zero defect (old model)
        ee_old, entrap_volume_old, _, m_old, _, _, _, _, _, _, _, _, _ = calculate_encapsulation(
            mu, sigma, size_range, 0, d, molarea, c, V, 1.0)  # Use spherical for old model
        
        print("\nW.R.T Existing model with no parameter for defect (spherical):\n\n")
        print(f"\n\t\tThe encapsulation efficiency for the system is {ee_old:.4f} %\n\n")
        print(f"\n\t\tThe entrapment volume is {sum(entrap_volume_old)/(c*V)*1000:.4f} μl/μmol\n")
        print(f"\n\t\tTotal Vesicles \t\t- {m_old:.4E}\n\n"
              f"\t\tTotal internal volume \t- {sum(entrap_volume_old):.4f} ml\n\n"
              f"\t\tEntrapment volume \t- {sum(entrap_volume_old)/(c*V)*1000:.4f} μl/μmol\n\n")
        
        # Create pie chart for old model
        plot_encapsulation(entrap_volume_old, ee_old, V)
        plt.title("Old Model's - Predicted EE")
        plt.savefig("old model-predicted.pdf", bbox_inches='tight')
        plt.close()
        
        # Compare preparation methods
        methods, method_ee_values = plot_preparation_method_comparison(
            mu, sigma, size_range, defect, d, molarea, c, V, aspect_ratio)
        
        # Optimize if expected_ee is provided
        if expected_ee is not None and abs(ee - expected_ee) >= REFERENCE_THRESHOLD:
            # Optimize defect to match expected_ee
            ee_opt, entrap_volume_opt, defect_opt, iterations = optimize_encapsulation(
                mu, sigma, size_range, defect, d, molarea, c, V, expected_ee, aspect_ratio)
            
            # Create pie chart for optimized model
            plot_encapsulation(entrap_volume_opt, ee_opt, V, expected_ee, True)
            
            # Determine max EE
            ee_max, defect_max = determine_max_ee(
                mu, sigma, size_range, defect_opt, d, molarea, c, V, ee_opt, aspect_ratio)
            
            # Generate comparison plots
            generate_comparison_plots(mu, sigma, size_range, defect, defect_opt, d, molarea, c, V, aspect_ratio)
        else:
            # Create pie chart for predicted model
            plot_encapsulation(entrap_volume, ee, V)
        
        # Print execution time
        print(f'\n\nThe total execution time is {time.time() - timer:.2f} seconds\n\n')
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
