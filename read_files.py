"""
Module for reading mass spectrometry feature data from various file formats.

This module provides functions to read mass features from different file formats
including Excel, TSV, MSP, and MGF files. It handles the extraction of key
information such as m/z values, retention times, and peak intensities.

Main functions/classes:
    - read_mgf: Reads features from MGF format files
    - read_msp: Reads features from MSP format files  
    - read_excel: Reads features from Excel files with specific column mapping
    - read_features: Main function that auto-detects format and reads features
    - collect_files: Collects all compatible files from a directory

Inputs:
    - Excel files with columns: Peak ID, RT [min], Precursor m/z, Height, MSMS spectrum (optional)
    - MGF/MSP files with standard mass spectrometry formats
    - TSV files with tab-separated values

Outputs:
    - List of dictionaries containing feature information (peak_id, mz, rt, intensity, etc.)
    - File metadata and dataset identification

Important arguments:
    - file_path: Path to the input file
    - input_dir: Directory containing multiple input files
"""
import os
import pandas as pd
import re
import warnings
from typing import List, Dict, Tuple, Any, Optional
import glob
import numpy as np

# Configuration for MS/MS array processing
MAX_MZ = 2000  # Maximum m/z value for array size
MIN_INTENSITY = 0.0  # Minimum intensity threshold


def parse_msms_string_to_arrays(msms_string: str, 
                               max_mz: int = MAX_MZ,
                               min_intensity: float = MIN_INTENSITY) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert MS/MS spectrum string to boolean and intensity arrays.
    
    IMPORTANT: Uses round() for nominal m/z conversion. This affects which integer bin
    each peak falls into (e.g., 149.7 → 150, not 149). Consider impact on results.
    
    Inputs:
        msms_string (str): MS/MS spectrum in format "149.02344 2466;150.02811 260;151.03277 0"
        max_mz (int): Maximum m/z value for array size
        min_intensity (float): Minimum intensity threshold
        
    Outputs:
        Tuple[np.ndarray, np.ndarray]: 
            - Boolean array where index = nominal m/z, value = peak present
            - Intensity array where index = nominal m/z, value = intensity
    """
    # Initialize arrays
    peak_present = np.zeros(max_mz, dtype=bool)
    intensities = np.zeros(max_mz, dtype=np.float32)
    
    if not msms_string or msms_string.strip() == "":
        return peak_present, intensities
    
    try:
        # Split by semicolon to get individual peaks
        peaks = msms_string.strip().split(';')
        
        for peak in peaks:
            if peak.strip():
                # Split by space to get m/z and intensity
                parts = peak.strip().split()
                if len(parts) >= 2:
                    mz = float(parts[0])
                    intensity = float(parts[1])
                    
                    # Convert to nominal m/z (integer) using rounding
                    # NOTE: Using round() instead of truncation - this impacts results!
                    # Examples: 149.7 → 150, 149.4 → 149, 149.5 → 150
                    # Alternative approaches: int(mz) for truncation, math.floor(mz) for floor
                    nominal_mz = int(round(mz))
                    
                    # Check bounds and intensity threshold
                    if 0 <= nominal_mz < max_mz and intensity >= min_intensity:
                        peak_present[nominal_mz] = True
                        # If multiple peaks at same nominal m/z, take maximum intensity
                        intensities[nominal_mz] = max(intensities[nominal_mz], intensity)
                        
    except Exception as e:
        warnings.warn(f"Error parsing MS/MS string '{msms_string}': {e}")
    
    return peak_present, intensities


def add_msms_arrays_to_feature(feature: Dict[str, Any], 
                              max_mz: int = MAX_MZ,
                              min_intensity: float = MIN_INTENSITY) -> Dict[str, Any]:
    """
    Add MS/MS boolean and intensity arrays to a feature dictionary.
    
    Inputs:
        feature (Dict[str, Any]): Feature dictionary from file reading
        max_mz (int): Maximum m/z value for array size
        min_intensity (float): Minimum intensity threshold
        
    Outputs:
        Dict[str, Any]: Feature dictionary with added 'msms_peaks' and 'msms_intensities'
    """
    # Extract MS/MS string from feature
    msms_string = feature.get('ms2', '') or feature.get('msms', '') or feature.get('MSMS spectrum', '')
    
    # Convert to arrays
    peak_present, intensities = parse_msms_string_to_arrays(msms_string, max_mz, min_intensity)
    
    # Add arrays to feature
    feature['msms_peaks'] = peak_present
    feature['msms_intensities'] = intensities
    feature['has_msms'] = np.any(peak_present)  # Quick check for MS/MS availability
    
    return feature


def read_mgf(file_path: str) -> List[Dict[str, Any]]:
    """
    Read mass features from an MGF (Mascot Generic Format) file and preprocess MS/MS data into arrays.
    
    Inputs:
        file_path (str): Path to the MGF file to read
    
    Outputs:
        List[Dict[str, Any]]: List of feature dictionaries containing:
            - title: Feature identifier
            - retention_time: RT in seconds  
            - precursor_mz: Precursor m/z value
            - charge: Ion charge state
            - signal_intensity: Peak intensity
            - fragment_spectrum: List of (mz, intensity) tuples for MS/MS
            - msms_peaks: Boolean array for fast similarity calculation
            - msms_intensities: Intensity array for fast similarity calculation
    """
    list_features = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        feature = {}
        for line in lines:
            if line.startswith("BEGIN IONS"):
                feature = {}
            elif line.startswith("TITLE"):
                feature['title'] = line.split('=')[1].strip()
            elif line.startswith("RTINSECONDS"):
                feature['retention_time'] = float(line.split('=')[1].strip())
                # Convert to minutes for consistency with other formats
                feature['rt'] = feature['retention_time'] / 60.0
            elif line.startswith("PEPMASS"):
                feature['precursor_mz'] = float(line.split('=')[1].split()[0])
                feature['mz'] = feature['precursor_mz']  # Standardize field name
            elif line.startswith("CHARGE"):
                feature['charge'] = int(line.split('=')[1].strip().replace('+', ''))
            elif line.startswith("Signal_intensity"):
                feature['signal_intensity'] = float(line.split('=')[1].strip())
                feature['intensity'] = feature['signal_intensity']  # Standardize field name
            elif line[0].isdigit():
                if 'fragment_spectrum' not in feature:
                    feature['fragment_spectrum'] = []
                mz, intensity = map(float, line.split())
                feature['fragment_spectrum'].append((mz, intensity))
            elif line.startswith("END IONS"):
                # Convert fragment spectrum to MS/MS string format for array processing
                if 'fragment_spectrum' in feature and feature['fragment_spectrum']:
                    msms_string = ';'.join([f"{mz} {intensity}" for mz, intensity in feature['fragment_spectrum']])
                    feature['ms2'] = msms_string
                else:
                    feature['ms2'] = ''
                
                # Add MS/MS arrays to feature
                feature = add_msms_arrays_to_feature(feature)
                
                list_features.append(feature)
    return list_features

def read_msp(file_path: str) -> List[Dict[str, Any]]:
    """
    Read mass features from an MSP format file and preprocess MS/MS data into arrays.
    
    Inputs:
        file_path (str): Path to the MSP file to read
    
    Outputs:
        List[Dict[str, Any]]: List of feature dictionaries with MS/MS arrays
    """
    list_features = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        feature = {}
        for line in lines:
            if line.startswith("Name:"):
                feature = {'fragment_spectrum': []}
            elif line.startswith("PrecursorMZ:"):
                feature['precursor_mz'] = float(line.split(':')[1].strip())
                feature['mz'] = feature['precursor_mz']  # Standardize field name
            elif line.startswith("RetentionTime:"):
                feature['retention_time'] = float(line.split(':')[1].strip())
                feature['rt'] = feature['retention_time']  # Standardize field name
            elif line.startswith("Signal_intensity"):
                feature['signal_intensity'] = float(line.split(':')[1].strip())
                feature['intensity'] = feature['signal_intensity']  # Standardize field name
            elif line[0].isdigit():
                mz, intensity = map(float, line.split())
                feature['fragment_spectrum'].append((mz, intensity))
            elif line.startswith("Num peaks:"):
                pass
            elif line == '\n':
                # Convert fragment spectrum to MS/MS string format for array processing
                if 'fragment_spectrum' in feature and feature['fragment_spectrum']:
                    msms_string = ';'.join([f"{mz} {intensity}" for mz, intensity in feature['fragment_spectrum']])
                    feature['ms2'] = msms_string
                else:
                    feature['ms2'] = ''
                
                # Add MS/MS arrays to feature
                feature = add_msms_arrays_to_feature(feature)
                
                list_features.append(feature)
    return list_features

def read_excel(file_path: str) -> List[Dict[str, Any]]:
    """
    Read features from an Excel file and preprocess MS/MS data into arrays.
    
    Parameters:
    -----------
    file_path : str
        Path to the Excel file
        
    Returns:
    --------
    list
        List of feature dictionaries with MS/MS arrays
    """
    print(f"Reading features from {file_path}...")
    
    # Determine the engine based on file extension
    _, ext = os.path.splitext(file_path)
    if ext.lower() in ['.xlsx', '.xlsm']:
        engine = 'openpyxl'
    else:  # .xls
        engine = 'xlrd'
    
    try:
        # Read the Excel file
        df = pd.read_excel(file_path, engine=engine)
        
        # Check required columns
        required_columns = ['Peak ID', 'Scan', 'RT (min)', 'Precursor m/z', 'Height']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            warnings.warn(f"Missing columns in {file_path}: {missing_columns}")
        
        # Extract features
        features = []
        for _, row in df.iterrows():
            try:
                feature = {
                    'peak_id': row.get('Peak ID', ''),
                    'scan': row.get('Scan', 0),
                    'rt': row.get('RT (min)', 0.0),
                    'mz': row.get('Precursor m/z', 0.0),
                    'intensity': row.get('Height', 0.0),
                    'ms2': row.get('MSMS spectrum', '')
                }
                
                # Add MS/MS arrays to feature
                feature = add_msms_arrays_to_feature(feature)
                
                features.append(feature)
            except Exception as e:
                warnings.warn(f"Error processing row in {file_path}: {e}")
        
        print(f"Extracted {len(features)} features from {file_path}")
        return features
    
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return []

def collect_files(directory: str, file_extension: str = ".xlsx") -> List[str]:
    """
    Collect all files with the specified extension from the directory.
    
    Parameters:
    -----------
    directory : str
        Directory to search for files
    file_extension : str
        File extension to filter by (default: .xlsx)
        
    Returns:
    --------
    list
        List of file paths
    """
    print(f"Collecting {file_extension} files from {directory}...")
    
    # Check if directory exists
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")
    
    # Collect files
    file_pattern = os.path.join(directory, f"*{file_extension}")
    files = glob.glob(file_pattern)
    
    # Also check for .xls files if looking for Excel files
    if file_extension.lower() == ".xlsx":
        xls_files = glob.glob(os.path.join(directory, "*.xls"))
        files.extend(xls_files)
    
    print(f"Found {len(files)} files")
    return files

def read_features(file_path: str) -> List[Dict[str, Any]]:
    """
    Read features from a file based on its extension.
    
    Parameters:
    -----------
    file_path : str
        Path to the file
        
    Returns:
    --------
    list
        List of feature dictionaries
    """
    _, ext = os.path.splitext(file_path)
    
    if ext.lower() in ['.xlsx', '.xls', '.xlsm']:
        return read_excel(file_path)
    else:
        raise ValueError(f"Unsupported file extension: {ext}")

def test_read_excel(directory: str) -> None:
    """
    Test reading Excel files from a directory.
    
    Parameters:
    -----------
    directory : str
        Directory containing Excel files
    """
    # Collect Excel files
    excel_files = collect_files(directory, ".xlsx")
    
    # Read features from each file
    all_features = []
    for excel_file in excel_files:
        list_features = read_excel(excel_file)
        all_features.append((excel_file, list_features))
    
    # Print summary
    print(f"Read {len(all_features)} Excel files")
    for file_path, features in all_features:
        print(f"  {file_path}: {len(features)} features")



