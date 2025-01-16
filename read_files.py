import os
import pandas as pd
import re

def read_mgf(file_path):
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
            elif line.startswith("PEPMASS"):
                feature['precursor_mz'] = float(line.split('=')[1].split()[0])
            elif line.startswith("CHARGE"):
                feature['charge'] = int(line.split('=')[1].strip().replace('+', ''))
            elif line.startswith("Signal_intensity"):
                feature['signal_intensity'] = float(line.split('=')[1].strip())
            elif line[0].isdigit():
                if 'fragment_spectrum' not in feature:
                    feature['fragment_spectrum'] = []
                mz, intensity = map(float, line.split())
                feature['fragment_spectrum'].append((mz, intensity))
            elif line.startswith("END IONS"):
                list_features.append(feature)
    return list_features

def read_msp(file_path):
    list_features = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        feature = {}
        for line in lines:
            if line.startswith("Name:"):
                feature = {'fragment_spectrum': []}
            elif line.startswith("PrecursorMZ:"):
                feature['precursor_mz'] = float(line.split(':')[1].strip())
            elif line.startswith("RetentionTime:"):
                feature['retention_time'] = float(line.split(':')[1].strip())
            elif line.startswith("Signal_intensity"):
                feature['signal_intensity'] = float(line.split(':')[1].strip())
            elif line[0].isdigit():
                mz, intensity = map(float, line.split())
                feature['fragment_spectrum'].append((mz, intensity))
            elif line.startswith("Num peaks:"):
                pass
            elif line == '\n':
                list_features.append(feature)
    return list_features

def read_excel(file_path):
    # Determine the file extension
    _, file_extension = os.path.splitext(file_path)
    
    # Choose the appropriate engine based on the file extension
    if file_extension.lower() in ['.xlsx', '.xlsm']:
        engine = 'openpyxl'
    elif file_extension.lower() == '.xls':
        engine = 'xlrd'
    else:
        raise ValueError(f"Unsupported file format: {file_extension}")

    try:
        df = pd.read_excel(file_path, engine=engine)
    except Exception as e:
        print(f"Error reading file {file_path}: {str(e)}")
        return []

    list_features = []
    for index, row in df.iterrows():
        try:
            feature = {
                'Peak ID': row['Peak ID'],
                'Scan': row['Scan'],
                'retention_time': row['RT (min)'],
                'precursor_mz': row['Precursor m/z'],
                'signal_intensity': row['Height'],
                'fragment_spectrum': parse_msms_spectrum(row['MSMS spectrum'])
            }
            list_features.append(feature)
        except KeyError as e:
            print(f"Warning: Missing column {e} in row {index+2} of file {file_path}. Skipping this row.")
        except Exception as e:
            print(f"Warning: Error processing row {index+2} of file {file_path}: {e}. Skipping this row.")
    return list_features

def parse_msms_spectrum(spectrum_string):
    if pd.isna(spectrum_string):
        return []
    
    # Use regular expression to split the peaks
    peaks = re.findall(r'(\d+\.\d+)\s+(\d+)', spectrum_string)
    
    spectrum = []
    for mz, intensity in peaks:
        try:
            spectrum.append((float(mz), float(intensity)))
        except ValueError:
            print(f"Warning: Could not parse peak '{mz} {intensity}' in spectrum: {spectrum_string}")
    
    return spectrum

def collect_files(directory):
    mgf_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.mgf')]
    msp_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.msp')]
    excel_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.xlsx') or f.endswith('.xls')]
    return mgf_files, msp_files, excel_files

def read_features(directory):
    mgf_files, msp_files, excel_files = collect_files(directory)
    all_list_features = []
    for mgf_file in mgf_files:
        list_features = read_mgf(mgf_file)
        all_list_features.append((mgf_file, list_features))
    for msp_file in msp_files:
        list_features = read_msp(msp_file)
        all_list_features.append((msp_file, list_features))
    for excel_file in excel_files:
        list_features = read_excel(excel_file)
        all_list_features.append((excel_file, list_features))
    return all_list_features



