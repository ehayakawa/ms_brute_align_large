"""
Module for fast MS/MS spectral similarity calculations.

This module provides ultra-fast cosine similarity calculations between MS/MS spectra
using precomputed boolean and intensity arrays. It implements the matching logic
defined in the design document for features with MS/MS data.

Main functions/classes:
    - fast_cosine_similarity: Core cosine similarity using precomputed arrays
    - calculate_spectral_similarity: Main interface for feature-to-feature comparison
    - has_msms_data: Quick check for MS/MS availability
    - get_msms_stats: Get statistics about MS/MS data in a feature

Inputs:
    - Feature dictionaries with precomputed 'msms_peaks' and 'msms_intensities' arrays
    - Configuration parameters for minimum shared peaks and similarity thresholds

Outputs:
    - Cosine similarity scores (0.0 to 1.0)
    - Boolean flags for MS/MS data availability
    - MS/MS statistics for analysis

Important arguments:
    - min_shared_peaks: Minimum number of shared peaks required (default: 3)
    - cosine_threshold: Minimum cosine similarity for valid matches (default: 0.0)
"""
import numpy as np
import logging
from typing import Dict, Any, Tuple, Optional

# Configure logger for this module
logger = logging.getLogger(__name__)

# Default configuration parameters
DEFAULT_MIN_SHARED_PEAKS = 3
DEFAULT_COSINE_THRESHOLD = 0.0


def has_msms_data(feature: Dict[str, Any]) -> bool:
    """
    Quick check if a feature has usable MS/MS data using precomputed flag.
    
    Inputs:
        feature (Dict[str, Any]): Feature dictionary with processed MS/MS data
        
    Outputs:
        bool: True if feature has MS/MS data, False otherwise
    """
    return feature.get('has_msms', False)


def fast_cosine_similarity(peaks1: np.ndarray, 
                          intensities1: np.ndarray,
                          peaks2: np.ndarray, 
                          intensities2: np.ndarray,
                          min_shared_peaks: int = DEFAULT_MIN_SHARED_PEAKS) -> Tuple[float, int]:
    """
    Calculate cosine similarity between two preprocessed spectra using vectorized operations.
    
    This is the core function that performs ultra-fast cosine similarity calculation
    using boolean masks and numpy vectorization.
    
    Inputs:
        peaks1 (np.ndarray): Boolean array for spectrum 1 (peak presence)
        intensities1 (np.ndarray): Intensity array for spectrum 1
        peaks2 (np.ndarray): Boolean array for spectrum 2 (peak presence)  
        intensities2 (np.ndarray): Intensity array for spectrum 2
        min_shared_peaks (int): Minimum number of shared peaks required
        
    Outputs:
        Tuple[float, int]: (cosine_similarity_score, number_of_shared_peaks)
    """
    # Find shared peaks using boolean AND operation (extremely fast)
    shared_peaks_mask = peaks1 & peaks2
    num_shared = np.sum(shared_peaks_mask)
    
    # Check minimum shared peaks requirement
    if num_shared < min_shared_peaks:
        return 0.0, num_shared
    
    # Extract intensity vectors only where both spectra have peaks
    # This is much faster than iterating through all m/z values
    vec1 = intensities1[shared_peaks_mask]
    vec2 = intensities2[shared_peaks_mask]
    
    # Calculate cosine similarity using vectorized operations
    dot_product = np.dot(vec1, vec2)
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)
    
    # Handle zero norm cases
    if norm1 == 0.0 or norm2 == 0.0:
        return 0.0, num_shared
    
    similarity = dot_product / (norm1 * norm2)
    
    # Ensure result is between 0 and 1
    similarity = max(0.0, min(1.0, similarity))
    
    return similarity, num_shared


def calculate_spectral_similarity(feature1: Dict[str, Any], 
                                 feature2: Dict[str, Any],
                                 min_shared_peaks: int = DEFAULT_MIN_SHARED_PEAKS,
                                 cosine_threshold: float = DEFAULT_COSINE_THRESHOLD) -> float:
    """
    Calculate spectral similarity between two features with preprocessed MS/MS data.
    
    This is the main interface function that should be used by graph construction
    and other modules for MS/MS similarity calculation.
    
    Inputs:
        feature1 (Dict[str, Any]): First feature with 'msms_peaks' and 'msms_intensities'
        feature2 (Dict[str, Any]): Second feature with 'msms_peaks' and 'msms_intensities'
        min_shared_peaks (int): Minimum number of shared peaks required
        cosine_threshold (float): Minimum cosine similarity for valid match
        
    Outputs:
        float: Cosine similarity score (0.0 to 1.0), 0.0 if below threshold or no MS/MS
    """
    # Quick check: both features must have MS/MS data
    if not (has_msms_data(feature1) and has_msms_data(feature2)):
        return 0.0
    
    # Get precomputed arrays
    peaks1 = feature1['msms_peaks']
    intensities1 = feature1['msms_intensities']
    peaks2 = feature2['msms_peaks']
    intensities2 = feature2['msms_intensities']
    
    # Calculate fast cosine similarity
    similarity, num_shared = fast_cosine_similarity(
        peaks1, intensities1, peaks2, intensities2, min_shared_peaks
    )
    
    # Apply threshold filter
    if similarity < cosine_threshold:
        logger.debug(f"Spectral similarity {similarity:.3f} below threshold {cosine_threshold}")
        return 0.0
    
    logger.debug(f"Spectral similarity: {similarity:.3f} "
                f"(shared peaks: {num_shared})")
    
    return similarity


def get_msms_stats(feature: Dict[str, Any]) -> Dict[str, Any]:
    """
    Get statistics about MS/MS data for a feature.
    
    Inputs:
        feature (Dict[str, Any]): Feature dictionary with processed MS/MS data
        
    Outputs:
        Dict[str, Any]: Dictionary with MS/MS statistics
    """
    if not has_msms_data(feature):
        return {
            'has_msms': False,
            'num_peaks': 0, 
            'max_mz': 0,
            'total_intensity': 0.0,
            'max_intensity': 0.0
        }
    
    peaks = feature['msms_peaks']
    intensities = feature['msms_intensities']
    
    # Get peak positions where peaks exist
    peak_positions = np.where(peaks)[0]
    peak_intensities = intensities[peaks]
    
    return {
        'has_msms': True,
        'num_peaks': len(peak_positions),
        'max_mz': int(np.max(peak_positions)) if len(peak_positions) > 0 else 0,
        'min_mz': int(np.min(peak_positions)) if len(peak_positions) > 0 else 0,
        'total_intensity': float(np.sum(peak_intensities)),
        'max_intensity': float(np.max(peak_intensities)) if len(peak_intensities) > 0 else 0.0,
        'mean_intensity': float(np.mean(peak_intensities)) if len(peak_intensities) > 0 else 0.0
    }


def batch_similarity_matrix(features: list,
                           min_shared_peaks: int = DEFAULT_MIN_SHARED_PEAKS,
                           cosine_threshold: float = DEFAULT_COSINE_THRESHOLD) -> np.ndarray:
    """
    Calculate pairwise similarity matrix for a batch of features (optional optimization).
    
    This function is useful for calculating similarities between many features efficiently,
    though for graph construction the pairwise approach is usually sufficient.
    
    Inputs:
        features (list): List of feature dictionaries with MS/MS arrays
        min_shared_peaks (int): Minimum shared peaks requirement
        cosine_threshold (float): Minimum similarity threshold
        
    Outputs:
        np.ndarray: Symmetric similarity matrix (n_features x n_features)
    """
    n_features = len(features)
    similarity_matrix = np.zeros((n_features, n_features), dtype=np.float32)
    
    # Fill upper triangle of matrix
    for i in range(n_features):
        for j in range(i + 1, n_features):
            sim = calculate_spectral_similarity(
                features[i], features[j], min_shared_peaks, cosine_threshold
            )
            similarity_matrix[i, j] = sim
            similarity_matrix[j, i] = sim  # Symmetric
    
    # Diagonal is 1.0 for features that have MS/MS data
    for i in range(n_features):
        if has_msms_data(features[i]):
            similarity_matrix[i, i] = 1.0
    
    return similarity_matrix


def validate_feature_arrays(feature: Dict[str, Any]) -> bool:
    """
    Validate that a feature has properly formatted MS/MS arrays.
    
    Inputs:
        feature (Dict[str, Any]): Feature dictionary to validate
        
    Outputs:
        bool: True if arrays are valid, False otherwise
    """
    try:
        # Check required fields exist
        if 'msms_peaks' not in feature or 'msms_intensities' not in feature:
            return False
        
        peaks = feature['msms_peaks']
        intensities = feature['msms_intensities']
        
        # Check array types and shapes
        if not isinstance(peaks, np.ndarray) or not isinstance(intensities, np.ndarray):
            return False
        
        if peaks.dtype != bool or intensities.dtype != np.float32:
            return False
        
        if peaks.shape != intensities.shape:
            return False
        
        # Check has_msms flag consistency
        has_msms_flag = feature.get('has_msms', False)
        actual_has_msms = np.any(peaks)
        
        if has_msms_flag != actual_has_msms:
            logger.warning(f"Inconsistent has_msms flag: {has_msms_flag} vs actual: {actual_has_msms}")
            return False
        
        return True
        
    except Exception as e:
        logger.error(f"Error validating feature arrays: {e}")
        return False


# Configuration functions
def set_default_min_shared_peaks(new_min_peaks: int) -> None:
    """Set the default minimum shared peaks requirement."""
    global DEFAULT_MIN_SHARED_PEAKS
    DEFAULT_MIN_SHARED_PEAKS = new_min_peaks
    logger.info(f"Set DEFAULT_MIN_SHARED_PEAKS to {new_min_peaks}")


def set_default_cosine_threshold(new_threshold: float) -> None:
    """Set the default cosine similarity threshold."""
    global DEFAULT_COSINE_THRESHOLD
    DEFAULT_COSINE_THRESHOLD = new_threshold
    logger.info(f"Set DEFAULT_COSINE_THRESHOLD to {new_threshold}")


def get_config() -> Dict[str, Any]:
    """Get current configuration."""
    return {
        'min_shared_peaks': DEFAULT_MIN_SHARED_PEAKS,
        'cosine_threshold': DEFAULT_COSINE_THRESHOLD
    }