"""
Module for handling 2D k-point grids.

This module provides functionality for generating and managing 2D k-point
grids needed for CBS calculations.
"""

import numpy as np
from typing import Tuple, List, Optional


class KPointGrid:
    """
    2D k-point grid for CBS calculations.
    
    Parameters
    ----------
    nk1 : int
        Number of k-points along first direction
    nk2 : int
        Number of k-points along second direction
    k1 : float
        Shift along first direction
    k2 : float
        Shift along second direction
    b1 : np.ndarray
        First reciprocal lattice vector
    b2 : np.ndarray
        Second reciprocal lattice vector
    """
    
    def __init__(
        self,
        nk1: int,
        nk2: int,
        k1: float,
        k2: float,
        b1: np.ndarray,
        b2: np.ndarray
    ):
        self.nk1 = nk1
        self.nk2 = nk2
        self.k1 = k1
        self.k2 = k2
        self.b1 = b1
        self.b2 = b2
        
        self.nkpts = nk1 * nk2
        
        # Generate grid
        self.xyk = np.zeros((self.nkpts, 2))
        self.wkpt = np.zeros(self.nkpts)
        
        self._generate_grid()
        
    def _generate_grid(self):
        """Generate the 2D k-point grid."""
        # Weight per k-point
        weight = 1.0 / self.nkpts
        
        ik = 0
        for i1 in range(self.nk1):
            for i2 in range(self.nk2):
                # Fractional coordinates
                xk1 = (i1 + self.k1) / self.nk1
                xk2 = (i2 + self.k2) / self.nk2
                
                # Store in crystal coordinates
                self.xyk[ik, 0] = xk1
                self.xyk[ik, 1] = xk2
                self.wkpt[ik] = weight
                
                ik += 1
    
    @classmethod
    def from_custom_kpoints(
        cls,
        kpoints: np.ndarray,
        weights: np.ndarray,
        b1: np.ndarray,
        b2: np.ndarray
    ):
        """
        Create a KPointGrid from custom k-point list.
        
        Parameters
        ----------
        kpoints : np.ndarray
            Array of k-point coordinates with shape (nkpts, 2)
            in crystal coordinates
        weights : np.ndarray
            Array of k-point weights with shape (nkpts,)
        b1 : np.ndarray
            First reciprocal lattice vector
        b2 : np.ndarray
            Second reciprocal lattice vector
            
        Returns
        -------
        KPointGrid
            K-point grid object with custom k-points
        """
        # Create a dummy instance with minimal grid
        instance = cls.__new__(cls)
        instance.nk1 = 0  # Indicate custom k-points
        instance.nk2 = 0
        instance.k1 = 0.0
        instance.k2 = 0.0
        instance.b1 = b1
        instance.b2 = b2
        
        # Set custom k-points
        instance.nkpts = kpoints.shape[0]
        instance.xyk = kpoints.copy()
        instance.wkpt = weights.copy()
        
        return instance
                
    def get_kpoint(self, ik: int) -> np.ndarray:
        """
        Get k-point coordinates in crystal units.
        
        Parameters
        ----------
        ik : int
            K-point index
            
        Returns
        -------
        np.ndarray
            K-point coordinates [k1, k2]
        """
        return self.xyk[ik]
        
    def get_weight(self, ik: int) -> float:
        """
        Get k-point weight.
        
        Parameters
        ----------
        ik : int
            K-point index
            
        Returns
        -------
        float
            K-point weight
        """
        return self.wkpt[ik]
        
    def get_all_kpoints(self) -> np.ndarray:
        """
        Get all k-points.
        
        Returns
        -------
        np.ndarray
            Array of shape (nkpts, 2) with k-point coordinates
        """
        return self.xyk
        
    def get_all_weights(self) -> np.ndarray:
        """
        Get all k-point weights.
        
        Returns
        -------
        np.ndarray
            Array of k-point weights
        """
        return self.wkpt
