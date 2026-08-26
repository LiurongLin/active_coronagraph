import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from basic import *
from phase_masks import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
from matplotlib import animation
from astropy.io import fits
from functools import partial
import os
from make_plot import *
from hcipy import *

mpl.rcParams['figure.dpi'] = 300


class Phase_masks:

    def __init__(self, fpm_sam, pupil_size):
        """
        Initialize the phase mask parameters.

        """
        self.fpm_sam = fpm_sam
        self.pupil_size = pupil_size

    def dual_zone(self, show=False, a1=0.515, a2=0.705, z1=0.47, z2=0.92):
        """
        Create a dual zone phase mask.

        Parameters:
        show (bool): Whether to display the resulting mask. Default is False.

        Returns:
        np.ndarray: A padded matrix representing the dual zone phase mask.
        """
        # Matrix size in pixels
        matrix_size = int(self.pupil_size * self.fpm_sam)

        # Radii in pixels
        small_radius = self.fpm_sam * a1
        large_radius = self.fpm_sam * a2

        # Phase values
        small_circle_value = 2 * np.pi * z1
        ring_value = 2 * np.pi * z2

        # Build coordinates centered at the geometric center (N-1)/2
        center = (matrix_size - 1) / 2.0
        y = np.arange(matrix_size) - center
        x = np.arange(matrix_size) - center
        yy, xx = np.meshgrid(y, x, indexing='ij')
        r = np.hypot(xx, yy)

        # Allocate and fill
        matrix = np.zeros((matrix_size, matrix_size), dtype=float)
        matrix[r <= small_radius] = small_circle_value
        mask_ring = (r > small_radius) & (r <= large_radius)
        matrix[mask_ring] = ring_value

        if show:
            plt.imshow(matrix, cmap='viridis', origin='upper')
            plt.colorbar(label='Phase (rad)')
            plt.title(f"Dual-zone phase mask ({matrix_size}×{matrix_size})")
            plt.savefig("dual_zone_mask.png", dpi=150, bbox_inches='tight')
            plt.show()

        return matrix

    def dual_zone_old(self, show=False, a1 = 0.475, a2 = 0.8,  z1 = 0.34, z2 = 0.74):
        """
        Create a dual zone phase mask.

        Parameters:
        show (bool): Whether to display the resulting mask. Default is False.

        Returns:
        np.ndarray: A padded matrix representing the dual zone phase mask.
        """
        # Constants
        matrix_size = self.pupil_size * self.fpm_sam  # Size of the matrix
        small_radius = self.fpm_sam * a1 # Radius of the smaller circle
        large_radius = self.fpm_sam * a2 # Radius of the larger circle (outer boundary of the ring)
        center = matrix_size // 2  # Center of the matrix (middle pixel)
        pi_value = np.pi

        # Values for each region
        small_circle_value = 2 * z1 * pi_value
        ring_value = 2 * z2 * pi_value

        # Create a matrix filled with zeros
        matrix = np.zeros((matrix_size, matrix_size))

        # Create grid indices
        y, x = np.ogrid[:matrix_size, :matrix_size]

        # Calculate the Euclidean distance of each point from the center
        distance_from_center = np.sqrt((x - center) ** 2 + (y - center) ** 2)

        # Assign value to the small circle region (inside the smaller radius)
        matrix[distance_from_center <= small_radius] = small_circle_value

        # Assign value to the ring region (between small circle and larger circle)
        matrix[(distance_from_center > small_radius) & (distance_from_center <= large_radius)] = ring_value

        return matrix


    def FQPM(self, center=0):

        matrix_size = self.pupil_size * self.fpm_sam  # Size of the matrix
        c = (matrix_size - 1) / 2.0 + center
        if np.isclose(center, -0.5):
            c = (matrix_size - 1) / 2.0

        y, x = np.indices((matrix_size, matrix_size), dtype=float)
        masq0 = np.zeros((matrix_size, matrix_size))
        masq0[(x < c) & (y < c)] = np.pi
        masq0[(x >= c) & (y >= c)] = np.pi

        if np.isclose(center, -0.5):
            row_idx = 500
            col_idx = 500
            if 0 <= row_idx < matrix_size:
                masq0[row_idx, :] = np.pi / 2.0
            if 0 <= col_idx < matrix_size:
                masq0[:, col_idx] = np.pi / 2.0

        return masq0

    def vortex(self, charge, rotate = False, center = 0):

        matrix_size = self.pupil_size * self.fpm_sam  # Size of the matrix
        m_vor_grid = make_uniform_grid([matrix_size, matrix_size], [matrix_size, matrix_size], center = center)
        phase = charge * (m_vor_grid.as_('polar').theta + np.pi).reshape([matrix_size, matrix_size])
        phase = np.angle(np.exp(1j * phase)) + np.pi

        if rotate:
            phase = charge * ((m_vor_grid.as_('polar').theta + np.pi).reshape([matrix_size, matrix_size]) + (np.pi / 4))
            phase = np.angle(np.exp(1j * phase)) + np.pi

        return phase

    def ACM(self, charge, rotate = False):

        matrix_size = self.pupil_size * self.fpm_sam  # Size of the matrix
        m_vor_grid = make_uniform_grid([matrix_size, matrix_size], [matrix_size, matrix_size])
        phase = charge * (m_vor_grid.as_('polar').theta + np.pi).reshape([matrix_size, matrix_size])
        phase = 2.4048 * np.cos(phase)
        phase = np.angle(np.exp(1j * phase))+np.pi

        # if rotate:
        #     phase = charge * ((m_vor_grid.as_('polar').theta + np.pi).reshape([dims, dims]) + (np.pi / 4))
        #     phase = np.angle(np.exp(1j * phase)) + np.pi


        return phase

    def roddier(self, roddier_d):
        roddier_d = roddier_d * self.fpm_sam
        matrix_size = self.pupil_size * self.fpm_sam
        pupil_grid = make_pupil_grid( matrix_size,  matrix_size)
        mask = make_circular_aperture(roddier_d)
        mask = evaluate_supersampled(mask, pupil_grid, 8)
        phase = np.reshape(mask, ( matrix_size,  matrix_size)) * np.pi

        return phase
