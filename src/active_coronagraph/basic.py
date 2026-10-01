import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
from matplotlib import animation
from astropy.io import fits
import os
from scipy import ndimage
import cmath
import math

def wfits(im, fname):
    """wfits - write im to file fname, automatically overwriting any old file"""
    from astropy.io import fits
    hea = fits.PrimaryHDU(im)
    hea.writeto(fname, overwrite=True)


def phase_2pi(z):
    phases = np.angle(z)  # Get the phase of each complex number
    phases[phases < 0] += 2 * np.pi  # Adjust phases to be in the range [0, 2*pi]
    return phases

def FFT(c,padsize):
    """FFT - carry out a complex Fourier transform (with optional padding)
            c - the input 2D Complex numpy array
            pad - integer multiplier for the padding/sampling
            Returns the `Complex` FFT padded array"""
    from numpy.fft import fft2,fftshift,ifft2,ifftshift
    # psfA = fftshift(fft2(ifftshift(padcplx(c,pad))))
    psfA = fftshift(fft2(ifftshift(pad_array(c, padsize))))
    return psfA

def IFFT(cb,pad=5):
    """IFFT - carry out the complex Fourier transform (with optional padding)
            and return the FFT padded array"""
    from numpy.fft import fft2,fftshift,ifft2,ifftshift
    psfB = fftshift(ifft2(ifftshift(padcplx(cb,pad))))
    return psfB

def make_complex_pupil(A, m_vor, m_vor_1, shift= True):
    '''create a complex pupil'''
    ph = np.zeros_like(A)
    X, Y = np.meshgrid(np.linspace(-np.pi, np.pi, np.shape(A)[0], endpoint=False),
                       np.linspace(-np.pi, np.pi, np.shape(A)[1], endpoint=False))

    if shift:
        A_shift = A * np.exp( 1j  * (ph + 0.5 * X + 0.5 * Y))

    Field = FFT(A_shift, 1)

    Field1 = FFT(IFFT(Field * m_vor_1, 1), 1)


    if shift:
        A_1 = IFFT(Field1*m_vor_1, 1) * np.exp( 1j  * (ph - 0.5 * X - 0.5 * Y))

    # A = IFFT(Field, 1)


    return A_1


def padcplx(c,pad=5):
    """padcplx - puts a `Complex` array into the centre of a zero-filled `Complex` array
               pad is an integer defining the padding multiplier for the output array """
    (nx, ny) = c.shape
    bignx = nx*pad
    bigny = ny*pad
    big_c = np.zeros((bignx,bigny),dtype=complex)

    dx = int((nx * (pad-1)) / 2)
    dy = int((ny * (pad-1)) / 2)

    big_c[dx:dx+nx,dy:dy+ny] = c
    return(big_c)


def pad_array(array, pad_size = 1000):
    """
    Pads a given 2D array with a specified value to make it 1000x1000 in size.

    Parameters:
    - array: 2D numpy array to pad
    - pad_value: Value to pad with (default is 0)

    Returns:
    - padded_array: 2D numpy array of size 1000x1000
    """
    pad_value = 0
    original_shape = array.shape
    if original_shape[0] > pad_size or original_shape[1] > pad_size:
        raise ValueError("Input array dimensions must be less than or equal to 1000x1000")

    # Calculate the required padding for each dimension
    pad_height = pad_size - original_shape[0]
    pad_width = pad_size - original_shape[1]

    # Calculate padding for top, bottom, left, and right
    pad_top = pad_height // 2
    pad_bottom = pad_height - pad_top
    pad_left = pad_width // 2
    pad_right = pad_width - pad_left

    # Apply the padding
    padded_array = np.pad(array, ((pad_top, pad_bottom), (pad_left, pad_right)), 'constant', constant_values=pad_value)

    return padded_array



def circle_mask(im, xc, yc, rcirc):
    """circle_mask - function that takes the input 2D array 'im' that evaluates the equation
            (x-x_c)^2 + (y-y_c)^2 < r^2 with circle center coordinates (x_c, y_c) and a radius 'r'
            as input parameters and return a mask array with the same shape as 'im'."""
    ny, nx = im.shape
    y,x = np.mgrid[0:nx,0:ny]
    r = np.sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc))
    return ( (r < rcirc))

def Isum(c):
    """Isum - routine that calculate the summed (real) intensity of a complex amplitude input image, c"""
    return (np.sum(c * c.conjugate()).real)

def zoom(im, trim):
    """ Zoom into the center of the image"""

    cx, cy = (np.array(im.shape) / 2).astype('int')
    return im[cx-trim:cx+trim, cy-trim:cy+trim]


def phi_ramp(im, npx, npy):
    """phi_ramp - make an X-ramp, making the 0 be in the middle,
            make it -0.5 to +0.5 and make it -pi to +pi. """
    ny, nx = im.shape
    ly = np.linspace(-0.5, 0.5, ny) * np.pi * npy * 2
    lx = np.linspace(-0.5, 0.5, nx) * np.pi * npx * 2

    x, y = np.meshgrid(lx, ly)
    return (x + y)


def xaxis_energy(im, xr):
    """xaxis_energy - takes a complex image im, cuts out the line of pixels from -xr to +xr
            across the middle of the image and returns the energy per pixel along that line"""
    ny, nx = im.shape
    xc = int((nx-1) / 2.)
    yc = int((ny-1) / 2.)
    sli = im[yc,xc-xr:xc+xr+1]
    return((np.abs(sli)**2))


def rings(im, x, y, r_rings):
    """rings(im, x, y, r_rings) - makes a mask for rings of different radii specified in a numpy array
       r_rings and central point (x,y): invalid values are -1
    - first, second, third... ring has value (0, 1, 2, ....)
    - ring 0 is from r_rings[0] to r_rings[1]
    - ring 1 is from r_rings[1] to r_rings[2] """

    im_rings = np.zeros_like(im) - 1.

    # make r_inner and r_outer
    r_inner = r_rings[0:-1]
    r_outer = r_rings[1:]

    r, t = r_theta(im_rings, x, y)
    for i, (rin, rout) in enumerate(zip(r_inner, r_outer)):
        im_rings[(r >= rin) * (r < rout)] = i

    r_middle = (r_inner + r_outer) / 2.
    return (im_rings, r_middle)


def r_theta(im, xc, yc):
    """r_theta - make a radius mask and return the radius rr and the angle phi for point (xc,yc)"""
    ny, nx = im.shape
    yp, xp = np.mgrid[0:ny,0:nx]
    yp = yp - yc
    xp = xp - xc
    rr = np.sqrt(np.power(yp,2.) + np.power(xp,2.))
    phi = np.arctan2(yp, xp)
    return(rr, phi)

def rms(im,rr,r_in,r_out):
    '''calculate the median in each annuli'''
    rms = []
    for i in range(len(r_in)):
        annuli = (rr * scale > r_in[i]) * (rr * scale < r_out[i])
        rms.append(np.std(med[annuli == True]))

    return rms

def make_annuli_profile(foldername, fitsfile, pic_index, nsamp):

    os.chdir(foldername)
    hdul = fits.open(fitsfile)
    data = hdul[0].data[pic_index]
    print(data.shape)
    dims = data.shape[0]

    annuli = np.arange(0, nsamp*10, nsamp/10)
    r_in = annuli[:-1]
    r_out = annuli[1:]
    r_mid = (r_in+r_out)/2/nsamp
    F_energy_med = []
    rr= r_theta(data, (dims-1)/2, (dims-1)/2 )[0]

    for i in range(len(r_mid)):
        mask = (rr>r_in[i])*(rr<r_out[i])
        F_energy_med.append(np.median(data[mask]))

    # displC(F_star/peak, f"final_focal_plane_coro_mask={name}_sample={nsamp}_{shift}.")
    fig,ax = plt.subplots(figsize = (8, 6))
    ax.plot(r_mid, F_energy_med)
    plt.ylabel("contrast/energy")
    plt.xlabel(r"$\lambda$/D")
    plt.yscale("log")
    # plt.title(f"Energy_annuli_med_{nsamp}_{name}_shift_{shift}_binary_{binary}")
    plt.savefig(f"{fitsfile}_annuli_profile.png")
    plt.show()
    os.chdir("..")

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)


def round_to_even(number):
    # Round to the nearest integer first
    rounded = round(number)

    # Check if the rounded number is even
    if rounded % 2 == 0:
        return rounded
    else:
        # If not even, round to the nearest even number
        # If the decimal part is 0.5, Python will round to the nearest even integer
        # If the number is odd, add or subtract 1 to make it even
        return rounded + 1 if number > rounded else rounded - 1


