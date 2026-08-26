from hcipy import *
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



def make_entrance_pupil(dim, d, nor = True):

    """Function to build a binary pupil"""

    im = np.zeros((dim, dim))
    ampli = circle_mask(im, (dim - 1) // 2, (dim - 1) // 2, d / 2)
    phase = np.zeros_like(im)
    A = ampli * np.exp(1j * phase)
    if nor is True:
        A /= np.sqrt(Isum(A) / 100)  # normalized the intensity of A to 100 intensity unit.
    else:
        A

    return A

def casted_spider(p1, p2, spider_width):
    return lambda grid: make_spider(p1, p2, spider_width)(grid).astype(np.float64)


def spider_arms(dims, pupil_grid, normalized = True):

    pupil_diameter = 4  # meter
    spider_width = 1.958e-2  # meter
    central_obscuration_ratio = 0.25
    spider_offset = np.array([0, 0.37251])  # meter

    if normalized:
        spider_width /= pupil_diameter
        spider_offset = [x / pupil_diameter * 100 for x in spider_offset]
        pupil_diameter = dims
        spider_width *= dims

    mirror_edge1 = (pupil_diameter / (2 * np.sqrt(2)), pupil_diameter / (2 * np.sqrt(2)))
    mirror_edge2 = (-pupil_diameter / (2 * np.sqrt(2)), pupil_diameter / (2 * np.sqrt(2)))
    mirror_edge3 = (pupil_diameter / (2 * np.sqrt(2)), -pupil_diameter / (2 * np.sqrt(2)))
    mirror_edge4 = (-pupil_diameter / (2 * np.sqrt(2)), -pupil_diameter / (2 * np.sqrt(2)))


    print(np.shape(spider_offset))
    # spider1 = make_spider(spider_offset, mirror_edge1, spider_width)(pupil_grid)
    # spider2 = make_spider(spider_offset, mirror_edge2, spider_width)(pupil_grid)
    # spider3 = make_spider([-v for v in spider_offset], mirror_edge3, spider_width)(pupil_grid)
    # spider4 = make_spider([-v for v in spider_offset], mirror_edge4, spider_width)(pupil_grid)

    spider1 = casted_spider(spider_offset, mirror_edge1, spider_width)
    spider2 = casted_spider(spider_offset, mirror_edge2, spider_width)
    spider3 = casted_spider([-v for v in spider_offset], mirror_edge3, spider_width)
    spider4 = casted_spider([-v for v in spider_offset], mirror_edge4, spider_width)

    spider1 = evaluate_supersampled(spider1, pupil_grid, 8)
    spider2 = evaluate_supersampled(spider2, pupil_grid, 8)
    spider3 = evaluate_supersampled(spider3, pupil_grid, 8)
    spider4 = evaluate_supersampled(spider4, pupil_grid, 8)


    return spider1 * spider2 * spider3 * spider4




def make_aperture(dim,oversample, obstruction = True, spider = False, obstruction_ratio = 0.25, lyot = False, lyot_fraction = None, nsamp = 10, nor = True, phase = None):
    """build a non binary pupil or lyot stop"""

    """if lyot is true, then the code is set to be building the lyot stop"""
    if lyot:
        dims = 100*nsamp
        dim = dim * lyot_fraction
    else:
        dims = dim

    pupil_grid = make_pupil_grid(dims,dims)

    '''
    if obstruction is set to be trphase = cmath.phase(z)
    if phase < 0:
        phase += 2 * math.pi
    return phaseue, we make add obstrcution of the second mirror '''
    if obstruction:
        aperture = make_obstructed_circular_aperture(dim, obstruction_ratio)

    else:
        aperture = make_circular_aperture(dim)

    '''
    use supersampling to evaluate this aperture to partially suppress sampling artefacts
    '''
    aperture = evaluate_supersampled(aperture, pupil_grid, oversample)
    if spider:
        aperture *= spider_arms(dims, pupil_grid)
        # plt.imshow(np.reshape(spider_arms(dims, pupil_grid), (dims, dims)))
        # plt.show()

    # imshow_field(aperture)
    # plt.title("aperture or lyot stop")
    # plt.colorbar()
    # plt.show()
    # print(np.shape(aperture)[0])
    A = np.reshape(aperture, (dims, dims))


    if lyot:
        cube = np.array([A])
        # wfits(cube, f"lyot_stop_large_cross_nsamp_{nsamp}_obstruction_{obstruction}_lyot_fraction_{lyot_fraction}.fits")
        return A

    if phase is None:
        phase = np.zeros_like(A)
    A = A * np.exp(1j * phase)

    if nor == True:
        return (A / np.sqrt(Isum(A) / 1))
    else:
        return A



def focal_mask_new(name, dim, nsamp, fpm_sam, wavelength, charge, opt_params=None, greyscale=None, cal_factor=1, fill_factor=False, noise_level = None):
    """
    Simulates a focal phase mask and rescales it to the native numerical simulation grid.

    Parameters:
        name (str)        : Name of the focal phase mask (e.g., "FQPM", "ACM", "vortex", etc.).
        dim (int)         : Total number of I/D of the point of view.
        nsamp (int)       : Sampling number of the numerical simulation grid.
        fpm_sam (int)     : Sampling number of the SLM grid.
        wavelength (float): Wavelength for phase calculation.
        mask_params (dict, optional): Dictionary containing mask-specific parameters:
            - "charge" (int)       : Charge for vortex/ACM mask.
            - "rotate" (bool)      : Whether to rotate the vortex mask.
            - "roddier_d" (float)  : Diameter for Roddier mask.
            - "plot_diff" (bool)   : Whether to plot phase difference.
        greyscale (int, optional)  : Number of grayscale levels (for quantization).
        cal_factor (float, optional): Calibration factor for phase scaling.
        fill_factor (bool, optional): If True, every 10th row/column of `m` is set to zero.

    Returns:
        np.array: The final complex focal phase mask.
    """

    if opt_params is None:
        opt_params = {"rotate":False, "roddier_d":1.06}  # Ensure it's a dictionary if not provided

    dims = dim * nsamp
    phase_mask = Phase_masks(fpm_sam=fpm_sam, pupil_size=100)

    # Generate the appropriate phase mask
    if name == "FQPM":
        phase = phase_mask.FQPM()
    elif name == "ACM":
        phase = phase_mask.ACM(charge)
    elif name == "vortex":
        phase = phase_mask.vortex(charge, opt_params.get("rotate"))
    elif name == "roddier":
        phase = phase_mask.roddier(opt_params.get("roddier_d"))
    elif name == "dual_zone":
        phase = phase_mask.dual_zone(a1=0.515, a2=0.705, z1=0.47, z2=0.92)
    else:
        raise ValueError(f"Unknown focal mask type: {name}")

    # Rescale phase mask to match numerical simulation grid
    scale_factor = nsamp / fpm_sam
    m_vor_phase_repeat = np.repeat(np.repeat(phase, scale_factor, axis=0), scale_factor, axis=1)

    # Ensure roddier mask is binary
    if name == "roddier":
        thres = 0.8
        m_vor_phase_repeat[m_vor_phase_repeat < thres] = 0
        m_vor_phase_repeat[m_vor_phase_repeat > thres] = np.pi

    # Apply grayscale quantization if required
    if greyscale is not None:
        bin_edges = np.linspace(0, 2 * np.pi * cal_factor, 2 ** greyscale, endpoint=True)
        m_vor_phase_repeat = bin_edges[np.digitize(m_vor_phase_repeat * cal_factor, bin_edges) - 1]

    # Final phase mask with normalization
    m = np.ones((dims, dims))

    m_vor_phase_repeat = m_vor_phase_repeat - np.pi

    m_vor_phase_repeat = np.angle(m * np.exp(1j * m_vor_phase_repeat)) + np.pi

    m_vor_phase_repeat = m_vor_phase_repeat * (1 / wavelength) - np.pi

    # Apply fill_factor: Set every 10th row and column of `m` to zero
    if fill_factor:
        print("yes")
        m[::10, :] = 0  # Set every 10th row to zero
        m[:, ::10] = 0  # Set every 10th column to zero

    # If plot_diff is enabled, compute and display the phase difference
    if opt_params.get("plot_diff", False):
        plt.figure(figsize=(10, 10))
        plt.imshow(m_vor_phase_repeat, cmap='jet')
        plt.colorbar()
        plt.title(f"Phase Mask: {name} (Charge={mask_params.get('charge', 'N/A')})")
        plt.show()

    '''add gaussian noise to the phase mask'''
    if noise_level is not None:
        noise = np.random.normal(0, noise_level * 2 * np.pi, m_vor_phase_repeat.shape)
        m_vor_phase_repeat = m_vor_phase_repeat + noise

    # Generate the final complex mask
    m_vor = m * np.exp(1j * m_vor_phase_repeat)

    # Save output
    # wfits([np.angle(m_vor)+np.pi], f"m_vor_{name}_{charge}_nsamp={nsamp}_fpm_sam={fpm_sam}_rotate={opt_params.get('rotate')}_greyscale={greyscale}_wl_{wavelength}_cal_{cal_factor}.fits")

    return m_vor






def make_lyot_stop(dim, fraction = 0.95, sec_fraction = 1.1, obstruction = True, nsamp = 10):

    """ simulate the lyot stop component in the coronagraph."""

    N_vor = make_aperture(dim, 8, obstruction= obstruction, obstruction_ratio=0.25*sec_fraction/fraction, lyot = True, lyot_fraction= fraction, nsamp = nsamp)

    return N_vor


def coro_high_sam(A, m , N, nsamp, fmp_sam, name, charge = None, lyot_sum = None, obstruction = True, shift = True, save = False, get_lyot = False):

    """
    simulated coronograph taking a telescope pupil 'A', transmissive masks 'm' and 'N'
    and a sampling number nsamp.

    A -- entrace pupil
    m -- focal phase mask
    nsamp -- the sampling number of the numerical simulation grid at the entrance pupil, with the unit [(numerical simulation pixel)/(I/D)]
    fpm_sam -- the sampling number of the SLM grid
    name -- the name of the focal phase mask
    charge -- the charge of the vortex focal phase mask. If the phase mask is not vortex then the charge will set to be None.
    lyot_sum -- the integration value of the lyot plane without focal phase mask but with lyot stop.
    obstruction -- if obstruction is True, a second mirror will be added.
    shift -- if shift is true, shift the phase of the entrance pupil by half of the pixel, to get everything symmetrical in the pupil.
    the center of the pixel in the entrance pupil.
    save -- If save is true, the image of the middle stages in the coronagraph simulation are saved as the fits file. Otherwise, no fits file
    will be saved.
    get_lyot -- If get_lyot is True, return the intensity integration value of the lyot plane.

    """


    if shift:
        """ shift half pixel in the first focal plane. """
        ph = np.zeros_like(A)

        X, Y = np.meshgrid(np.linspace(-np.pi, np.pi, np.shape(A)[0], endpoint=False),
                           np.linspace(-np.pi, np.pi, np.shape(A)[1], endpoint=False))
        A = A * np.exp(1j * (ph - 0.5 / nsamp * X - 0.5 / nsamp * Y))


    """first focal plane"""
    A = FFT(A, nsamp*100)
    # psf_sum = np.sum(zoom(np.abs(A) ** 2, 50))
    # print(f"pdf_sum_{psf_sum}")
    # if save:
    #     displC(A, f"first_focal_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fmp_sam}_obstruction_{obstruction}.")


    """roddier energy throughput"""
    # mask = np.angle(m)
    # mask[mask<3] = 0
    # mask[mask>3] = 1
    # A_ene = A * mask
    # ratio = np.sum(np.abs(A_ene)**2)/np.sum(np.abs(A)**2)
    # print(f"ratio_{ratio}")
    """plane after focal mask"""
    A = m*A

    A = IFFT(A, 1)

    if save:
        if lyot_sum is not None:
            print(os.getcwd())
            """ save the normalized intensity of the lyot plane """
            cube = np.array([np.abs(A) ** 2 / lyot_sum])
            lyot_stop_plane = (np.abs(A) ** 2 / lyot_sum)[400:600,400:600]
            wfits(lyot_stop_plane, f"lyot_plane_profile_{nsamp}_{charge}_{fmp_sam}.fits")
            # plt.imshow(lyot_stop_plane)
            # plt.show()


        #     displC((np.abs(A) ** 2 / lyot_sum), f"lyot_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fmp_sam}_obstruction_{obstruction}_normalized.")
        #
        # else:
        #     displC(A, f"lyot_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fmp_sam}_obstruction_{obstruction}.")

    """plane after lyot stop"""
    A_sum = np.sum(np.abs(A)**2)
    A = N*A
    res_ene = np.sum(np.abs(A) ** 2)/A_sum
    if save:
        if lyot_sum is not None:
            after_lyot_stop_plane = (np.abs(A) ** 2 / lyot_sum)[400:600, 400:600]
            wfits(after_lyot_stop_plane, f"after_lyot_plane_profile_{nsamp}_{charge}_{fmp_sam}.fits")
            # plt.imshow(after_lyot_stop_plane)
            # plt.show()

    #         cube = np.array([np.abs(A) ** 2 / lyot_sum])
    #         wfits(cube, f"after_lyot_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fmp_sam}_obstruction_{obstruction}_normalized.fits")
    #     else:
    #         displC(A, f"after_lyot_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fmp_sam}_obstruction_{obstruction}.")

    if get_lyot:
        return np.sum((np.abs(A))**2)

    else:
        A = FFT(A,nsamp*100)
        return A, res_ene


def simple_coro(dim, name, wavelength, nsamp, fpm_sam, lyot_stop = None, offset = None, charge = None, binary = False, obstruction = True, spider = False, get_lyot = True,
                rotate = False, greyscale = True, roddier_d = 1.06, cal_factor = 1, phase = None, noise_level = None, fill_factor = False):

    """

    Apply the entrance pupil, focal phase mask and lyot stop into the coronagraph simulator.

    dim  -- the total number of I/D of the point of view.
    name -- the name of the focal phase mask
    nsamp -- the sampling number of the numerical simulation grid at the entrance pupil, with the unit [(numerical simulation pixel)/(I/D)]
    fpm_sam -- the sampling number of the SLM grid.
    charge -- the charge of the vortex focal phase mask. If the phase mask is not vortex then the charge will set to be None.
    lyot_stop -- the type of lyot stop, e.g. small, large or zernike.
    offset -- If offset is not None, a phase ramp is added to the entrance pupil to simulate the planet signal.
    binary -- if binary is True, the entrance pupil will be simulated with only binary value, otherwise, there will be smooth change at the edge
    of the pupil.
    obstruction -- if obstruction is True, a second mirror will be added.
    get_lyot  -- if get_lyot is True, the lyot plane with chosen focal phase mask is normalized with the sum of lyot plane without focal phase mask. And the planes before
    and after lyot are saved.

    """


    '''build a binary or non binary pupil'''


    if binary:
        A = make_entrance_pupil(dim)
    else:
        A = make_aperture(round_to_even(dim/wavelength), 8, obstruction = obstruction, spider = spider, obstruction_ratio=0.25, phase = phase)


    if offset is not None:
        '''adding a phase ramp to the entrance pupil to simulate the planet signal'''
        A_1 = A * np.exp(1j * phi_ramp(A, offset[0]*(round_to_even(dim/wavelength)/dim), offset[1]*(round_to_even(dim/wavelength)/dim)))
    else:
        A_1 = A

    lyot_stop_para = {"small_cross": [0.98, 1.2], "large_cross": [0.95, 1.3], "zernike_cross": [0.98, 1.05], "None": [1,1.2],
                      "None_1_1": [1,1], "None_1_1.1": [1,1.1], "None_1_1.3": [1,1.3], "None_1_1.4": [1, 1.4], "None_1_1.5": [1, 1.5]}

    """ set different sizes of lyot stop."""
    N_vor = make_lyot_stop(round_to_even(dim/wavelength), fraction = lyot_stop_para[lyot_stop][0], sec_fraction = lyot_stop_para[lyot_stop][1], obstruction = obstruction, nsamp = nsamp)



    '''if true, save the fits files of each step.'''
    save = False
    '''the integration intensity after lyot stop without coro, set the initial value to None.'''
    lyot_sum = None
    shift = True
    if get_lyot:
        m_vor = np.ones_like(N_vor)
        """ get the integration intensity after lyot stop without phase mask but with lyot stop. """
        lyot_sum = coro_high_sam(A, m_vor, N_vor, nsamp, fpm_sam, name = name, charge = charge, obstruction = obstruction, lyot_sum = lyot_sum, shift = shift, save = save, get_lyot = get_lyot)
        # print(lyot_sum)
        save = True

    """ set the focal phase mask to the chosen one."""
    # m_vor = focal_mask_new(name, dim, nsamp, fpm_sam, wavelength, charge, rotate = rotate, greyscale = greyscale, roddier_d = roddier_d, cal_factor = cal_factor)
    m_vor =  focal_mask_new(name, dim, nsamp, fpm_sam, wavelength, charge = charge, greyscale = greyscale, cal_factor = cal_factor, fill_factor = fill_factor, noise_level = noise_level)

    '''if get_lyot is True and lyot_sum is not None, then save will be True at this stage and the plane images of each stages 
    except the Final Science Plane, and planes before and after lyot stop are normalized with the integration of lyot plane values without phase mask '''
    final_focal = coro_high_sam(A_1, m_vor, N_vor, nsamp, fpm_sam, name = name, charge = charge, obstruction = obstruction, lyot_sum = lyot_sum, shift = shift, save = save)[0]


    if not get_lyot:

        m_vor = np.ones_like(m_vor)
        if fill_factor:
            m_vor[::10, :] = 0  # Set every 10th row to zero
            m_vor[:, ::10] = 0  # Set every 10th column to zero

        if offset != None:
            N_vor_nor = make_lyot_stop(int(dim / wavelength), fraction=lyot_stop_para["None_1_1"][0],
                                   sec_fraction=lyot_stop_para["None_1_1"][1], obstruction=obstruction)
            wo_mask_sum = np.sum(np.abs(coro_high_sam(A, m_vor, N_vor_nor, nsamp, fpm_sam, name = name, charge = charge, obstruction = obstruction, lyot_sum = lyot_sum, shift = shift, save = False)[0])**2)
            cube = np.array([np.sum(np.abs(final_focal)**2)]/wo_mask_sum)
            np.savetxt("throughput.txt", [cube])
            return np.max(cube)
        else:
            wo_mask_peak = np.max(np.abs(coro_high_sam(A_1, m_vor, N_vor, nsamp, fpm_sam, name = name, charge = charge, obstruction = obstruction, lyot_sum = lyot_sum, shift = shift, save = False)[0])**2)
            cube = np.array(np.abs(final_focal)**2/wo_mask_peak)

            """ save the normalized intensity in the final focal plane in the fits file. """
            wfits([cube], f"final_focal_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fpm_sam}_obstruction_{obstruction}.fits")
            # displC(cube[0], "final")
            return (np.abs(final_focal) ** 2 / wo_mask_peak)


    # else:
    #     m_vor = focal_mask(name, dim, nsamp, charge)
    #     A_star = np.copy(A)
    #     # propagation through the system with coronagraph
    #     B_star, C_star, D_star, E_star, F_star = coro(A_star, m_vor, N_vor, nsamp)
    #
    #     # propagation through the system without coronagraph
    #     m = np.ones_like((m_vor))
    #     B_star_1, C_star_1, D_star_1, E_star_1, F_star_1 = coro(A_star, m, N_vor, nsamp)
    #     print(np.sum(np.abs(E_star_1)**2))
    #     displC(E_star_1, f"after_lyot_plane_coro_mask={name}_sample={nsamp}_obstruction_{obstruction}_no_coro.")
    #
    #     if single:
    #         if charge == None:
    #             name = name
    #         else:
    #             name = f"{name}_{charge}"
    #         displC(A, f"entrance_pupil_{name}_{nsamp}_binary_{binary}.")
    #         displC(zoom(B_star, 4), f"first_focal_plane_coro_mask={name}_sample={nsamp}_zoom_in_{shift}.")
    #         displC(B_star, f"first_focal_plane_coro_mask={name}_sample={nsamp}_obstruction_{obstruction}.")
    #         displC(D_star, f"lyot_plane_coro_mask={name}_sample={nsamp}_obstruction_{obstruction}.")
    #         displC(E_star, f"after_lyot_plane_coro_mask={name}_sample={nsamp}_obstruction_{obstruction}.")
    #         # peak = np.max(np.real(F_star_1*np.conjugate(F_star_1)))
    #         peak = np.max(np.abs(F_star_1))
    #         print(peak)
    #         # print(peak)
    #         displC(zoom(F_star/peak,100), f"final_focal_plane_coro_mask={name}_sample={nsamp}_obstruction_{obstruction}_without_square.")
    #         # cube = np.array([np.real(F_star*np.conjugate(F_star))/peak])
    #         cube = np.array([np.abs(F_star)**2/peak**2])
    #         wfits(cube, f"final_focal_plane_coro_mask={name}_sample={nsamp}_obstruction_{obstruction}.fits")
    #
    #
    #         #makeing the annuli profile
    #         # annuli = np.arange(0, nsamp*10, nsamp/10)
    #         # r_in = annuli[:-1]
    #         # r_out = annuli[1:]
    #         # r_mid = (r_in+r_out)/2/nsamp
    #         # F_energy_med = []
    #         # rr= r_theta(np.abs(F_star), (dim*nsamp-1)/2, (dim*nsamp-1)/2 )[0]
    #         # peak = np.max(np.abs(F_star_1)**2)
    #         # for i in range(len(r_mid)):
    #         #     mask = (rr>r_in[i])*(rr<r_out[i])
    #         #     F_energy_med.append(np.median((np.abs(F_star)**2)[mask]))
    #         #
    #         # # displC(F_star/peak, f"final_focal_plane_coro_mask={name}_sample={nsamp}_{shift}.")
    #         # fig,ax = plt.subplots(figsize = (8, 6))
    #         # ax.plot(r_mid, F_energy_med/peak)
    #         # plt.ylabel("contrast/energy")
    #         # plt.xlabel(r"$\lambda$/D")
    #         # plt.yscale("log")
    #         # plt.title(f"Energy_annuli_med_{nsamp}_{name}_shift_{shift}_binary_{binary}")
    #         # plt.savefig(f"Energy_annuli_med_{nsamp}_{name}_shift_{shift}_binary_{binary}.png")
    #         # plt.show()
    #
    #     else:
    #         return F_star, F_star_1


def res_ene(dim, name, wavelength, nsamp, fpm_sam, lyot_stop=None, offset=None, charge=None, binary=False,
                    obstruction=True, get_lyot=True, rotate=False, greyscale=True, roddier_d=1.06):

        """

        Apply the entrance pupil, focal phase mask and lyot stop into the coronagraph simulator.

        dim  -- the total number of I/D of the point of view.
        name -- the name of the focal phase mask
        nsamp -- the sampling number of the numerical simulation grid at the entrance pupil, with the unit [(numerical simulation pixel)/(I/D)]
        fpm_sam -- the sampling number of the SLM grid.
        charge -- the charge of the vortex focal phase mask. If the phase mask is not vortex then the charge will set to be None.
        lyot_stop -- the type of lyot stop, e.g. small, large or zernike.
        offset -- If offset is not None, a phase ramp is added to the entrance pupil to simulate the planet signal.
        binary -- if binary is True, the entrance pupil will be simulated with only binary value, otherwise, there will be smooth change at the edge
        of the pupil.
        obstruction -- if obstruction is True, a second mirror will be added.
        get_lyot  -- if get_lyot is True, the lyot plane with chosen focal phase mask is normalized with the sum of lyot plane without focal phase mask. And the planes before
        and after lyot are saved.

        """

        '''build a binary or non binary pupil'''

        if binary:
            A = make_entrance_pupil(dim)
        else:
            A = make_aperture(round_to_even(dim / wavelength), 8, obstruction=obstruction, obstruction_ratio=0.25)

        if offset is not None:
            '''adding a phase ramp to the entrance pupil to simulate the planet signal'''
            A_1 = A * np.exp(1j * phi_ramp(A, offset[0] * (round_to_even(dim / wavelength) / dim),
                                           offset[1] * (round_to_even(dim / wavelength) / dim)))
        else:
            A_1 = A

        lyot_stop_para = {"small_cross": [0.98, 1.2], "large_cross": [0.95, 1.3], "zernike_cross": [0.98, 1.05],
                          "None": [1, 1.2],
                          "None_1_1": [1, 1], "None_1_1.1": [1, 1.1], "None_1_1.3": [1, 1.3], "None_1_1.4": [1, 1.4],
                          "None_1_1.5": [1, 1.5]}

        """ set different sizes of lyot stop."""
        N_vor = make_lyot_stop(round_to_even(dim / wavelength), fraction=lyot_stop_para[lyot_stop][0],
                               sec_fraction=lyot_stop_para[lyot_stop][1], obstruction=obstruction, nsamp=nsamp)

        '''if true, save the fits files of each step.'''
        save = False
        '''the integration intensity after lyot stop without coro, set the initial value to None.'''
        lyot_sum = None
        shift = True
        if get_lyot:
            m_vor = np.ones_like(N_vor)
            """ get the integration intensity after lyot stop without phase mask but with lyot stop. """
            lyot_sum = coro_high_sam(A, m_vor, N_vor, nsamp, fpm_sam, name=name, charge=charge, obstruction=obstruction,
                                     lyot_sum=lyot_sum, shift=shift, save=save, get_lyot=get_lyot)
            # print(lyot_sum)
            save = True

        """ set the focal phase mask to the chosen one."""
        m_vor =  focal_mask_new(name, dim, nsamp, fpm_sam, wavelength, charge = charge, greyscale = greyscale)

        '''if get_lyot is True and lyot_sum is not None, then save will be True at this stage and the plane images of each stages 
        except the Final Science Plane, and planes before and after lyot stop are normalized with the integration of lyot plane values without phase mask '''
        res_ene = coro_high_sam(A_1, m_vor, N_vor, nsamp, fpm_sam, name=name, charge=charge,
                                    obstruction=obstruction, lyot_sum=lyot_sum, shift=shift, save=save)[1]


        return res_ene



















