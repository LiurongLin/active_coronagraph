

from .basic import *
from .coronagraphs import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
from matplotlib import animation
from astropy.io import fits
from functools import partial
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams,rc
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm




def make_contrast_plot_new():

    default_name = ["vortex", "vortex", "vortex", "vortex", "FQPM", "roddier", "dual_zone"]
    default_charge = [2, 4, 6, 8, None, None, None]
    fpm_sam_list = [100, 100]
    sam_list = [10, 100]

    # Initialize figure for subplots
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))  # 2 rows, 4 columns (7 plots + 1 empty)
    axes = axes.ravel()

    # Loop through the names and charges to generate file paths and plot data
    for i, (name, charge) in enumerate(zip(default_name, default_charge)):
        ax = axes[i]
        for (sam, fpm_sam) in zip(fpm_sam_list, sam_list):
            # Generate file path
            charge_str = f"_{charge}"
            contrast = np.loadtxt(
                f"ideal_coro_2rd_mirror_False/lyot_None_1_1_{name}{charge_str}_lambda_1.0_fpm_sam={fpm_sam}_binary_False_obstruction_False/"
                f"final_focal_plane_coro_mask={name}{charge_str}_sample={sam}_fpm_sam={fpm_sam}_obstruction_False_wl_1.0.txt"
            )

            # Load data
            # Replace this with your actual data loading method, e.g., np.loadtxt(file_path)
            offset = np.arange(0, 10, 0.1)

            # Plot the data
            ax.plot(offset[1:], contrast, label=f"fpm_sam={fpm_sam}", linewidth = 5)

        # Customize subplot
        ax.set_title(f"{name} (charge={charge})")
        ax.set_xlabel("offset")
        ax.set_ylabel("contrast")
        ax.set_ylim(1e-10,1e-2)
        ax.set_yscale("log")
        ax.legend()
        ax.grid(True)

    # Remove the last subplot if unused
    fig.delaxes(axes[-1])

    # Adjust layout and show plot
    plt.tight_layout()
    plt.savefig("contrast_summary_diff_sam.png")
    plt.show()

def make_plot_contrast_vs_res():
    sample_list = [2, 2 ** 2, 2 ** 3, 2 ** 4, 2 ** 5, 2 ** 6, 78, 88, 98]
    name = ["vortex_2", "vortex_4", "FQPM"]
    os.chdir("contrast_vs_res_txt")
    binary = [True, False]
    shift = [True, False]
    fig, ax = plt.subplots(2,3, figsize = (15,10))
    custom_ylim = (1e-14, 1e-1)
    plt.setp(ax, ylim=custom_ylim)
    # fig1, ax1 = plt.subplots(2, 3, figsize=(15, 10))
    for z in range(3):
        for i in shift:
            for j in binary:
                contrast_arr = np.loadtxt(f"contrast_vs_res_{name[z]}_01_09_binary_{j}_shift_{i}.txt")
                if i == True:
                    if j == True:
                        ax[0][z].plot(sample_list, contrast_arr[0], color = "r", label = r"binary = True, 1 $\lambda/D$")
                        ax[0][z].set_yscale('log')
                        ax[0][z].set_xlabel("sample number")
                        ax[0][z].set_ylabel("contrast")
                        ax[0][z].plot(sample_list, contrast_arr[1], color="b", label=r"binary = True, 10 $\lambda/D$")
                        ax[0][z].set_yscale('log')
                    else:
                        ax[0][z].plot(sample_list, contrast_arr[0], color="r", linestyle = "dashed", label = f"binary = False, 1 $\lambda/D$")
                        ax[0][z].set_yscale('log')
                        ax[0][z].plot(sample_list, contrast_arr[1], color="b", linestyle="dashed",
                                      label=f"binary = False, 10 $\lambda/D$")
                        ax[0][z].set_yscale('log')
                else:
                    if j == True:
                        ax[1][z].plot(sample_list, contrast_arr[0], color = "r", label = r"binary = True, 1 $\lambda/D$")
                        ax[1][z].set_yscale('log')
                        ax[1][z].set_xlabel("sample number")
                        ax[1][z].set_ylabel("contrast")
                        ax[1][z].plot(sample_list, contrast_arr[1], color="b", label=r"binary = True, 10 $\lambda/D$")
                        ax[1][z].set_yscale('log')
                    else:
                        ax[1][z].plot(sample_list, contrast_arr[0], color="r", linestyle = "dashed", label = f"binary = False, 1 $\lambda/D$")
                        ax[1][z].set_yscale('log')
                        ax[1][z].plot(sample_list, contrast_arr[1], color="b", linestyle="dashed",
                                      label=f"binary = False, 10 $\lambda/D$")
                        ax[1][z].set_yscale('log')
        ax[0][z].set_title(f"{name[z]}_shift")
        ax[1][z].set_title(f"{name[z]}_no_shift")

        # ax[0][z].legend()
        # ax[1][z].legend()
    plt.tight_layout()
    plt.legend()
    plt.savefig("contrast_vs_res_sum.png")
    plt.show()

    os.chdir("..")
    print(os.getcwd())


def make_pupil_diff_image():

    binary = make_entrance_pupil(100)
    binary_cut = binary[20][:]

    non_binary = make_aperture(100,8)
    non_binary_cut = non_binary[20][:]

    diff = binary - non_binary
    diff_cut = diff[20][:]

    pupil = [binary, non_binary, diff]
    cut = [binary_cut, non_binary_cut, diff_cut]

    name = ["binary", "non_binary", "difference"]

    ratio = 10
    f = plt.figure(figsize=(15,6))
    gs = plt.GridSpec(10+8, ratio*3+8)

    for i in range(3):
        ax_joint = f.add_subplot(gs[8:, ratio*(i):ratio*(i+1)])
        ax_marg_x = f.add_subplot(gs[6, ratio*(i)+2:ratio*(i+1)-1])
        if i > 0:
            plt.setp(ax_marg_x.get_yticklabels(), visible=False)
            plt.setp(ax_marg_x.get_yticklabels(minor=True), visible=False)

            plt.setp(ax_joint.get_yticklabels(), visible=False)
            plt.setp(ax_joint.get_yticklabels(minor=True), visible=False)

        ax_joint.imshow(np.abs(pupil[i]))
        ax_joint.axhline(y = 20 , color = "r")
        ax_marg_x.set_title(name[i])
        # binary_cut = np.fliplr(np.abs(binary)).diagonal()
        # print(binary_cut)
        ax_marg_x.plot(np.abs(cut[i]))
        ax_marg_x.set_xlim(0,20)

        print("1")

    plt.savefig("binary_non_binary.png")
    plt.show()


def make_mask_image():
    fig, ax = plt.subplots(4, 2, figsize=(10, 15))

    mask_name = ["vortex", "vortex", "FQPM", "roddier"]
    charge = [2, 4, None, None]

    for i in range(4):
        mask = focal_mask(mask_name[i], 100, 32, charge[i])
        if i == 3:
            mask = zoom(mask, 100)
        amp = ax[i][0].imshow(np.abs(mask), vmax=np.pi, vmin=-np.pi)
        fig.colorbar(amp, ax=ax[i][0])
        ang = ax[i][1].imshow(np.angle(mask), vmax=np.pi, vmin=-np.pi)
        fig.colorbar(ang, ax=ax[i][1])

    plt.savefig("4_focal_mask.png")

    plt.show()


def make_annuli_profile(foldername, fitsfile, pic_index, nsamp):
    os.chdir("ideal_coro_2rd_mirror_True/")

    os.chdir(foldername)
    hdul = fits.open(fitsfile)
    data = hdul[0].data[pic_index]
    print(data.shape)
    dims = data.shape[0]

    annuli = np.arange(0, nsamp * 10, nsamp / 10)
    r_in = annuli[:-1]
    r_out = annuli[1:]
    r_mid = (r_in + r_out) / 2 / nsamp
    F_energy_med = []
    rr = r_theta(data, (dims - 1) / 2, (dims - 1) / 2)[0]

    for i in range(len(r_mid)):
        mask = (rr > r_in[i]) * (rr < r_out[i])
        F_energy_med.append(np.median(data[mask]))

    # displC(F_star/peak, f"final_focal_plane_coro_mask={name}_sample={nsamp}_{shift}.")
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(r_mid, F_energy_med)
    np.savetxt(f"{fitsfile[:-4]}.txt", F_energy_med)
    plt.ylabel("contrast/energy")
    plt.xlabel(r"$\lambda$/D")
    plt.yscale("log")
    # plt.title(f"Energy_annuli_med_{nsamp}_{name}_shift_{shift}_binary_{binary}")
    plt.savefig(f"{fitsfile}_annuli_profile.png")
    plt.show()

    os.chdir("../..")


def annuli_mean(name, nsamp, fpm_sam, lyot_stop, charge, obstruction, wl = None, rotate = False, greyscale = True, broadband = False, noise_level= None,
                folder_name = None, cus_name = None):

    """
    create a list of annuli with the same size of the ring equal to 1/10 of the numerical sampling number but with increasing inner circle
    and outer circle radius.
    get the median of intensities inside the rings and store them in a list.

    name -- the name of the focal phase mask
    nsamp -- the sampling number of the numerical simulation grid at the entrance pupil, with the unit [(numerical simulation pixel)/(I/D)]
    fpm_sam -- the sampling number of the SLM grid.
    charge -- the charge of the vortex focal phase mask. If the phase mask is not vortex then the charge will set to be None.
    lyot_stop -- the type of lyot stop, e.g. small, large or zernike.

    save the list of intensity median into a text file.
    """
    if cus_name:
        hdul = fits.open(cus_name)
    else:
        root_dir = f"ideal_coro_2rd_mirror_{obstruction}"
        if noise_level is not None:
            root_dir = f"{root_dir}_noise_{noise_level:g}"
        if folder_name is None:
             os.chdir(f"{root_dir}/")
        else:
            os.chdir(f"{root_dir}{folder_name}/")

        if rotate:
            os.chdir(f"lyot_{lyot_stop}_{name}_{charge}_lambda_{wl}_sam_{nsamp}_fpm_sam={fpm_sam}_binary_False_obstruction_{obstruction}_rotate_True")
        elif greyscale != None:
            if broadband:
                os.chdir(f"lyot_{lyot_stop}_{name}_{charge}_lambda_1.0_fpm_sam={fpm_sam}_binary_False_obstruction_{obstruction}_greyscale_{greyscale}")
            else:
                os.chdir(f"lyot_{lyot_stop}_{name}_{charge}_lambda_{wl}_fpm_sam={fpm_sam}_binary_False_obstruction_{obstruction}_greyscale_{greyscale}")

        else:
            if broadband:
                os.chdir(
                    f"lyot_{lyot_stop}_{name}_{charge}_lambda_1.0_fpm_sam={fpm_sam}_binary_False_obstruction_{obstruction}")
            else:
                os.chdir(
                    f"lyot_{lyot_stop}_{name}_{charge}_lambda_{wl}_fpm_sam={fpm_sam}_binary_False_obstruction_{obstruction}")

        if broadband:
            hdul = fits.open(f"final_focal_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fpm_sam}_obstruction_{obstruction}_bb.fits")
        else:
            hdul = fits.open(f"final_focal_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fpm_sam}_obstruction_{obstruction}.fits")

    data = hdul[0].data[0]
    hdul.close()
    print(data.shape)
    dims = data.shape[0]

    annuli = np.arange(0, nsamp * 15, nsamp/10)

    r_in = annuli[:-1]
    r_out = annuli[1:]
    r_mid = (r_in + r_out) / 2 / nsamp
    F_energy_med = []
    rr = r_theta(data, (dims - 1) / 2, (dims - 1) / 2)[0]

    for i in range(len(r_mid)):
        mask = (rr > r_in[i]) * (rr < r_out[i])
        # fig, ax = plt.subplots(3,1, figsize = (15,5))
        # ax[0].imshow(rr > r_in[i])
        # ax[1].imshow(rr < r_out[i])
        # ax[2].imshow(mask)
        # plt.suptitle(f"{r_in[i]}_{r_out[i]}")
        # plt.show()
        F_energy_med.append(np.median(data[mask]))
    if cus_name:
        np.savetxt("final_psf_contrast_ghost.txt", F_energy_med)
    else:
        if broadband:
            np.savetxt(
                f"final_focal_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fpm_sam}_obstruction_{obstruction}_bb.txt",
                F_energy_med)
        else:
            np.savetxt(f"final_focal_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fpm_sam}_obstruction_{obstruction}_wl_{wl}.txt", F_energy_med)
        print(f"final_focal_plane_coro_mask={name}_{charge}_sample={nsamp}_fpm_sam={fpm_sam}_obstruction_{obstruction}.txt")

        os.chdir("../..")
    return F_energy_med

def lyot_plane_mean(name, charge, fpm_sam):

    dir = (f"ideal_coro_2rd_mirror_False/lyot_None_{name}_{charge}_sam_100_fpm_sam={fpm_sam}_binary_False_obstruction_False_greyscale_8/")
    os.chdir(dir)

    data = np.loadtxt(f"after_lyot_plane_profile_100_{charge}_{fpm_sam}.txt")

    annuli = np.arange(1, np.shape(data)[0]/2+1,1)

    r_in = annuli[:-1]
    r_out = annuli[1:]
    r_mid = (r_in + r_out) / 2
    F_energy_med = []
    rr = r_theta(data, np.shape(data)[0] / 2, np.shape(data)[0] / 2)[0]

    for i in range(len(r_mid)):
        mask = (rr > r_in[i]) * (rr < r_out[i])
        # fig, ax = plt.subplots(3,1, figsize = (15,5))
        # ax[0].imshow(rr > r_in[i])
        # ax[1].imshow(rr < r_out[i])
        # ax[2].imshow(mask)
        # plt.suptitle(f"{r_in[i]}_{r_out[i]}")
        # plt.show()
        F_energy_med.append(np.median(data[mask]))

    np.savetxt(
        f"after_lyot_stop_plane_coro_mask={name}_{charge}_sample=100_fpm_sam={fpm_sam}_obstruction_False.txt",
        F_energy_med)
    print(f"after_lyot_stop_plane_coro_mask={name}_{charge}_sample=100_fpm_sam={fpm_sam}_obstruction_False.txt")

    os.chdir("../..")



def make_lyot_plane_plot(name, charge, poster = True):
    plt.style.use('dark_background')
    fig, ax = plt.subplots(1,6, figsize = (30,7), sharex=True, sharey=True)
    rc('font', weight='bold', size = 30)
    plt.rcParams['savefig.dpi'] = 500
    fig.subplots_adjust(wspace=0.001)
    title_list = ["Vortex n=2", "Vortex n=4", "Vortex n=6", "Vortex n=8", "FQPM", "Roddier"]

    for i in range(len(name)):
        dir = f"ideal_coro_2rd_mirror_True/lyot_None_1_1_{name[i]}_{charge[i]}_sam_100_fpm_sam=100_binary_False_obstruction_True_greyscale_8"
        os.chdir(dir)
        lsp = np.loadtxt(f"lyot_plane_profile_100_{charge[i]}_100.txt")

        ax[i].imshow(lsp, vmax = 1e-4, vmin = 1e-7, cmap = "inferno")
        # if charge[i] == None:
        #     ax[i].set_title(f"{name[i]}", weight="bold", fontsize=30)
        # else:
        #     ax[i].set_title(f"{name[i]}_{charge[i]}", weight = "bold", fontsize = 30)
        circle1 = plt.Circle(xy=(np.shape(lsp)[0]//2-0.5, np.shape(lsp)[0]//2-0.5), radius=25/2, color='red', linewidth=2, fill=False)
        ax[i].add_patch(circle1)
        circle2 = plt.Circle(xy=(np.shape(lsp)[0] // 2 - 0.5, np.shape(lsp)[0] // 2 -0.5), radius=25*1.2/2, color='darkgreen', linewidth=2,
                             fill=False)
        ax[i].add_patch(circle2)
        circle3 = plt.Circle(xy=(np.shape(lsp)[0] // 2 -0.5, np.shape(lsp)[0] // 2-0.5), radius=25*1.3/2, color='darkorange', linewidth=2,
                             fill=False)
        ax[i].add_patch(circle3)

        # ax[i].xaxis.set_tick_params(labelsize=30)
        xticklabels = plt.get(ax[i], 'xticklabels')
        plt.setp(xticklabels, fontsize=30, weight='bold')
        if poster == False:
            ax[i].set_title(title_list[i], fontsize = 35, weight = "bold")


        ax[i].set_xlim(10, 140)
        ax[i].set_ylim(10, 140)
        os.chdir("../..")
    yticklabels = plt.get(ax[0], 'yticklabels')
    plt.setp(yticklabels, fontsize=30, weight='bold')
    # ax[0].yaxis.set_tick_params(labelsize=30)
    # yticklabels = getp(gca(), 'yticklabels')
    # setp(yticklabels, fontsize=14, weight=‘bold’)

    ax[2].set_xlabel("Pixel", weight = "bold", fontsize = 30)
    # fig.supxlabel("pixel", weight = "bold", fontsize = 30)
    # plt.yscale("log")

    # plt.tight_layout()
    plt.savefig("lyot_plane_poster.png", transparent = True)
    plt.show()

def make_lyot_plane_diff(name, charge, fpm_sam_1, fpm_sam_2):

    dir = (f"ideal_coro_2rd_mirror_False/lyot_None_{name}_{charge}_sam_100_fpm_sam={fpm_sam_1}_binary_False_obstruction_False_greyscale_8/")
    os.chdir(dir)

    data_1 = np.loadtxt(f"after_lyot_plane_profile_100_{charge}_{fpm_sam_1}.txt")

    os.chdir("../..")

    dir = (f"ideal_coro_2rd_mirror_False/lyot_None_{name}_{charge}_sam_100_fpm_sam={fpm_sam_2}_binary_False_obstruction_False_greyscale_8/")
    os.chdir(dir)

    data_2 = np.loadtxt(f"after_lyot_plane_profile_100_{charge}_{fpm_sam_2}.txt")

    os.chdir("../..")

    plt.figure(figsize = (10,10))
    plt.imshow(data_2-data_1)
    print(np.sum(data_2-data_1))
    plt.xlabel("pixel")
    plt.ylabel("pixel")
    plt.colorbar()
    plt.savefig(f"after_lyot_stop_diff_{name}_{charge}_{fpm_sam_2}_{fpm_sam_1}.png")
    plt.show()




def make_contrast_plot(obstruction, lyot_stop):

    plt.style.use("seaborn-v0_8-talk")
    name = ["vortex", "vortex", "vortex", "FQPM", "roddier"]
    tp_i = {"large_cross": 0.801, "small_cross":0.904, "zernike_cross": 0.922}
    charge = [2, 4, 6, None, None]
    offset = np.arange(0, 10, 0.1)
    nsamp = 10

    os.chdir(f"ideal_coro_2rd_mirror_{obstruction}/")


    fig, ax = plt.subplots(3,1, figsize = (8, 15))
    fig.subplots_adjust(hspace=0)
    for i in range(len(name)):

        os.chdir(f"lyot_{lyot_stop}_{name[i]}_{charge[i]}_sam_{nsamp}_binary_False_obstruction_{obstruction}")
        contrast = np.loadtxt(f"final_focal_plane_coro_mask={name[i]}_{charge[i]}_sample={nsamp}_obstruction_{obstruction}.txt")
        tp = np.loadtxt(f"{lyot_stop}_tp_{name[i]}_{charge[i]}_obstruction={obstruction}.txt")

        ax[0].plot(offset[1:], contrast, label = f"{name[i]}_{charge[i]}")
        ax[0].set_ylabel("contrast")
        ax[0].set_ylim(1e-12,1)
        ax[0].set_yscale("log")
        ax[1].plot(offset[1:], tp[1:], label = f"{name[i]}_{charge[i]}")
        ax[1].set_ylabel("throughput")
        ax[1].set_ylim(0,1)
        ax[2].plot(offset[1:], contrast/tp[1:], label = f"{name[i]}_{charge[i]}")
        ax[2].set_ylim(1e-10, 1)
        ax[2].set_yscale("log")
        ax[2].set_xlabel(r"$\lambda$/D")
        ax[2].set_ylabel("contrast/throughput")
        os.chdir("..")

    plt.legend()
    plt.tight_layout()
    plt.savefig(f"contrast_tp_lyot_{lyot_stop}_obstruction={obstruction}.png")
    plt.show()

    print(os.getcwd())

    os.chdir("..")


def initialize_plot_styles():

    plt.style.use("seaborn-v0_8-talk")
    rc('font', weight='bold', size=30)

    #


def set_plot_parameters(data):
    """
        Set the global parameters for the plot.
    """
    name = ["vortex", "vortex", "vortex", "vortex", "FQPM", "roddier", "dual_zone", "ACM", "ACM", "ACM"]
    title_list = ["Vortex, n=2", "Vortex, n=4", "Vortex, n=6", "Vortex, n=8", "FQPM", "Roddier", "Dual Zone", "ACM, n=2", "ACM, n=4", "ACM, n=6"]
    yaxis_dic = {"contrast": "log(contrast)", "SNR": r"$\eta_{p}$/$\sqrt{\eta_{s}}$", "tp": "throughput"}
    charge = [2, 4, 6, 8, None, None, None, 2, 4, 6]
    if data == "contrast":
        offset = np.arange(0, 10, 0.1)
    else:
        offset = np.arange(0, 6, 0.1)

    return name, title_list, charge, offset, yaxis_dic[data]


def determine_parameters(diff_lyot, diff_fpm_sam, diff_nsamp, greyscale, diff_bb):
    """
    Determine the parameters based on the input flags.
    """
    # Default parameters
    para_len = 1
    lyot_stop = ["None_1_1"]
    fpm_sam = ["10"]
    nsamp = ["100"]
    greyscale_list = [64]


    if diff_lyot:
        lyot_stop = ["None_1_1", "None_1_1.1", "None", "None_1_1.3", "None_1_1.4", "None_1_1.5"]
        # lyot_stop = ["None_1_1",  "None",  "None_1_1.4"]
        para_len = len(lyot_stop)
    elif diff_fpm_sam:
        fpm_sam = ["10", "100"]
        para_len = len(fpm_sam)
    elif greyscale:
        greyscale_list = ["8", "64"]
        para_len = len(greyscale_list)
    elif diff_nsamp:
        nsamp = ["10", "100"]
        para_len = len(nsamp)
    elif diff_bb:
        para_len = 2

    # Ensure all lists have the same length
    fpm_sam = fpm_sam * (para_len // len(fpm_sam))
    nsamp = nsamp * (para_len // len(nsamp))
    greyscale_list = greyscale_list * (para_len // len(greyscale_list))
    lyot_stop = lyot_stop * (para_len // len(lyot_stop))

    return para_len, lyot_stop, fpm_sam, nsamp, greyscale_list

def plot_data(ax, title, offset, data, tp, contrast, SNR, label, color):
    if data == "tp":
        ax.plot(offset, tp, label=label, color = color)
        ax.set_ylim(0, 1)
        ax.set_ylabel("throughput")
    elif data == "SNR":
        ax.plot(offset, np.log10(SNR), label=label, linewidth = 4, color = color)
    elif data == "contrast":
        ax.plot(offset[1:], np.log10(contrast[:99]), label=label, linewidth = 4, color = color)
        ax.set_ylim(-10, -2)
        ax.set_yticks(np.arange(-10, -2, 1))
        ax.set_xticks(np.arange(0, 10, 2))


    ax.set_title(title, fontsize=40, weight="bold")
    ax.grid(True, color="white")
    ax.xaxis.set_tick_params(labelsize=30)
    ax.yaxis.set_tick_params(labelsize=30)


def save_plot(ang, data, nsamp, fpm_sam, obstruction, greyscale, diff_bb, broadband, diff_lyot):
    """
    Save the plot to a file.
    """
    if data == "SNR":
        plt.savefig(f"SNR_diff_ls_nsamp_{nsamp}_fpm_sam={fpm_sam}_ang_{ang}.png")
    elif data == "tp":
        plt.savefig(f"tp_diff_ls_nsamp_{nsamp}_fpm_sam={fpm_sam}.png")
    elif data == "contrast" and diff_bb:
        plt.savefig(f"contrast_sam_{nsamp}_fpm_sam={fpm_sam}_new_obstruction_{obstruction}_greyscale_{greyscale}_poster_bb.png", transparent=True)
    elif data == "contrast" and broadband:
        plt.savefig(f"contrast_sam_{nsamp}_fpm_sam={fpm_sam}_new_obstruction_{obstruction}_greyscale_{greyscale}_poster_bb_diff_lyot_{diff_lyot}.png")
    elif data == "contrast":
        plt.savefig(f"contrast_sam_{nsamp}_fpm_sam={fpm_sam}_new_obstruction_{obstruction}_greyscale_{greyscale}_poster_diff_lyot_{diff_lyot}.png")
    plt.show()


def broadband_w_o_obstruction_contrast():
    name, title_list, charge, offset, yaxis = set_plot_parameters()
    plt.style.use('dark_background')
    fig, ax = plt.subplots(1, 6, figsize=(30, 10), sharex=True, sharey=True)
    for i in range(len(name)):
        num_noise = np.loadtxt(
            f"ideal_coro_2rd_mirror_False/lyot_large_cross_{name[i]}_{charge[i]}_sam_100_fpm_sam=100_binary_False_obstruction_False/final_focal_plane_coro_mask={name[i]}_{charge[i]}_sample=100_fpm_sam=100_obstruction_False.txt")
        contrast_obstruction = np.loadtxt(
            f"ideal_coro_2rd_mirror_True/lyot_None_1_1_{name[i]}_{charge[i]}_lambda_1_fpm_sam=10_binary_False_obstruction_True_greyscale_8/final_focal_plane_coro_mask={name[i]}_{charge[i]}_sample=10_fpm_sam=10_obstruction_True_bb.txt")
        if name[i] != "roddier" and "dual_zone":
            contrast_obstruction = contrast_obstruction - num_noise[:len(contrast_obstruction)]  # subtract numerical noise
        contrast_without_obstruction = np.loadtxt(
            f"ideal_coro_2rd_mirror_False/lyot_None_1_1_{name[i]}_{charge[i]}_lambda_1_fpm_sam=10_binary_False_obstruction_False_greyscale_8/final_focal_plane_coro_mask={name[i]}_{charge[i]}_sample=10_fpm_sam=10_obstruction_False_bb.txt")
        if name[i] != "roddier":
            contrast_without_obstruction = contrast_without_obstruction - num_noise[:len(contrast_without_obstruction)]  # subtract numerical noise

        plot_data(ax[i], title_list[i], offset, "contrast", None, contrast_obstruction, None, 'with obstrution')
        plot_data(ax[i], title_list[i], offset, "contrast", None, contrast_without_obstruction, None, 'without obstrution')

    fig.supxlabel(f"{yaxis}", fontsize=30, weight="bold")
    fig.supxlabel(r"offset[$\lambda$/D]", fontsize=30, weight="bold")
    plt.legend(loc="best", fontsize=30)
    plt.tight_layout()


    plt.savefig('broadband_w_o_obstruction_contrast.png')
    plt.show()


def make_summary_plot(data, ang, diff_lyot=False, diff_fpm_sam=False, diff_nsamp=False, diff_bb= False, obstruction=True, greyscale=True, broadband=False):
    """
    Main function to create TP, SNR, and contrast plots.
    """
    plt.style.use('dark_background')
    dim = 100
    initialize_plot_styles()
    name, title_list, charge, offset, yaxis = set_plot_parameters(data)
    para_len, lyot_stop, fpm_sam, nsamp, greyscale_list = determine_parameters(diff_lyot, diff_fpm_sam, diff_nsamp, greyscale, diff_bb)
    print(greyscale_list)
    fig, ax = plt.subplots(1, len(name), figsize=(35, 10), sharex=True, sharey=True)
    broadband_list = [True, False]
    lyot_size_list = [1, 1.1, 1.2, 1.3, 1.4, 1.5]
    contrast_ind = np.arange(0,60,1)
    # colormap = cm.hot
    # colors = [colormap(i) for i in np.linspace(0.2, 0.95, para_len)]

    colors = [
        '#00FFFF',  # Cyan
        '#FF00FF',  # Magenta
        '#FFFF00',  # Yellow
        '#FF0000',  # Red
        '#00FF00',  # Green
        '#0000FF',  # Blue
        '#FFA500',  # Orange
        '#90EE90',  # Light Green
        '#ADD8E6',  # Light Blue
        '#D3D3D3',  # Light Gray
        '#FFC0CB',  # Pink
        '#800080',  # Purple
        '#40E0D0',  # Turquoise
        '#00FF00'  # Lime (same as Green for variation)
    ]

    for i in range(len(name)):
        # Change to the directory for the specific obstruction
        if name[i] != "dual_zone":
            num_noise = np.loadtxt(f"ideal_coro_2rd_mirror_False/lyot_None_1_1_{name[i]}_{charge[i]}_lambda_1.0_fpm_sam=100_binary_False_obstruction_False/final_focal_plane_coro_mask={name[i]}_{charge[i]}_sample=100_fpm_sam=100_obstruction_False_wl_1.0.txt")
        os.chdir(f"ideal_coro_2rd_mirror_{obstruction}/")
        for k in range(para_len):
            print(nsamp)
            # Create the label for the plot

            label = {
                'greyscale': f"{greyscale_list[k]} bits",
                'diff_fpm_sam': rf"{fpm_sam[k]} pixel/ ($\lambda$/D)",
                'diff_nsamp': f"array size = {(dim*int(nsamp[k]))**2:.1e} pixel",
                # 'diff_bb': f"broadband: {broadband_list[k]}",
                'diff_lyot':f"D\'s/Ds = {lyot_size_list[k]}"
            }.get(
                'greyscale' if greyscale else
                'diff_fpm_sam' if diff_fpm_sam else
                'diff_nsamp' if diff_nsamp else
                # 'diff_bb' if diff_bb else
                'diff_lyot' if diff_lyot else
                'default', 'Default'
            )
            if name[i] == "FQPM":
                tp = np.loadtxt(
                    f"lyot_{lyot_stop[k]}_{name[i]}_{charge[i]}_sam_{nsamp[k]}_fpm_sam={fpm_sam[k]}_binary_False_obstruction_{obstruction}_greyscale_8_bb_ang_45.txt") if data != "contrast" else None
            else:
                tp = np.loadtxt(
                    f"lyot_{lyot_stop[k]}_{name[i]}_{charge[i]}_sam_{nsamp[k]}_fpm_sam={fpm_sam[k]}_binary_False_obstruction_{obstruction}_greyscale_8_bb_ang_{ang}.txt") if data != "contrast" else None
            if greyscale and greyscale_list[k] == 64:
                path = f"lyot_{lyot_stop[k]}_{name[i]}_{charge[i]}_lambda_1.0_fpm_sam={fpm_sam[k]}_binary_False_obstruction_{obstruction}"
            elif diff_bb:
                path = f"lyot_None_1_1_{name[i]}_{charge[i]}_lambda_1_fpm_sam=10_binary_False_obstruction_{obstruction}_greyscale_8"
            elif greyscale:
                path = f"lyot_{lyot_stop[k]}_{name[i]}_{charge[i]}_sam_{nsamp[k]}_fpm_sam={fpm_sam[k]}_binary_False_obstruction_{obstruction}_greyscale_{greyscale_list[k]}"
            elif diff_lyot:
                if broadband:
                    path = f"lyot_{lyot_stop[k]}_{name[i]}_{charge[i]}_lambda_1.0_fpm_sam=10_binary_False_obstruction_{obstruction}_greyscale_8"
                else:
                    path = f"lyot_{lyot_stop[k]}_{name[i]}_{charge[i]}_lambda_1.0_fpm_sam={fpm_sam[k]}_binary_False_obstruction_{obstruction}_greyscale_8"
            else:
                if broadband:
                    path = f"lyot_None_1_1_{name[i]}_{charge[i]}_lambda_1.0_fpm_sam=10_binary_False_obstruction_{obstruction}_greyscale_8"
                else:
                    path = f"lyot_None_1_1_{name[i]}_{charge[i]}_lambda_1.0_fpm_sam={fpm_sam[k]}_binary_False_obstruction_{obstruction}"
            os.chdir(path)

            # Load the data

            if data != "tp":
                print(os.getcwd())
                if (diff_bb and k == 0) or broadband:
                    print(os.getcwd())
                    contrast = np.loadtxt(
                        f"final_focal_plane_coro_mask={name[i]}_{charge[i]}_sample=100_fpm_sam={fpm_sam[k]}_obstruction_{obstruction}_bb.txt")
                else:
                    contrast = np.loadtxt(
                        f"final_focal_plane_coro_mask={name[i]}_{charge[i]}_sample=100_fpm_sam={fpm_sam[k]}_obstruction_{obstruction}_wl_1.0.txt")
                if name[i] != "roddier" and "dual_zone":
                    contrast = contrast - num_noise[:len(contrast)] # subtract numerical noise
            else:
                contrast = None
            SNR = tp / np.sqrt(contrast[contrast_ind]) if data == "SNR" else None


            # Plot the data
            plot_data(ax[i], title_list[i], offset, data, tp, contrast, SNR, label, colors[k])
            # if broadband and k == 0:
            #     pass
            # else:
            os.chdir("..")
        os.chdir("..")

    # Set super labels and save the plot
    # fig.supylabel("log(contrast)", fontsize=30, x=.01, weight="bold")
    fig.supylabel(f"{yaxis}", fontsize=30, weight="bold")
    fig.supxlabel(r"offset[$\lambda$/D]", fontsize=30, weight="bold")
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=30)
    plt.tight_layout()
    save_plot(ang, data, nsamp, fpm_sam, obstruction, greyscale, diff_bb, broadband, diff_lyot)


def make_TP_SNR_plot_archived(data, nsamp, diff_lyot = False, fpm_sam = None, diff_fpm_sam = False, obstruction = True, greyscale = True, grey_scale_diff = False):

    """
    making plots of throughput of different offset, signal to noise
    and contrast at different offsets.
    data -- "tp": throughput, "contrast": contrast, "SNR" : signal to noise.
    nsamp -- the sampling number of the numerical simulation grid at the entrance pupil, with the unit [(numerical simulation pixel)/(I/D)]
    fpm_sam -- the sampling number of the SLM grid.
    diff_lyot -- if diff_lyot is True, the final image includes all the simulated lyot stop options
    diff_fpm -- if diff_fpm_sam is True, the final image includes all the simulated focal phase mask sampling options
    obstruction -- if obstruction is True, a second mirror will be added.

    A image of desired parameter is plotted out.
    """
    plt.style.use("seaborn-v0_8-talk")
    name = ["vortex", "vortex", "vortex", "vortex", "FQPM", "roddier"]
    title_list = ["Vortex, n=2", "Vortex, n=4", "Vortex, n=6", "Vortex, n=8", "FQPM", "Roddier"]
    tp_i = {"large_cross": 0.801, "small_cross": 0.904, "zernike_cross": 0.922}
    charge = [2, 4, 6, 8, None, None]
    offset = np.arange(0, 10, 0.1)
    greyscale_list = [64, 8]

    if diff_lyot:
        lyot_stop = ["large_cross", "small_cross", "zernike_cross"]
        fpm_sam = np.repeat(fpm_sam, len(lyot_stop))
        para_len = len(lyot_stop)

    elif diff_fpm_sam:
        fpm_sam = ["10", "100"]
        para_len = len(fpm_sam)
        gretscale_list = ["64"]
        lyot_stop = np.repeat("None", len(fpm_sam))


    elif greyscale:
        greyscale_list = [8, 64]
        fpm_sam = ["10", "10"]
        para_len = len(greyscale_list)
        lyot_stop = np.repeat("None", len(fpm_sam))


    else:
        para_len = 1
        lyot_stop = "large_cross"


    # color = {"large_cross": "y", "small_cross": "b", "zernike_cross": "r"}

    rc('font', weight='bold',size = 30)
    fig, ax = plt.subplots(1, 6, figsize=(30, 10), sharex = True, sharey = True)
    for i in range(6): #for multiple grey scales.

        if data == "tp":
            ax[j][i].set_ylim(0, 1)

        os.chdir(f"ideal_coro_2rd_mirror_{obstruction}/")

        for k in range(para_len):
            if greyscale:
                label = f"{greyscale_list[k]} bits"
            else:
                label = rf"{fpm_sam[k]} pixel/ ($\lambda$/D)"
            # if i == 2:
            #     os.chdir(f"lyot_{lyot_stop[k]}_{name[i]}_{charge[i]}_sam_{nsamp}_fpm_sam={fpm_sam[k]}_binary_False_obstruction_{obstruction}_rotate_True")
            # else:
            if greyscale:
                if greyscale_list[k] == 64:
                    os.chdir(f"lyot_{lyot_stop[k]}_{name[i]}_{charge[i]}_sam_{nsamp}_fpm_sam={fpm_sam[k]}_binary_False_obstruction_{obstruction}")
                else:
                    os.chdir(f"lyot_{lyot_stop[k]}_{name[i]}_{charge[i]}_sam_{nsamp}_fpm_sam={fpm_sam[k]}_binary_False_obstruction_{obstruction}_greyscale_{greyscale_list[k]}")
            else:
                os.chdir(f"lyot_{lyot_stop[k]}_{name[i]}_{charge[i]}_sam_{nsamp}_fpm_sam={fpm_sam[k]}_binary_False_obstruction_{obstruction}")

            print(os.getcwd())
            # if i == 2:
            #     ax.set_title(f"{name[i]}_{charge[i]}_miniFQPM")
            # else:
            # if j == 0:

            ax[i].set_title(f"{title_list[i]}", fontsize=40, weight="bold")

            if data != "contrast":
                tp = np.loadtxt(f"{lyot_stop[k]}_tp_{name[i]}_{charge[i]}_fpm_sam={fpm_sam[k]}_obstruction={obstruction}.txt")

            if data != "tp":
                contrast = np.loadtxt(f"final_focal_plane_coro_mask={name[i]}_{charge[i]}_sample={nsamp}_fpm_sam={fpm_sam[k]}_obstruction_{obstruction}.txt")
                if data == "SNR":
                    SNR = tp/np.sqrt(contrast[contrast_ind])
                print(len(contrast))

            if data == "tp":
                ax[j][i].plot(offset, tp, label = label)
                ax[j][i].set_ylabel("throughput")
                # ax.set_ylim(0,1)

            elif data == "SNR":
                ax[j][i].plot(offset[1:], np.log10(SNR), label = label)
                # ax.set_ylim(-1,4)
                ax[j][i].set_ylabel("log(SNR)")

            elif data == "contrast":
                ax[i].plot(offset[1:], np.log10(contrast), label = label, linewidth = 6)

                if obstruction == True:
                    ax[i].set_ylim(-10, -1)
                else:
                    ax[i].set_ylim(-10, -2)

                # ax[j][i].set_ylabel("log(contrast)")

            os.chdir("..")
        ax[i].grid(True, color = "black")
        ax[i].xaxis.set_tick_params(labelsize=30)
        ax[i].yaxis.set_tick_params(labelsize=30)
        # ax[j][i].legend()
        os.chdir("..")


    fig.supylabel("log(contrast)", fontsize = 30, x = .01, weight = "bold")
    fig.supxlabel(r"offset[$\lambda$/D]", fontsize = 30, weight = "bold")
    plt.legend(loc="best", fontsize = 30)
    plt.tight_layout()
    if data == "SNR":
        plt.savefig(f"SNR_diff_ls_nsamp_{nsamp}_fpm_sam={fpm_sam}.png")
    elif data == "tp":
        plt.savefig(f"tp_diff_ls_nsamp_{nsamp}_fpm_sam={fpm_sam}.png")
    elif data == "contrast":
        print(os.getcwd())
        plt.savefig(f"contrast_diff_ls_{nsamp}_fpm_sam={fpm_sam}_new_obstruction_{obstruction}_greyscale_{greyscale}_poster.png", transparent = True)
    plt.show()


def plot_rotate():

    os.chdir(f"ideal_coro_2rd_mirror_False/")
    ''' Plot the final contrast profiles of vortex 4 and 8 with and without Mini FQPM in the center.'''

    plt.style.use("seaborn-v0_8-talk")
    name = ["vortex_4", "vortex_8"]
    fpm_sam = [5,10,20,100]
    vortex_num = 2
    fig, ax = plt.subplots(vortex_num, vortex_num, figsize=(10, 10))
    for i in range(vortex_num):
        for k in range(len(fpm_sam)):
            label = f"fpm_sam = {fpm_sam[k]}"
            os.chdir(f"lyot_large_cross_{name[i]}_sam_100_fpm_sam={fpm_sam[k]}_binary_False_obstruction_False")
            contrast = np.loadtxt(f"final_focal_plane_coro_mask={name[i]}_sample=100_fpm_sam={fpm_sam[k]}_obstruction_False.txt")
            offset = np.arange(0, 10, 0.1)
            ax[i][0].plot(offset[1:], np.log10(contrast), label = label)
            ax[i][0].set_title(f"{name[i]}")
            ax[i][0].set_xlabel("offset")
            ax[i][0].set_ylabel("contrast")
            ax[i][0].set_ylim(-13,0)
            ax[i][0].grid(True)
            ax[i][0].legend()
            os.chdir("..")

            os.chdir(f"lyot_large_cross_{name[i]}_sam_100_fpm_sam={fpm_sam[k]}_binary_False_obstruction_False_rotate_True")
            contrast = np.loadtxt(f"final_focal_plane_coro_mask={name[i]}_sample=100_fpm_sam={fpm_sam[k]}_obstruction_False.txt")
            offset = np.arange(0, 10, 0.1)
            ax[i][1].plot(offset[1:], np.log10(contrast), label = label)
            ax[i][1].set_title(f"{name[i]}_miniFQPM")
            ax[i][1].set_xlabel("offset")
            ax[i][1].set_ylabel("contrast")
            ax[i][1].set_ylim(-13, 0)
            ax[i][1].grid(True)
            ax[i][1].legend()
            os.chdir("..")

    plt.tight_layout()
    plt.savefig("vortex4_8_miniFQPM.png")
    plt.show()


def plot_greyscale():

    os.chdir(f"ideal_coro_2rd_mirror_False/")
    ''' Plot the final contrast profiles of vortex 2 roddier and FQPM with and without grey scale setting.'''

    plt.style.use("seaborn-v0_8-talk")
    name = ["vortex_2", "FQPM_None", "roddier_None"]
    grey_scale = [64, 8, 6]
    grey_scale_t = [64, 8, 6]
    color_l = ["b", "g", "r"]
    fig, ax = plt.subplots(len(name), 1, figsize=(5, 15))
    for i in range(len(name)):
        for j in range(len(grey_scale)):
            if grey_scale[j] == 64:
                os.chdir(f"lyot_large_cross_{name[i]}_sam_100_fpm_sam=10_binary_False_obstruction_False")
            else:
                os.chdir(f"lyot_large_cross_{name[i]}_sam_100_fpm_sam=10_binary_False_obstruction_False_greyscale_{grey_scale[j]}")
            contrast = np.loadtxt(f"final_focal_plane_coro_mask={name[i]}_sample=100_fpm_sam=10_obstruction_False.txt")
            label = f"greyscale_{grey_scale_t[j]}"
            offset = np.arange(0, 10, 0.1)
            ax[i].plot(offset[1:], np.log10(contrast), label = label, color = color_l[j])
            ax[i].set_title(f"{name[i]}")
            ax[i].set_xlabel("offset")
            ax[i].set_ylabel("contrast")
            ax[i].set_ylim(-10, 0)
            ax[i].grid(True)
            ax[i].legend()
            os.chdir("..")

            if grey_scale[j] == 64:
                os.chdir(f"lyot_large_cross_{name[i]}_sam_100_fpm_sam=5_binary_False_obstruction_False")
            else:
                os.chdir(f"lyot_large_cross_{name[i]}_sam_100_fpm_sam=5_binary_False_obstruction_False_greyscale_{grey_scale[j]}")
            contrast = np.loadtxt(f"final_focal_plane_coro_mask={name[i]}_sample=100_fpm_sam=5_obstruction_False.txt")
            label = f"greyscale_{grey_scale_t[j]}"
            offset = np.arange(0, 10, 0.1)
            ax[i].plot(offset[1:], np.log10(contrast), label=label, color = color_l[j], linestyle = "--")
            ax[i].set_title(f"{name[i]}")
            ax[i].set_xlabel("offset")
            ax[i].set_ylabel("contrast")
            ax[i].set_ylim(-10, 0)
            ax[i].grid(True)
            ax[i].legend()
            os.chdir("..")

    plt.tight_layout()
    plt.savefig("vortex_2_roddier_FQPM_greyscale_fmp_10_5.png")
    plt.show()





def make_lyot_diff_im(folder_name, before_lyot, after_lyot, num_noise, name, charge, zoom_pix, second_mirr=False):

    """
    plot out the images before and after lyot plane.
    folder_name -- the folder that contains the fits file used for making the image
    before_lyot -- the file name of the image before lyot plane.
    after_lyot -- the file name of the image after lyot plane.
    zoom_pix -- the zoom in size, with the unit of pix
    second_mirr -- if True the secondary mirror is added, if False there is no secondary mirror

    """
    if name == "vortex":
        name = (f"{name}_{charge}")
    limit = {"vortex_2":[1e-8, 1e-9], "vortex_4":[1e-7, 1e-8], "FQPM":[0.1,-0.1], "roddier":[1e-5, 1e-6], "vortex_6":[1e-9, 1e-10]}
    plt.style.use("seaborn-v0_8-talk")

    os.chdir(f"ideal_coro_2rd_mirror_{second_mirr}/")
    hdul_num_noise = fits.open(num_noise)
    print(os.getcwd())
    os.chdir(folder_name)

    hdul_before = fits.open(before_lyot)
    plane_bl = hdul_before[0].data[0]

    hdul_after = fits.open(after_lyot)

    if name == "roddier":
        plane_al = hdul_after[0].data[0]
    else:
        plane_al = hdul_after[0].data[0] - hdul_num_noise[0].data[0]
    shape = np.shape(zoom(plane_bl, zoom_pix))[0]

    noise_mask = make_aperture(shape, 8, obstruction=False, nsamp=1, lyot_fraction=86 / 200)
    plt.imshow(noise_mask)
    plt.show()

    # col = 3
    # col_size = 19
    # if second_mirr:
    col = 2
    col_size = 12

    fig, ax = plt.subplots(1, col, figsize=(col_size, 5))

    bl = ax[0].imshow(zoom(plane_bl, zoom_pix), vmax=0.00016,
                      extent=[-0.5 * shape, 0.5 * shape, -0.5 * shape, 0.5 * shape], aspect="auto")
    ax[0].set_ylabel("pixel")
    ax[0].set_xlabel("pixel")
    ax[0].set_title("plane before lyot stop")
    plt.colorbar(bl, ax=ax[0])

    if second_mirr:

        al = ax[1].imshow(zoom(plane_al, zoom_pix), vmax=0.00016,
                          extent=[-0.5 * shape, 0.5 * shape, -0.5 * shape, 0.5 * shape], aspect="auto")
        ax[1].set_ylabel("pixel")
        ax[1].set_xlabel("pixel")
        ax[1].set_title("plane after lyot stop")
        plt.colorbar(al, ax=ax[1])

    else:

        al_noise = ax[1].imshow(zoom(plane_al, zoom_pix), vmax = limit[name][0], vmin = limit[name][1],
                                extent=[-0.5 * shape, 0.5 * shape, -0.5 * shape, 0.5 * shape], aspect="auto")
        ax[1].set_ylabel("pixel")
        ax[1].set_xlabel("pixel")
        ax[1].set_title("plane after lyot stop")

        plt.colorbar(al_noise, ax=ax[1])

    print(np.sum(zoom(plane_al, zoom_pix)))




    plt.tight_layout()

    plt.savefig("b_a_lyot_normalized.png")
    plt.show()

    os.chdir("../..")



def star_offset(dim, obstruction, name, nsamp, fpm_sam, charge, lyot_stop, get_animation = True):

    """
    shift the phase of the entrance pupil to simulate the planet signal

    dim  -- the total number of I/D of the point of view.
    name -- the name of the focal phase mask
    nsamp -- the sampling number of the numerical simulation grid at the entrance pupil, with the unit [(numerical simulation pixel)/(I/D)]
    fpm_sam -- the sampling number of the SLM grid.
    charge -- the charge of the vortex focal phase mask. If the phase mask is not vortex then the charge will set to be None.
    lyot_stop -- the type of lyot stop, e.g. small, large or zernike.
    obstruction -- if obstruction is True, a second mirror will be added.
    get_animation -- if get_animation is True, the animation of shifting signal in the final science plane is saved.
    """

    print(os.getcwd())
    tp_s = {"small_cross": 0.904, "large_cross": 0.801, "zernike_cross": 0.922}

    os.chdir(f"ideal_coro_2rd_mirror_{obstruction}/")
    dir = f"lyot_{lyot_stop}_{name}_{charge}_sam_{nsamp}_fpm_sam={fpm_sam}_binary_False_obstruction_{obstruction}"
    os.chdir(dir)

    ims = []
    throughput = []
    fig = plt.figure()

    if nsamp < 10:
        ang = np.arange(0,10, 1/nsmap)/((2)**(0.5))
    else:
        ang = np.arange(0, 10, .1)/((2) ** (0.5))

    for i in ang:

        # A = make_aperture(dim, 8, obstruction = obstruction, obstruction_ratio=0.25)
        # if i == 0:
        #     A = A
        # else:
        #     A = A*np.exp(1j*phi_ramp(A,i,0))
        #
        # m_vor = focal_mask(name, dim, nsamp, charge)
        # lyot_stop_para = {"small_cross": [0.98, 1.2], "large_cross": [0.95, 1.3], "zernike_cross": [0.98, 1.05]}
        #
        # """ set different sizes of lyot stop."""
        # N_vor = make_lyot_stop(dim, nsamp, fraction=lyot_stop_para[lyot_stop][0],
        #                        sec_fraction=lyot_stop_para[lyot_stop][1], obstruction=obstruction)
        #
        # final_focal = coro_high_sam(A, m_vor, N_vor, nsamp, name = name, charge = charge, obstruction = obstruction, lyot_sum = lyot_sum, shift = True, save = save)
        # if i == 0:
        #     final_focal_1 = coro_high_sam(A, np.ones_like(m_vor), N_vor, nsamp, name, obstruction=True, save=False)
        #     peak = np.max(np.abs(final_focal_1) ** 2)


        # displC(zoom(final_focal, 320), f"final_focal_shift_{i}.")
        # im.set_array(zoom(np.abs(final_focal), 320))
        final_focal = simple_coro(dim, name, nsamp, fpm_sam, lyot_stop = lyot_stop, offset = i, charge = charge, binary = False, obstruction = obstruction, get_lyot = False)
        tp = np.max(final_focal)
        # print(tp)
        throughput.append(tp)

        if get_animation:
            im = plt.imshow(np.log(zoom(final_focal, 200)), vmax=-1, vmin=-5)
        # plt.show()

        # if i == 0:
        #     ax.imshow(np.log(zoom(final_focal, 200)), vmax=-1, vmin=-5)
            # cax = plt.axes((0.85, 0.1, 0.075, 0.8))
            # plt.colorbar(im, cax=cax)
            # plt.show()

            ims.append([im])

    if get_animation:
        ani = animation.ArtistAnimation(fig, ims, interval=200, blit=True, repeat_delay=1000)
        f = f"final_{name}_{charge}_obstruction={obstruction}_animation.gif"
        writergif = animation.PillowWriter(fps = 1)
        ani.save(f, writer = writergif)
    else:
        np.savetxt(f"{lyot_stop}_tp_{name}_{charge}_nsamp_{nsamp}_fpm_sam={fpm_sam}_obstruction={obstruction}.txt", np.array(throughput)*tp_s[lyot_stop])

    os.chdir("../..")


def plot_tp(lyot_stop, name, charge, obstruction, ang):
     tp = np.loadtxt(f"{lyot_stop}_tp_{name}_{charge}_obstruction={obstruction}.txt")
     ang = ang
     fig, ax = plt.subplots()
     ax.plot(ang, tp)
     plt.show()







def annuli_energy_summary(folder, obstruction = True, nsamp = 100):

    plt.style.use("seaborn-v0_8-talk")

    os.chdir(folder)

    print(os.getcwd())
    vortex_2 = np.loadtxt(f"lyot_large_cross_vortex_2_sam_100_binary_False_obstruction_{obstruction}/final_focal_plane_coro_mask=vortex_sample=100_obstruction_{obstruction}..txt")
    vortex_4 = np.loadtxt(f"lyot_large_cross_vortex_4_sam_100_binary_False_obstruction_{obstruction}/final_focal_plane_coro_mask=vortex_sample=100_obstruction_{obstruction}..txt")
    FQPM = np.loadtxt(f"lyot_large_cross_FQPM_None_sam_100_binary_False_obstruction_{obstruction}/final_focal_plane_coro_mask=FQPM_sample=100_obstruction_{obstruction}..txt")
    roddier = np.loadtxt(f"lyot_large_cross_roddier_None_sam_100_binary_False_obstruction_{obstruction}/final_focal_plane_coro_mask=roddier_sample=100_obstruction_{obstruction}..txt")

    annuli = np.arange(0, nsamp * 10, nsamp / 10)
    r_in = annuli[:-1]
    r_out = annuli[1:]
    r_mid = (r_in + r_out) / 2 / nsamp

    fig = plt.figure(figsize=(10,6))
    plt.plot(r_mid, vortex_2, label = "vortex_2")
    plt.plot(r_mid, vortex_4, label="vortex_4")
    plt.plot(r_mid, FQPM, label="FQPM")
    plt.plot(r_mid, roddier, label="roddier")
    plt.yscale("log")
    plt.xlabel(r"$\lambda/D$")
    plt.ylabel("normalized intensity")
    plt.legend()
    plt.savefig(f"final_intensity_2rd_mirror_{obstruction}.png")
    plt.show()



def final_intensity_sum(nsmap):

    plt.style.use("seaborn-v0_8-talk")

    fig, ax = plt.subplots(2, 4, figsize = (20,9))
    name_list = ["vortex_2", "vortex_4", "FQPM_None", "roddier_None"]
    obstruction_list = [False, True]
    title = ["without secondary mirror", "with secondary mirror"]

    for j in range (len(name_list)):
        for i in range(len(obstruction_list)):
            folder_name = f"ideal_coro_2rd_mirror_{obstruction_list[i]}/lyot_large_cross_{name_list[j]}_sam_{nsamp}_binary_False_obstruction_{obstruction_list[i]}"
            file_name = f"final_focal_plane_coro_mask={name_list[j]}_sample={nsamp}_obstruction_{obstruction_list[i]}.fits"
            file_name = f"{folder_name}/{file_name}"
            hdul = fits.open(file_name)
            final_int = hdul[0].data[0]
            int = ax[i][j].imshow(zoom(final_int, nsamp*10), extent= [-10,10,-10,10], aspect="auto")
            if j == 1 and i == 1:
                int = ax[i][j].imshow(zoom(final_int, nsamp * 10), vmin = 0, vmax = 0.0035, extent=[-10, 10, -10, 10], aspect="auto")

            # if i == 0:
            #     ax[i][j].set_title(title[i])
            # ax[i][j].set_xlabel(r"$\lambda$/D")
            # ax[i][j].set_ylabel(r"$\lambda$/D")
            plt.colorbar(int)

    plt.tight_layout()
    plt.savefig(f"final_intensity_sum_{nsamp}.png")
    plt.show()


def make_fpm_image_poster(fpm_sam, greyscale):

    plt.style.use("seaborn-v0_8-talk")

    plt.rcParams['axes.facecolor'] = 'black'
    plt.rcParams['figure.facecolor'] = 'black'
    plt.rcParams['text.color'] = 'white'
    plt.rcParams['axes.labelcolor'] = 'white'
    plt.rcParams['xtick.color'] = 'white'
    plt.rcParams['ytick.color'] = 'white'
    plt.rcParams['legend.facecolor'] = 'black'
    plt.rcParams['legend.edgecolor'] = 'white'

    start_fraction = 0.1  # Start of the segment
    end_fraction = 0.65  # End of the segment

    # Get the 'hot' colormap
    hot_colormap = plt.get_cmap('hot')

    # Create a new colormap from the segment of the 'hot' colormap
    new_colormap = hot_colormap(np.linspace(start_fraction, end_fraction, 256))

    # Define the new colormap
    custom_colormap = LinearSegmentedColormap.from_list('custom_hot', new_colormap)

    name_list = ["vortex_2", "vortex_4", "vortex_6", "vortex_8", "FQPM_None", "roddier_None", "dual_zone_None"]
    title_list = ["Vortex, n=2", "Vortex, n=4", "Vortex, n=6", "Vortex, n=8", "FQPM", "Roddier", "Dual Zone"]

    rc('font', weight='bold')
    fig, ax = plt.subplots(2,7, sharex= True, sharey = True, figsize = (38,10), constrained_layout=True)
    for j in range(len(fpm_sam)):
        for i in range(len(name_list)):

            filename = f"m_vor_{name_list[i]}_nsamp=100_fpm_sam={fpm_sam[j]}_rotate=False_greyscale={greyscale}.fits"

            hdul = fits.open(filename)
            if name_list[i] == "dual_zone_None":
                mask = hdul[0].data
            else:
                mask = hdul[0].data[1]
            mask = ax[j][i].imshow(zoom(mask,100), extent = (-1,1,-1,1), vmax = np.pi, vmin = -np.pi, cmap = "twilight")
            if j == 0:
                ax[j][i].set_title(title_list[i], fontsize=35, weight="bold")
            # ax[j][i].set_xlabel(r"$\lambda$/D")
            # ax[j][i].set_ylabel(r"$\lambda$/D")
            ax[j][i].xaxis.set_tick_params(labelsize=30)
            ax[j][i].yaxis.set_tick_params(labelsize=30)

    fig.supxlabel(r"$\lambda$/D", fontsize=30, weight="bold")
    # fig.supylabel(r"$\lambda$/D", fontsize=30, weight="bold")
    plt.tight_layout()
    cbar = plt.colorbar(mask, ax=ax, pad=0.01, shrink = 0.8)
    cbar.set_label("Phase [Radian]", fontsize = 30, weight = "bold")
    cbar.ax.tick_params(labelsize=30)

    plt.savefig(f"mask_poster_fmp_sam={fpm_sam}.png", transparent = True, dpi = 200)
    plt.show()

def make_fpm_image(fpm_sam, greyscale, mode="slides"):
    """
    Generate and save an FPM image poster.

    Parameters:
    - fpm_sam (list): List of FPM sampling values.
    - greyscale (bool): Whether to use greyscale images.
    - mode (str): "slides" for dark background, "paper" for white background.
    """
    if mode == "slides":
        plt.style.use("seaborn-v0_8-dark")
        bg_color = 'black'
        text_color = 'white'
    else:
        plt.style.use("seaborn-v0_8-talk")
        bg_color = 'white'
        text_color = 'black'

    plt.rcParams['axes.facecolor'] = bg_color
    plt.rcParams['figure.facecolor'] = bg_color
    plt.rcParams['text.color'] = text_color
    plt.rcParams['axes.labelcolor'] = text_color
    plt.rcParams['xtick.color'] = text_color
    plt.rcParams['ytick.color'] = text_color
    plt.rcParams['legend.facecolor'] = bg_color
    plt.rcParams['legend.edgecolor'] = text_color

    name_list = ["vortex_2", "vortex_4", "vortex_6", "vortex_8", "FQPM_None", "roddier_None", "dual_zone_None"]
    title_list = ["Vortex, n=2", "Vortex, n=4", "Vortex, n=6", "Vortex, n=8", "FQPM", "Roddier", "Dual Zone"]

    rc('font', weight='bold')
    fig, ax = plt.subplots(2, 7, sharex=True, sharey=True, figsize=(38, 10), constrained_layout=True)
    # ax = np.atleast_2d(ax)  # Ensure ax is always a 2D array

    for j in range(len(fpm_sam)):
        for i in range(len(name_list)):
            filename = f"m_vor_{name_list[i]}_nsamp=100_fpm_sam={fpm_sam[j]}_rotate=False_greyscale={greyscale}.fits"

            try:
                hdul = fits.open(filename)
                mask = hdul[0].data
                print(title_list[i])
                print(np.shape(mask))
                    # if name_list[i] == "dual_zone_None" else hdul[0].data[1]
                hdul.close()
            except FileNotFoundError:
                print(f"Warning: {filename} not found.")
                continue

            im = ax[j][i].imshow(zoom(mask, 100), extent=(-1, 1, -1, 1), vmax = 2*np.pi, vmin = 0, cmap="twilight")
            if j == 0:
                ax[j][i].set_title(title_list[i], fontsize=35, weight="bold")
            ax[j][i].xaxis.set_tick_params(labelsize=30)
            ax[j][i].yaxis.set_tick_params(labelsize=30)

    fig.supxlabel(r"$\lambda$/D", fontsize=30, weight="bold")

    cbar = plt.colorbar(im, ax=ax, pad=0.01, shrink=0.8)
    cbar.set_label("Phase [Radian]", fontsize=30, weight="bold")
    cbar.ax.tick_params(labelsize=30)

    plt.savefig(f"mask_poster_fpm_sam={fpm_sam}.png", transparent=(mode == "slides"), dpi=200)
    plt.show()


def make_lyot_plane_profile(name, charge, greyscale):
    dir = "ideal_coro_2rd_mirror_False"
    os.chdir(dir)
    fpm_sam_list = [5,10,100]

    plt.figure(figsize=(10, 10))
    for i in range(len(fpm_sam_list)):

        if greyscale == 64:
            folder = f"lyot_None_{name}_{charge}_sam_100_fpm_sam={fpm_sam_list[i]}_binary_False_obstruction_False"
        else:
            folder = f"lyot_None_{name}_{charge}_sam_100_fpm_sam={fpm_sam_list[i]}_binary_False_obstruction_False_greyscale_{greyscale}"

        filename = f"lyot_plane_profile_100_{charge}_{fpm_sam_list[i]}.txt"
        filename = f"{folder}/{filename}"
        print(filename)
        profile = np.loadtxt(filename)


        plt.plot(profile, label = f"fpm_sam = {fpm_sam_list[i]}")

    plt.xlabel("Pixel")
    plt.ylabel("Intensity Ratio")
    # plt.ylim(0.0001216, 0.0001220)
    plt.yscale("log")
    plt.legend()
    plt.show()


def make_throughtput_summary(name, charge, fpm_sam, greyscale, contrast_ind, lyot_stop, obstruction = False, same_fpm = False):

    plt.style.use("seaborn-v0_8-talk")
    dir = f"ideal_coro_2rd_mirror_{obstruction}"
    offset = np.arange(0,6,0.2)
    rc('font', weight='bold', size = 30)
    color_list = ["red", "darkgreen", "darkorange"]
    title_list = ["Vortex n=2", "Vortex, n=4", "Vortex, n=6", "FQPM", "Roddier"]
    lyot_stop_label = [1,1.2,1.3]

    if same_fpm:
        fig, ax = plt.subplots(1, 5, figsize=(25, 8), sharex = True, sharey = True)

        os.chdir(dir)
        for n in range(len(lyot_stop)):
            for j in range(len(name)):
                folder_name = f"lyot_{lyot_stop[n]}_{name[j]}_{charge[j]}_sam_100_fpm_sam={fpm_sam[0]}_binary_False_obstruction_{obstruction}_greyscale_{greyscale}"
                os.chdir(folder_name)
                contrast = np.loadtxt(
                    f"final_focal_plane_coro_mask={name[j]}_{charge[j]}_sample=100_fpm_sam={fpm_sam[0]}_obstruction_{obstruction}.txt")
                contrast = contrast[contrast_ind]
                os.chdir("..")

                tp = np.loadtxt(
                    f"lyot_{lyot_stop[n]}_{name[j]}_{charge[j]}_sam_100_fpm_sam={fpm_sam[0]}_binary_False_obstruction_{obstruction}_greyscale_{greyscale}.txt")
                SNR = tp / np.sqrt(contrast)

                if lyot_stop[n] == "None":
                    lyot_stop_l = "1_1.2"
                else:
                    lyot_stop_l = lyot_stop[n][5:]
                ax[j].plot(offset, SNR, label=rf"D$s'$/Ds= {lyot_stop_label[n]}", color = color_list[n], linewidth = 6)
                ax[j].set_title(f"{title_list[j]}", weight = "bold", fontsize = 40)
                ax[j].set_yscale("log")
                ax[j].grid(True, color = "black")

                ax[j].set_ylim(1e1, 1e4)
                ax[j].xaxis.set_tick_params(labelsize=30)
                ax[j].yaxis.set_tick_params(labelsize=30)

        ax[3].legend(loc = "best", fontsize = 25)
        ax[2].set_xlabel(r"offset[$\lambda$/D]", weight = "bold", fontsize = 30)
        ax[0].set_ylabel(r"$\eta_{p}$/$\sqrt{\eta_{s}}$", x = 0.01, weight = "bold", fontsize = 30)
        plt.tight_layout()
        plt.savefig(f"SNR_newplot_{lyot_stop}.png", transparent = True)
        plt.show()

    else:
        fig, ax = plt.subplots(1,2, figsize = (10,5))

        os.chdir(dir)
        for j in range(len(name)):
            for i in range(2):

                folder_name = f"lyot_{lyot_stop}_{name[j]}_{charge[j]}_sam_100_fpm_sam={fpm_sam[i]}_binary_False_obstruction_{obstruction}_greyscale_{greyscale}"
                os.chdir(folder_name)
                contrast = np.loadtxt(f"final_focal_plane_coro_mask={name[j]}_{charge[j]}_sample=100_fpm_sam={fpm_sam[i]}_obstruction_{obstruction}.txt")
                contrast = contrast[contrast_ind]
                os.chdir("..")

                tp = np.loadtxt(f"lyot_{lyot_stop}_{name[j]}_{charge[j]}_sam_100_fpm_sam={fpm_sam[i]}_binary_False_obstruction_{obstruction}_greyscale_{greyscale}.txt")
                SNR = tp/np.sqrt(contrast)

                ax[i].plot(offset, SNR, label = f"fpm = {name[j]}_{charge[j]}")
                ax[i].set_title(f"fpm_sam = {fpm_sam[i]}")
                ax[i].set_yscale("log")
                ax[i].grid(True)
                ax[i].set_xlabel(r"offset[$\lambda$/D]", fontsize = 30)
                ax[i].set_ylabel("SNR", fontsize = 30)
                ax[i].set_ylim(1e-2,1e4)
                ax[i].xaxis.set_tick_params(labelsize=30)
                ax[i].yaxis.set_tick_params(labelsize=30)
                if j == 3:
                    ax[i].legend()

        # fig.legend(loc = "right")
        plt.tight_layout()
        plt.savefig(f"SNR_newplot_{lyot_stop}.png")
        plt.show()

def vis_grey_scale(method):

    rc('font', weight='bold')

    if method == "1":
        array1 = np.arange(0,1,0.01)
        print(array1)
        bit_list = [2,4,8]
        for i in range(len(bit_list)):

            bin = np.linspace(0, 1, 2 ** bit_list[i], endpoint=True)

            inds = np.digitize(array1, bin)
            array2 = bin[inds]
            plt.plot(array1, array2, label = "{bit_list[i]}_bits", linewidth = 4)

        plt.legend()
        plt.show()

    if method == "2":

        grey_list = ["2","4","6","8","64"]
        for i in grey_list:

            filename = f"m_vor_vortex_8_nsamp=100_fpm_sam=100_rotate=False_greyscale={i}.fits"
            hdul = fits.open(filename)
            mask = zoom(hdul[0].data[1],50)

            plt.plot(mask[90], label = f"{i} bit", linewidth = 3)
            plt.xlabel("Pixel", weight = "bold")
            plt.ylabel("Phase Shift Value", weight = "bold")
            plt.xlim(0,5)
            plt.ylim(0,1)

        plt.legend()
        plt.savefig(f"vis_greyscale_{grey_list}.png")
        plt.show()


    if method == "3":

        fig, ax = plt.subplots(figsize = (15,15))
        rc("font", size = 30)
        grey_list = ["4","6","8","64"]
        for i in grey_list:
            filename = f"m_vor_vortex_8_nsamp=100_fpm_sam=100_rotate=False_greyscale={i}.fits"
            hdul = fits.open(filename)
            phase_mask = zoom(hdul[0].data[1], 300)


            rr = r_theta(phase_mask, (np.shape(phase_mask)[0] - 1) / 2, (np.shape(phase_mask)[0] - 1) / 2)[0]
            mask = (rr > 200) * (rr < 201)
            mask_1 = np.zeros_like(mask)
            mask_1[:int(np.shape(mask)[0]/2),int(np.shape(mask)[0]/2):] = 1
            mask_1 = mask*mask_1
            mask_2 = mask*mask_1*(-4)
            # ax.imshow((phase_mask*(~mask_1)+mask_2)[0:190,:], vmax = np.pi)
            # ax[0].set_ylim(190,0)


            phase_mask = phase_mask[mask_1]
            ax.plot(phase_mask*360/np.pi, label=f"{i} bits", linewidth=5)
            ax.set_xlabel("Pixel", weight="bold", fontsize = 30)
            ax.set_ylabel("Phase Shift Value [Degree]", weight="bold", fontsize = 30)
            ax.legend()
            ax.xaxis.set_tick_params(labelsize=30)
            ax.yaxis.set_tick_params(labelsize=30)
            plt.xlim(16,70)
            plt.ylim(70,360)


        fig.subplots_adjust(wspace=0.02)
        plt.tight_layout()
        plt.savefig(f"vis_greyscale_{grey_list}_method_{method}_zoom.png", transparent = True)
        plt.show()


def displC(c, title, only_phase = False, trim=0):
    """displC - display a Complex number c as four plots

               The top two plots are (Real, Imaginary) quantities
               The bottom two plots are (Amplitude, Phase)

               Optionally cut out the central square  with a size of 'trim x trim' pixels"""
    c2 = np.copy(c)
    if (trim > 0):  # if the user specifies a trim value, cut out the centre of the image
        (nx, ny) = c.shape
        dx = int((nx - trim) / 2)
        dy = int((nx - trim) / 2)
        c2 = c[dx:dx + trim, dy:dy + trim]

    # set up the plot panels
    fig = plt.figure(figsize=(10, 8))
    axre = fig.add_subplot(221)
    axim = fig.add_subplot(222)
    axamp = fig.add_subplot(223)
    axpha = fig.add_subplot(224)
    # plot out the panels
    im = axre.imshow(c2.real)
    im = axim.imshow(c2.imag)
    im = axamp.imshow(np.abs(c2))
    im1 = axpha.imshow(np.angle(c2))

    axre.set_title('Real')
    axim.set_title('Imag')
    axamp.set_title('Amplitude')
    axpha.set_title('Phase')

    if only_phase:
        cube = np.array([np.angle(c2)])
    else:
        cube = np.array([np.abs(c2),np.angle(c2)])
    wfits(cube, title + "fits")


    # fig.subplots_adjust(right=0.8)
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    # fig.colorbar(im1, cax=cbar_ax)
    # plt.savefig(title + "png")
    plt.show()


def fpm_diff():
    # Load the first FITS file
    file1 = 'ideal_coro_2rd_mirror_True/lyot_None_1_1_vortex_2_lambda_1_fpm_sam=10_binary_False_obstruction_True_greyscale_8/m_vor_vortex_2_nsamp=10_fpm_sam=10_rotate=False_greyscale=8.fits'
    hdul1 = fits.open(file1)
    data1 = hdul1[0].data

    # Load the second FITS file
    file2 = 'ideal_coro_2rd_mirror_True/lyot_None_1_1_vortex_2_lambda_0.9_fpm_sam=10_binary_False_obstruction_True_greyscale_8/m_vor_vortex_2_nsamp=10_fpm_sam=10_rotate=False_greyscale=8.fits'
    hdul2 = fits.open(file2)
    data2 = hdul2[0].data

    # Compute the difference
    difference = data2 - data1

    # Close the FITS files
    hdul1.close()
    hdul2.close()

    # Plot the difference
    plt.figure(figsize=(10, 8))
    plt.imshow(difference)
    plt.colorbar()
    plt.title(r'Difference between vortex 2, 0.9 $\lambda_{0}$ and $\lambda_{0}$')
    plt.xlabel('X Pixel')
    plt.ylabel('Y Pixel')

    # Save the figure as a PNG file
    plt.savefig('difference_between_vortex_2_fpm_0.9_1_files.png')

    # Show the plot
    plt.show()


def plot_psf_plot(name, charge, lyot = False):

    folder_name = "ideal_coro_2rd_mirror_True"
    os.chdir(folder_name)
    # Array of folder names
    folder_names = ["None_1_1", "None_1_1.1", "None", "None_1_1.3", "None_1_1.4"]
    title_list = ["D's/Ds = 1", "D's/Ds = 1.1", "D's/Ds = 1.2", "D's/Ds = 1.3", "D's/Ds = 1.4"]

    # Base path for the folder
    base_path = "lyot_{0}_{1}_{2}_lambda_1_fpm_sam=10_binary_False_obstruction_True_greyscale_8"

    # File name within each folder
    if lyot:
        print("1")
        file_name = "after_lyot_plane_profile_10_{0}_10.fits"
    else:
        file_name = "final_focal_plane_coro_mask={0}_{1}_sample=10_fpm_sam=10_obstruction_True.fits"


    # Initialize variables to track the global min and max values
    global_min = np.inf
    global_max = -np.inf

    # First pass to determine global min and max values
    for folder in folder_names:
        # Construct the full path to the FITS file
        if lyot:
            full_path = os.path.join(base_path.format(folder, name, charge), file_name.format(charge))
        else:
            full_path = os.path.join(base_path.format(folder, name, charge), file_name.format(name, charge))

        # Read the FITS file
        with fits.open(full_path) as hdul:
            if lyot:
                data = hdul[0].data
            else:
                data = np.log10(hdul[0].data[0])[450:550, 450:550]
        # Update global min and max

        global_min = min(global_min, np.min(data))
        global_max = max(global_max, np.max(data))


    # Create a figure with 5 subplots
    fig, axes = plt.subplots(1, 5, figsize=(25, 5))

    # Second pass to plot the data with consistent colorscale
    for i, folder in enumerate(folder_names):
        # Construct the full path to the FITS file
        if lyot:
            full_path = os.path.join(base_path.format(folder, name, charge), file_name.format(charge))
        else:
            full_path = os.path.join(base_path.format(folder, name, charge), file_name.format(name, charge))

        # Read the FITS file
        with fits.open(full_path) as hdul:
            if lyot:
                data = hdul[0].data
            else:
                data = np.log10(hdul[0].data[0])[450:550, 450:550]

        # Plot the data
        ax = axes[i]
        if lyot:
            im = ax.imshow(data, origin='lower', vmax=global_max, vmin=global_min)
        else:
            im = ax.imshow(data, origin='lower', vmax = -2, vmin = -7)
        ax.set_title(title_list[i])

        # Set tick labels to "spex" where 1 spex = 10 pixels
        if lyot == False:
            num_ticks = data.shape[0] // 10
            ticks = np.arange(0, num_ticks * 10, 10)
            tick_labels = np.arange(-(num_ticks // 2), num_ticks // 2 )

            ax.set_xticks(ticks)
            ax.set_xticklabels(tick_labels)
            ax.set_yticks(ticks)
            ax.set_yticklabels(tick_labels)

            ax.set_xlabel(r'$\lambda$/D')
            ax.set_ylabel(r'$\lambda$/D')

    cbar = fig.colorbar(im, ax=axes[4], orientation='vertical', fraction=0.02, pad=0.04)
    cbar.set_label('Intensity')

    # Adjust layout and show the plot
    plt.tight_layout()
    if lyot:
        plt.savefig(f"after_lyot_{name}_{charge}.png")
    else:
        plt.savefig(f"final_psf_{name}_{charge}.png")
    plt.show()

def plot_tp(name, charge):

    data = f"ideal_coro_2rd_mirror_True/lyot_None_1_1_{name}_{charge}_sam_10_fpm_sam=10_binary_False_obstruction_True_greyscale_8_bb.txt"
    plt.plot(data)
    plt.show()




































