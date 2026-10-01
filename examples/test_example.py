import numpy as np
import matplotlib.pyplot as plt

from active_coronagraph.main_functions_polished import (
    broadband_combine_final_psf,
    run,
)

wavelengths = [0.9, 0.95, 1, 1.05, 1.10]
weights = np.ones(len(wavelengths))
weights[[0, -1]] = 0.5

for wl in wavelengths:
    run(
        dim=100,
        name="vortex",
        wavelength=wl,
        nsamp=10,
        fpm_sam=10,
        charge=2,
        lyot_stop="None_1_1",
        binary=False,
        obstruction=False,
        get_lyot=False,
        greyscale=8,
        folder_name="_broadband_example",
    )

broadband_psf = broadband_combine_final_psf(
    wavelength=wavelengths,
    fpm_sam=10,
    nsamp=10,
    weights=weights,
    name="vortex",
    charge=2,
    obstruction=False,
    lyot_stop="None_1_1",
    folder_name="_broadband_example",
)

image = np.squeeze(broadband_psf)
plot_data = np.log10(np.clip(image, 1e-12, None))

fig, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(plot_data, origin="lower", cmap="inferno")
ax.set_title("20% broadband vortex charge 2 normalized intensity")
ax.set_xlabel("x [pixel]")
ax.set_ylabel("y [pixel]")
cbar = fig.colorbar(im, ax=ax)
cbar.set_label("log10 normalized intensity")
fig.tight_layout()
fig.savefig("test_example_broadband_vortex_image.png", dpi=180)
plt.close(fig)
