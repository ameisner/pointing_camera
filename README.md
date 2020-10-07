This pointing camera reduction/analysis pipeline is implemented in pure Python. It has a relatively small number of dependencies that must be installed (see requirements.txt).

It also requires Gaia catalog files accessed via an environment variable called PC_GAIA_DIR and a static bad pixel mask accessed via an environment variable called POINTING_CAMERA_META that points to a directory for ancillary calibration products.

Here is an example invocation for running the full pipeline:

    python -u pc_proc.py /global/cfs/cdirs/desi/users/ameisner/pointing_camera/nino/20191103.234228.00498_03000.fits --outdir 20191103

This produces the following output files:

    20191103/20191103.234228.00498_03000-catalog.fits
    20191103/20191103.234228.00498_03000-sky.fits
    20191103/20191103.234228.00498_03000-zp.png
    20191103/20191103.234228.00498_03000-detrended.fits
    20191103/20191103.234228.00498_03000-zeropoints.fits

* The -catalog output is a source catalog with centroids and fluxes.
* The -sky output is a summary table of sky brightness measurements.
* The -zp output is a checkplot of instrumental magnitudes versus color-corrected Gaia G.
* The -detrended output is a detrended version of the raw pointing camera image.
* The -zeropoints output is a summary table of zeropoints (for several apertures and different image regions).
