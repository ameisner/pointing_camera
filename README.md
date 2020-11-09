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

# full help for running the pipeline

    pointing_camera/py/pointing_camera> python pc_proc.py --help
    usage: pc_proc.py [-h] [--outdir OUTDIR] fname_in

    run the pointing camera reduction pipeline on an exposure

    positional arguments:
      fname_in         pointing camera raw image file name

    optional arguments:
      -h, --help       show this help message and exit
      --outdir OUTDIR  directory to write outputs in
      --dont_write_detrended
                            don't write detrended image
      --skip_checkplot      don't create a checkplot
      --nightly_subdir      create output subdirectories per observing night
