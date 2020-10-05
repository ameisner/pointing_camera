import astropy.io.fits as fits
import os
import fitsio

def write_image_level_outputs(exp, outdir):

    print('Attempting to write image level outputs')

    assert(os.path.exists(outdir))

    outname = (os.path.split(exp.fname_im))[-1]

    outname = outname.replace('.fits', '-detrended.fits')

    outname = os.path.join(outdir, outname)
    
    outname_tmp = outname + '.tmp'

    fitsio.write(outname_tmp, exp.detrended)
    os.rename(outname_tmp, outname)

    
