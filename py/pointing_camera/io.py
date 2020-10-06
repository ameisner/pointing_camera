import astropy.io.fits as fits
import os
import fitsio
import pointing_camera.common as common

def write_image_level_outputs(exp, outdir):

    print('Attempting to write image level outputs')

    assert(os.path.exists(outdir))

    outname = (os.path.split(exp.fname_im))[-1]

    outname = outname.replace('.fits', '-detrended.fits')

    outname = os.path.join(outdir, outname)
    
    outname_tmp = outname + '.tmp'

    assert(not os.path.exists(outname))
    assert(not os.path.exists(outname_tmp))
    
    fitsio.write(outname_tmp, exp.detrended.astype('float32'))
    os.rename(outname_tmp, outname)

def write_sky_summary(sky, exp, outdir):

    print('Attempting to write sky summary table')

    assert(os.path.exists(outdir))

    outname = (os.path.split(exp.fname_im))[-1]

    outname = outname.replace('.fits', '-sky.fits')

    outname = os.path.join(outdir, outname)
    
    outname_tmp = outname + '.tmp'

    assert(not os.path.exists(outname))
    assert(not os.path.exists(outname_tmp))
    
    sky.write(outname_tmp, format='fits')
    os.rename(outname_tmp, outname)

def load_static_badpix():
    par = common.pc_params()

    fname = os.path.join(os.environ[par['meta_env_var']],
                         par['static_mask_filename'])

    assert(os.path.exists(fname))

    mask = fits.getdata(fname)

    return mask
