import astropy.io.fits as fits
import os
import pointing_camera.common as common
import matplotlib.pyplot as plt

def write_image_level_outputs(exp, outdir):

    print('Attempting to write image level outputs')

    assert(os.path.exists(outdir))

    outname = (os.path.split(exp.fname_im))[-1]

    outname = outname.replace('.fits', '-detrended.fits')

    outname = os.path.join(outdir, outname)
    
    outname_tmp = outname + '.tmp'

    assert(not os.path.exists(outname))
    assert(not os.path.exists(outname_tmp))

    hdu = fits.PrimaryHDU(exp.detrended.astype('float32'))
    hdu.writeto(outname_tmp)
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

def write_source_catalog(cat, exp, outdir):

    print('Attempting to write source catalog')

    assert(os.path.exists(outdir))

    outname = (os.path.split(exp.fname_im))[-1]

    outname = outname.replace('.fits', '-catalog.fits')

    outname = os.path.join(outdir, outname)

    outname_tmp = outname + '.tmp'

    assert(not os.path.exists(outname))
    assert(not os.path.exists(outname_tmp))

    cat.write(outname_tmp, format='fits')
    os.rename(outname_tmp, outname)

def write_zeropoints_table(zps, exp, outdir):

    print('Attempting to write zeropoints table')

    assert(os.path.exists(outdir))

    outname = (os.path.split(exp.fname_im))[-1]

    outname = outname.replace('.fits', '-zeropoints.fits')

    outname = os.path.join(outdir, outname)

    outname_tmp = outname + '.tmp'

    assert(not os.path.exists(outname))
    assert(not os.path.exists(outname_tmp))

    zps.write(outname_tmp, format='fits')
    os.rename(outname_tmp, outname)

def load_static_badpix():
    par = common.pc_params()

    fname = os.path.join(os.environ[par['meta_env_var']],
                         par['static_mask_filename'])

    assert(os.path.exists(fname))

    mask = fits.getdata(fname)

    return mask

def save_zp_checkplot(exp, outdir):

    print('Attempting to save zeropoint check plot')

    if not plt.fignum_exists(1):
        return

    assert(os.path.exists(outdir))

    outname = (os.path.split(exp.fname_im))[-1]

    outname = outname.replace('.fits', '-zp.png')
    outname_tmp = 'tmp.' + outname

    outname = os.path.join(outdir, outname)
    outname_tmp = os.path.join(outdir, outname_tmp)

    assert(not os.path.exists(outname))
    assert(not os.path.exists(outname_tmp))

    plt.savefig(outname_tmp, dpi=200, bbox_inches='tight')
    os.rename(outname_tmp, outname)
