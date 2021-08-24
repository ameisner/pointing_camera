#!/usr/bin/env python

import argparse
import astropy.io.fits as fits
import numpy as np
import os
import glob
from astropy.table import Table
import imageio
from pointing_camera.exposure import PC_exposure
import pointing_camera.util as util
from scipy.stats import scoreatpercentile
import pointing_camera.common as common

desidir = '/global/cfs/cdirs/desi/spectro/data'
lostfound = '/global/cfs/cdirs/desi/spectro/staging/lost+found'
pcdir = '/global/cfs/cdirs/desi/users/ameisner/pointing_camera/el_nino'

def rebin(a, shape):
    # https://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def guide_cube_flist(night):
    '''Get list of guide cubes'''

    flist = []
    for basedir in [desidir, lostfound]:
        pattern = os.path.join(basedir, night, '?'*8, 'guide-????????.fits.fz')
        flist += glob.glob(pattern)

    # any reason to sort the file list?

    # try to be careful about possible duplicates in both spectro/data and
    # staging/lost+found?

    return flist

def guide_cube_mjd_ranges(night):
    '''get a tabulation of guide cube file names and corresponding MJD ranges'''
    # check both spectro/data and staging/lost+found

    flist = guide_cube_flist(night)

    if len(flist) == 0:
        return None

    extnames = ['GUIDE' + str(gfa_loc) + 'T' for gfa_loc in [0, 2, 3, 5, 7, 8]]

    results = []
    for i, f in enumerate(flist):
        print('working on ', i+1, ' of ', len(flist), ' : ' + f)
        try:
            hdul = fits.open(f)
        except:
            print('corrupt raw guide cube??' + f)
            continue

        for extname in extnames:
            if extname in hdul:
                tab = hdul[extname].data
                # discard acquisition image ; could crash if guide 'cube'
                # has exactly one frame
                tab = tab[1:]
                expid = (os.path.split(f)[-1]).replace('guide-', '').replace('.fits.fz', '')
                result = (f, np.min(tab['MJD-OBS']), np.max(tab['MJD-OBS']),
                          extname, expid)
                results.append(result)
                break

    t = Table()
    t['FNAME'] = [r[0] for r in results]
    t['MJDMIN'] = [r[1] for r in results]
    t['MJDMAX'] = [r[2] for r in results]
    t['EXTNAME'] = [r[3] for r in results]
    t['EXPID'] = [r[4] for r in results]

    return t

def pointing_camera_index(night):
    '''get a tabulation of pointing camera exposures and their MJD values'''

    year = night[0:4]
    month = night[4:6]
    day = night[6:8]

    # pointing camera raw data convention
    monthdir = year + '-' + month
    nightdir = year + '-' + month + '-' + day
    
    pattern = os.path.join(pcdir, monthdir, nightdir, '*.fits')
    
    flist = glob.glob(pattern)

    if len(flist) == 0:
        return None
    
    flist.sort()

    results = []
    for i, f in enumerate(flist):
        print('working on ', i+1, ' of ', len(flist), ' : ' + f)
        try:
            h = fits.getheader(f)
        except:
            print('corrupt raw pointing camera image??' + f)
            continue

        if 'MJD-OBS' in h:
            tai_utc_offs = 37.0/(24.0*3600.0) # in days
            result = (f, h['MJD-OBS'] - tai_utc_offs)
            results.append(result)
        else:
            print(f + ' does not have MJD-OBS??')

    t = Table()
    t['FNAME'] = [r[0] for r in results]
    t['MJD'] = [r[1] for r in results]

    return t

def one_pc_rendering(fname):
    '''make a low-resolution rendering of one pointing camera image'''

    # read in
    exp = PC_exposure(fname)

    # detrend
    util.detrend_pc(exp)

    # downsample via averaging
    binfac = 8

    sh = exp.detrended.shape
    im = rebin(exp.detrended, (int(sh[0]/binfac), int(sh[1]/binfac)))
    sh = im.shape
    
    # figure out the stretch (should this be held constant across frames??)
    limits = scoreatpercentile(np.ravel(im), [1, 99])

    im[im < limits[0]] = limits[0]
    im[im > limits[1]] = limits[1]

    # overplot approximate DESI FOV to guide the eye
    par = common.pc_params()
    
    ybox = np.arange(sh[0]*sh[1], dtype=int) // sh[1]
    xbox = np.arange(sh[0]*sh[1], dtype=int) % sh[1]

    xbox = xbox.astype('float')
    ybox = ybox.astype('float')

    x_center = (sh[1] // 2) - 0.5*((sh[1] % 2) == 0)
    y_center = (sh[0] // 2) - 0.5*((sh[0] % 2) == 0)

    xbox -= x_center
    ybox -= y_center

    xbox = xbox.reshape(sh)
    ybox = ybox.reshape(sh)

    dist = np.sqrt(np.power(xbox, 2) + np.power(ybox, 2))

    mask = np.abs(dist - par['desi_radius_pix']/binfac) < 0.5

    im[mask] = limits[1]

    im -= limits[0]

    return im

def desi_exp_movie(mjdrange, expid, pc_index, fps=5, outdir='.'):
    '''make an animation of pointing camera images during one DESI exposure'''

    print('Working on DESI exposure ', expid)

    keep = (pc_index['MJD'] >= mjdrange[0]) & (pc_index['MJD'] < mjdrange[1])

    if np.sum(keep) == 0:
        print('no pointing camera exposures for exposure ' + expid)
        return

    pc_index = pc_index[keep]
    pc_index.sort('MJD')

    ims = []
    for row in pc_index:
        im = one_pc_rendering(row['FNAME'])
        ims.append(im)

    # what happens iin the case of exactly one pointing camera image?

    outname = expid + '-pointing_camera.gif'
    outname = os.path.join(outdir, outname)

    # use imageio to create the GIF, specifying an FPS value
    imageio.mimsave(outname, ims, fps=10)

def movies_1night(night, outdir='.'):
    '''generate all the per DESI exposure animations for one night'''

    pc = pointing_camera_index(night)
    cubes = guide_cube_mjd_ranges(night)
    
    # if no guide cubes or no pointing camera images on requested night
    # then give up

    if (pc is None) or (cubes is None):
        print('no pointing camera animations to generate')
        return

    for cube in cubes:
        # at some point could try to take into account the
        # *ending* MJD (rather than starting MJD) of the final guide frame
        mjdrange = (cube['MJDMIN'], cube['MJDMAX'])
        desi_exp_movie(mjdrange, cube['EXPID'], pc, outdir=outdir)

if __name__ == "__main__":
    descr = 'pointing camera movies for one night worth of DESI guiding'

    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('night', type=str, nargs=1,
                        help="observing night in YYYYMMDD format")

    parser.add_argument('--outdir', type=str, default='.',
                        help="output directory")

    parser.add_argument('--nightly_subdir', default=False,
                        action='store_true',
                        help="make nightly subdirectory for GIFs")

    args = parser.parse_args()

    night = args.night[0]

    assert(len(night) == 8)

    outdir = args.outdir

    assert(os.path.exists(outdir))
    
    if args.nightly_subdir:
        outdir = os.path.join(outdir, night)
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    movies_1night(night, outdir=outdir)
