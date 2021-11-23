#!/usr/bin/env python

import argparse
import glob
import time
import os
import pointing_camera.pc_proc as pipeline
import numpy as np

default_data_dir = '/global/cfs/cdirs/desi/users/ameisner/pointing_camera/el_nino'

default_outdir = '/global/cfs/cdirs/desi/users/ameisner/pointing_camera/summer_2021/v0006'

def _proc_many_exp(flist, indstart, nproc, outdir=default_outdir):
    n_all = len(flist)

    assert(indstart < n_all)

    indend = min(indstart + nproc, n_all)

    for i in range(indstart, indend):
        f = flist[i]

        assert(os.path.exists(f))

        try:
            pipeline.pc_proc(f, outdir=outdir, dont_write_detrended=True,
                             skip_checkplot=False, nightly_subdir=True,
                             send_redis=False, one_aper=True)
        except:
            print('PROCESSING FAILURE: ' + f)

def _load_flist():
    fname = '/global/homes/a/ameisner/pointing_camera/bin/summer_2021_wcs_solutions.txt'

    flist = np.genfromtxt(fname, dtype='U')

    flist = [f.replace('.wcs', '.fits') for f in flist]

    seed = 99
    np.random.seed(seed)

    np.random.shuffle(flist)

    return flist

if __name__ == "__main__":
    descr = 'process a chunk of the summer 2021 El Nino data set'

    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('indstart', type=int, nargs=1,
                        help="starting index within list of files")

    parser.add_argument('nproc', type=int, nargs=1,
                        help="number of files to process")

    parser.add_argument('--outdir', type=str, default=default_outdir,
                        help="base output directory")

    args = parser.parse_args()

    indstart = args.indstart[0]
    nproc = args.nproc[0]

    flist = _load_flist()

    _proc_many_exp(flist, indstart, nproc, outdir=args.outdir)
