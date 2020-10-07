#!/usr/bin/env python

import argparse
from datetime import datetime
import os
import time
from exposure import PC_exposure
import pointing_camera.util as util
import pointing_camera.io as io

def pc_proc(fname_in, outdir=None):

    print('Starting pointing camera reduction pipeline at: ' +
          str(datetime.utcnow()) + ' UTC')

    t0 = time.time()

    write_outputs = (outdir is not None)
    
    try:
        print('Running on host: ' + str(os.environ.get('HOSTNAME')))
    except:
        print('Could not retrieve hostname!')

    assert(os.path.exists(fname_in))

    exp = PC_exposure(fname_in)

    util.detrend_pc(exp)

    sky = util.sky_summary_table(exp)

    cat = util.pc_phot(exp)

    if write_outputs:
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        io.write_image_level_outputs(exp, outdir)
        io.write_sky_summary(sky, exp, outdir)
        io.write_source_catalog(cat, exp, outdir)

    dt = time.time() - t0
    print('pointing camera reduction pipeline took ' + '{:.2f}'.format(dt) +
          ' seconds')
    print('pointing camera reduction pipeline completed at: ' +
          str(datetime.utcnow()) + ' UTC') 

if __name__ == "__main__":
    descr = 'run the pointing camera reduction pipeline on an exposure'

    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument('fname_in', type=str, nargs=1,
                        help='pointing camera raw image file name')

    parser.add_argument('--outdir', default=None, type=str,
                        help='directory to write outputs in')
    
    args = parser.parse_args()
    
    pc_proc(args.fname_in[0], outdir=args.outdir)
