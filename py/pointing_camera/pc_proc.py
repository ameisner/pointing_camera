#!/usr/bin/env python

import argparse
from datetime import datetime
import os
import time
from exposure import PC_exposure

def pc_proc(fname_in):

    print('Starting pointing camera reduction pipeline at: ' +
          str(datetime.utcnow()) + ' UTC')

    t0 = time.time()

    try:
        print('Running on host: ' + str(os.environ.get('HOSTNAME')))
    except:
        print('Could not retrieve hostname!')

    assert(os.path.exists(fname_in))

    exp = PC_exposure(fname_in)

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

    args = parser.parse_args()
    
    pc_proc(args.fname_in[0])
