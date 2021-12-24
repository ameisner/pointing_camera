#!/usr/bin/env python

import pointing_camera.desi as desi
from pointing_camera.util import get_obs_night
import os
import time
from datetime import datetime

print('PID = ' + str(os.getpid()))

out_basedir = os.environ['PC_MOVIE_OUTDIR']
workers = 8

thresh = 700 # 7am KPNO time

last_time = 10000
while True:
    time.sleep(60)
    now = datetime.now()

    current_time = now.strftime("%H:%M:%S")
    print("Current Time = ", current_time)

    current_time = current_time.replace(':', '')
    current_time = int(current_time[0:4])

    print(last_time, current_time)

    date_string_local = now.strftime("%Y/%m/%d")
    time_string_local = now.strftime("%H:%M:%S") 

    night = get_obs_night(date_string_local, time_string_local)
    outdir = os.path.join(out_basedir, night)

    if (current_time >= thresh) and (last_time < thresh):
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        desi.all_movies_1night(night, outdir=outdir, nmp=workers)

    last_time = current_time
