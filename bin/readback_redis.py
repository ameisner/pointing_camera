
import redis
import os
import time

host = os.environ['REDISHOST']
port = int(os.environ['REDISPORT'])
db = int(os.environ['REDISDBNUM'])
key = os.environ['REDISKEY']

r = redis.Redis(host=host, port=port, db=db)


delay_seconds = 5

for i in range(300):
    # try hgetall later
    time.sleep(delay_seconds)
    data = r.hgetall(key)

    print(data[b'timestamp'].decode(), float(data[b'zp_adu_per_s']), float(data[b'sky_adu_per_s']), float(data[b'mjd_obs']), i)

