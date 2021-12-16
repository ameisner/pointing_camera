import DOSlib.exposure as exp
import psycopg2 as psycopg
import psycopg2.extras

sql = 'SELECT id, mjd_obs, exptime, reqra, reqdec, skyra, skydec, targtra, targtdec FROM exposure WHERE exposure.night = 20210611'

conn =  psycopg.connect(exp.dsn)

cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

cursor.execute(sql)

data = cursor.fetchall()

cursor.close()

conn.close()
