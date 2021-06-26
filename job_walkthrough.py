
length = 200
kbt = 1
amplitude = -3
tension = 20
num_polymers = 10
runtime = pow(10,7)
sample_freq = pow(10,4)
sheer_start = 0
sheer_end = 3
DF = 10
ID = 10000
disorder_level = 0.1

import os
os.system(f"python submit_scattered_jobs.py {length} {kbt} {amplitude} {tension} {num_polymers} {runtime} {sample_freq} {sheer_start} {sheer_end} {DF} {ID} {disorder_level}")