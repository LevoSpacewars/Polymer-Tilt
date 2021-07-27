import sys
length = 200                        # number of polymers that
kbt = 1                             # energy of the system
amplitude = -3                      # amplitude of the potential:: SEE README
tension = 20                        # tension force that pulls on polymers vertically (only along top particle since bottom is held fixed along vertical)
num_polymers = 10                   # number of polymers that populate a system
runtime = pow(10,7)                 # amount of steps to be taken for entire simulation
sample_freq = pow(10,4)             # frequency which simulation records and writes positioal data to trajectory.gsd
sheer_start = 0                     # starting sheerforce value
sheer_end = 3                       # ending sheerforce value
DF = 10                             # amount of sheer force steps inbetween sheer_start and sheer_end
ID = 10000                          # MUST BE UNIQUE:: might be useful to implement something to auto assign id
disorder_level = 0                  # indepenedent of modal or nonmodal disorder: 0-1 for disorder <= Periodic Amplitude;
ModalDisorder = 0                   # this should always either be 0 or 1 (true or false)
numbersims = 1
startID = ID
if len(sys.argv) > 1:
    for index in range(len(sys.argv)):
        if index == 1:
            numbersims = int(sys.argv[index])
        elif index == 2:
            disorder_level = float(sys.argv[index])
        elif index == 3:
            ModalDisorder = int(sys.argv[index])
        elif index == 4:
            startID = int(sys.argv[index])


import os
s = " "
for i in range(numbersims):
    ID = startID + i
    cmd = "python submit_scattered_jobs.py {} {} {} {} {} {} {} {} {} {} {} {} {}".format(length, kbt, amplitude, tension, num_polymers, runtime, sample_freq, sheer_start, sheer_end, DF, ID, disorder_level, ModalDisorder)
    print(cmd)
    os.system(cmd)
