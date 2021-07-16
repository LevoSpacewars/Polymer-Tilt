#here we are going to submit parallel jobs with group ID
#first give name to all the inputs
#then assign sysargs to a dictionary


args = ["length", "kbt", "amplitude", "tension", "npolymers", "runtime", "samplerate", "sheer_start","sheer_end","DF","ID","disorder_level","ModalDisorder"]
dic = {}
import sys
import os
import random
for i in range(1,len(sys.argv)):
    dic.update({args[i-1]: sys.argv[i], (i-1):sys.argv[i]})


# Now calculate sheer array
sheers = []
tstart = float(dic["sheer_start"])
tend = float(dic["sheer_end"])
dif = tend - tstart
conv = dif / float(dic["DF"])
rseed = random.randint(1,1000000)
for i in range(int(dic["DF"])):
    sheers.append(conv * i + tstart)

print("probing sheers: " + str(sheers))

#now let's get all the values in order
argsouput = []
for j in range(int(dic["DF"])):
    sargs = ""
    for i in range(7):
        sargs += " " + dic[i]
    sargs += " " + str(sheers[j]) + " " + str(tend) + " " + dic["ID"] + " " + str(rseed) + " " + str(dic['disorder_level'] + " " + str(dic['ModalDisorder']))

    argsouput.append(sargs)
print(argsouput)
#here we need to create the individual batch files that are going to call the Engine_fast_runs.py [args]

settings = "#!/bin/bash \n#SBATCH --job-name=Poltmer_Run_scattered     ### Job Name \n#SBATCH --partition=gpu	  ### Quality of Service (like a queue in PBS) \n#SBATCH --time=0-10:00:00     ### Wall clock time limit in Days-HH:MM:SS \n#SBATCH --nodes=1             ### Node count required for the job \n#SBATCH --ntasks-per-node=1   ### Nuber of tasks to be launched per Node \n#SBATCH --gres=gpu:1          ### General REServation of gpu:number of gpus \n#SBATCH --account=softmatter    ### Account used for job submission \n"
cmd = "python3 /home/apatapof/softmatter/Polymer-Tilt/Engine_fast_runs.py "
for i in range(int(dic["DF"])):
    filename = dic["ID"] + "_" + str(sheers[i]) + ".batch"
    print("creating: " + filename)
    cmda = cmd + argsouput[i] + " $SLURM_JOB_GPUS"
    f = open(filename, "w")
    f.write(settings)
    f.write(cmda)
    f.close()
    print("submitting job for: " + filename)
    os.system("sbatch " + filename)
