#compileAllRuns or compile a specific file
import os
import sys
SIMPATH = "/home/apatapof/softmatter/runs/"
path = SIMPATH
settings = "#!/bin/bash\n#SBATCH --no-requeue   ### Partition (like a queue in PBS)\n#SBATCH --job-name=CompileData      ### Job Name\n#SBATCH --output=Hi.out         ### File in which to store job output\n#SBATCH --error=Hi.err          ### File in which to store job error messages\n#SBATCH --time=0-01:00:00	### Wall clock time limit in Days-HH:MM:SS\n#SBATCH --nodes=1               ### Number of nodes needed for the job\n#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node\n#SBATCH --mem=40000M              ### memory limit per node, in MB\n#SBATCH --account=softmatter	  ### Account used for job submission \n"
command = "/home/apatapof/softmatter/Polymer-Tilt/DataHandler/Executable /home/apatapof/softmatter/runs/"

dirs =[]
print("starting")
dir_content = os.listdir(path)
if len(sys.argv) > 1:
    for item in dir_content:
        if os.path.isdir(path + item) and os.path.exists(path + item + "/trajectory.gsd") and sys.argv[1] in item:
            dirs.append(item)
        print(path  + item + "/trajectory.gsd")
else:
    for item in dir_content:
        if os.path.isdir(path + item) and os.path.exists(path + item + "/trajectory.gsd"):
            dirs.append(item)
        print(path  + item + "/trajectory.gsd")
print(dirs)
for items in dirs:
    name = "compile_" + str(items) + ".batch"
    f = open("compile_" + str(items) + ".batch",'a')
    f.write(settings)
    cmd = command + str(items)
    print(cmd)
    f.write(cmd)
    f.close()
    os.system("sbatch " + name)
    print(items)


