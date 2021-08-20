README:
This readme is seperated into 3 parts
	1. folder layout with desc
	2. introduction to files for making figures
	3. sending jobs for making data points
	4. adding functionality to DataHanlder
--------------- 1 ---------------------

Polymer-Tilt: Main-Directory
compiledRuns: Directory with processed run data
	All processed run data is found here

DataHandler: C++ code and excecutable for analysing trajectory files
singlePolymer: Frozen out single polymer system (it's not important but could be useful for debugging)


--------------- 2 ---------------------

The notebooks figure2, DisorderGraph, and HeatmapGeneraion are examples for generating particular figures

The utility.py has functions within it for finding the correct data and extracting it (see file for function descs).

--------------- 3 ---------------------
steps for submitting jobs at parameter values

1: $source /home/apatapof/softmatter/pydev/bin/activate
2: $cd /home/apatapof/softmatter/jobs
3: modify the job_walkthrough.py for parameters and note the arguments 
4:
	4a: $python3 job_walkthrough.py 50 0.1 0 1000       # this says to submit the same job 50 times at 0.1 disorder where the amplitude changes for the disorder and to start the folder indexing at 1000

	4b: $python3 job_walkthrough.py 			    # will run soley on the parameters you've assigned within the file. 
5: when simulations are finished you need to analyse the trajectory files that have been put into "/home/apatapof/softmatter/runs"
note: you must submit these jobs from within the jobs/ directory
	5a:you can do this by: $python3 rangecompile.py STARTINDEX FINALINDEX 			# this will analyze all of the runs between start and final index
	5b:you can do this by: $python3 spcomile.py INDEX 					# this will analyze the runs with the INDEX name
	5c:you can do this by: $python3 spcomile.py 						# this will analyze all the runs in the /runs/ folder
6: You'll notice that in the runs folder each INDEX has varied _sheer_a.b at the end. They now need to get merged to be a single INDEX folder in the compiledRuns folder

NOTE: you must now be in the "/home/apatapof/softmatter/runs" directory 
	$ python3 ../Polymer-Tilt/scatter_handler.py

this will take a while but it will merge all the analyzed data into a single directory.

--------------- 4 ----------------------

Data analysis happens in the compiler class.
add a function in the compiler.h and define it in compiler.cpp.
then run the linker_script.py, which will rebuild the binaries.

---------------END----------------------