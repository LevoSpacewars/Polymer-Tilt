import numpy as np
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import pandas
import matplotlib.cm as cmx
from scipy import stats
import sys
import os
import glob
import subprocess
import matplotlib.backends.backend_pdf


def makeprofiledist(rawfilepath):
    heatmap = pandas.read_csv(rawfilepath, header=1)
    name = rawfilepath.split(':')[1][0:-4]
    print(name)
    map = np.array(heatmap)
    x = []
    y = []
    for i in range(len(map)):
        x.append(map[i][0])
        y.append(map[i][1])
    fig = plt.figure()
    bins,a,b,d= plt.hist2d(x,y,bins=(100,100))
    plt.title("Theta: " +str(name))
    pdf.savefig(fig)
    plt.close()







print("RunDirectory:" + sys.argv[1])

dfgdsgfg


runDir = sys.argv[1]

files = os.listdir(runDir)
runDirs = []
for obj in files:
    if(os.path.isdir(obj)):
        runDirs.append(obj)
#print(runDirs)
for i in range(len(runDirs)):
    pdf = matplotlib.backends.backend_pdf.PdfPages(runDirs[i] + "/ProfileDensity.pdf")
    files = glob.glob(runDirs[i] + '/DensityData_raw*.txt')
    files_ordered = []
    order = []
    org_order=[]
    for j in range(len(files)):
        name = float(files[j].split(':')[1][0:-4])
        files_ordered.append(name)
    org_order=files_ordered
    files_ordered = sorted(files_ordered, key = lambda x:float(x))
    for j in range(len(files)):
        order.append(org_order.index(files_ordered[j]))
    
    for j in range(len(order)):
        print(files_ordered[j])
	print(files[j])
	print(files[order[j]])
	print("\n\n")

    for j in range(len(files)):
        makeprofiledist(files[order[j]])
    pdf.close()
   

