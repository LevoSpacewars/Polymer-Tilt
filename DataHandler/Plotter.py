import numpy as np;
from matplotlib import pyplot;
import os

DEFAULT_PARENT_PATH = "/home/alexander/Documents/Polymer-Tilt"

path = input("Please paste path to parent Dir:")
if path == "":
    path = DEFAULT_PARENT_PATH

traj_paths = []
objs = os.listdir(path)
for x in objs:
    p = path + "/" + x
    if os.path.isdir(p):
        if os.path.lexists(p + "/trajectory.gsd") and os.path.lexists(p + "/data.txt"):
            traj_paths.append(p)

print(traj_paths)




