import sys
import os
a = int(sys.argv[1])
b = int(sys.argv[2])

rlist = list(range(a,b + 1))

for element in rlist:
    os.system("python spcompile.py " + str(element))