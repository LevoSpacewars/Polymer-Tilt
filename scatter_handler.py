import os
import sys
nameid = sys.argv[1]
path = "./"
if (len(sys.argv) > 2):
    path = sys.argv[3]

# first check for if the directly nameid exists, so that we can delete it for recompilation
for objects in os.listdir(path):
    if objects == nameid:
        print("Name-ID exists")
        print("deleting contents of " + nameid)
        os.system("rm -r " + nameid)

#know make a list of all directories with the nameid within the title
runidp = []
for objects in os.listdir(path):
    if nameid in objects:
        print("Name-ID exists")
        print("adding " + objects + " to list")
        runidp.append(path + objects)
print(runidp)


# now let's order the sheer values, so we can then just append the values from the data.txt together
values = list(float(item.split("_")[-1]) for item in runidp) 
values.sort() 
values = list("_" + str(item) for item in values)
print(values)
# compare the order
dirs = []

#sort runidp
for i in values:
    dirs.append(list(filter(lambda x : i in x, iter(runidp)))[0])



#mkdir and extact data in order of all the sheers
os.system("mkdir " + path + nameid)
data = [[],[],[],[],[],[]]
for directory in (dirs):
    x = 0
    f = open(directory + "/data.txt", "r")
    for line in f.readlines():
        data[x].append(float(line.split(",")[1]))
        x = x + 1

#now that all the data.txts have been extracted, write data.txt into created directory
keys = ["dx","udx","length","ulength","output","uouput"]
f = open(path + nameid + "/data.txt",'w')
lines = []
for i in range(len(keys)):
    lines.append(keys[i] + "," + str(data[i]).strip("[").strip("]").replace(" ", "") + "\n")

f.writelines(lines)
f.close()

