import os
import sys




def getDirectoryDict(directories):
    runidp = {}
    while(len(raw_dirs) != 0):
        name = raw_dirs.pop()
        key = name.split('_')[0]
        if key not in runidp:
            runidp[key] = []
        runidp[key].append(name)
    print(runidp)
    print("/")
    return runidp


def renameFiles(key, input):
    for element in input[key]:

        p = str(path) + str(element)
        sheer = round(float(element.split('_')[-1]),1)
        for file in os.listdir(p):
            if "heatmap" in file:
                name = "heatmap_" + str(sheer) + ".txt"
                os.system("mv " + p + "/" + file + " " + p + "/" + name)
            if "potential" in file:
                name = "potential_" + str(sheer) + ".png"
                os.system("mv " + p + "/" + file + " " + p + "/" + name)
            if "profiledensity" in file.lower():
                name = "ProfileDensity_" + str(sheer) + ".txt"
                os.system("mv " + p + "/" + file + " " + p + "/" + name)
            if "profiledata" in file.lower():
                name = "profiledata_" + str(sheer) + ".txt"
                os.system("mv " + p + "/" + file + " " + p + "/" + name)
            


def handleHeatmaps(key, input, destination):

    for element in input[key]:

        fn = str(path) + str(element) + "/heatmap.txt"
        p = str(path) + str(element)
        for file in os.listdir(p):
            if "heatmap" in file:
                sheer = round(float(element.split('_')[-1]),1)
                os.system("cp " + p + "/" + file + " " + destination)
            if "potential" in file:
                os.system("cp " + p + "/" + file + " " + destination)
            if "disorder_info" in file:
                os.system("cp " + p + "/" + file + " " + destination)
            if "ProfileDensity" in file:
                os.system("cp " + p + "/" + file + " " + destination)
            if "profileData" in file:
                os.system("cp " + p + "/" + file + " " + destination)



def sortDictArrays(input):
    for key in input.keys():
        values = list(float(item.split("_")[-1]) for item in input[key])
        values.sort()
        values = list(key + "_sheer_" + str(item) for item in values)
        input[key] = values
    print(input)
    return input


def getDirs(path):
    raw_dirs = []


    for filename in os.listdir(path):
        print(filename)
        if '_' in filename and os.path.isfile(filename + '/data.txt'):
                print("Name-ID exists")
                print("adding " + filename + " to list")
                raw_dirs.append(filename)
    return raw_dirs


def removeWriteDirectory(path, writeDirectoryName):
    for objects in os.listdir(path):
        if objects == writeDirectoryName:
            print("Name-ID exists")
            print("deleting contents of " + writeDirectoryName)
            os.system("rm -r " + writeDirectoryName)


def makeWriteDirectory(path,nameid):
    os.system("mkdir " + path +"/" + nameid)
    return nameid


def getData(path, input, key):
    output = [[],[],[],[],[],[],[]] #NOTE: this is not equal to [[]]*7

    for element in input[key]:
        fn = str(path) + str(element) + "/data.txt"
        print(fn)
        file = open(fn,'r')


        lines = file.readlines()
        for line_index in range(len(lines)):
            file_element = float(lines[line_index].split(',')[-1])
            output[line_index].append(file_element)
        file.close()

    return output


def writeData(dirpath, data):
    keys = ["dx","udx","length","ulength","output","uouput","dx2"]
    file = open(dirpath + "/data.txt",'w')
    for i in range(len(data)):
        line = str(keys[i])
        for element in data[i]:
            line += "," + str(element)
        line +="\n"
        file.write(line)
    file.close()


def changeDF(filepath,df):
    file = open(filepath, 'r')
    lines = file.readlines()
    for i in range(len(lines)):
        if 'df' in lines[i]:
            lines[i] = lines[i].split('=')[0] + "=" + str(df) + "\n"
            break
    file.close()

    file = open(filepath, 'w')
    file.writelines(lines)
    file.close()

def changeRange(filepath, init_sheer, final_sheer):
    file = open(filepath, 'r')
    lines = file.readlines()
    for i in range(len(lines)):
        if 'sheerForceRange' in lines[i]:
            lines[i] = lines[i].split('=')[0] + "=" + str(init_sheer) + "," + str(final_sheer) + "\n"
            break
    file.close()

    file = open(filepath, 'w')
    file.writelines(lines)
    file.close()

def getSheerRange(sorted_dict,key):
    min = -1
    max = -1
    items = sorted_dict[key]
    min = float(items[0].split('_')[-1])
    max = float(items[-1].split('_')[-1])
    return min, max

compdir = "compiledRuns"
nameid = None
path = "./"
if (len(sys.argv) > 1):
    nameid = sys.argv[1]
    if len(sys.argv) > 2:
        path = sys.argv[2]

# first check for if the directly nameid exists, so that we can delete it for recompilation
if nameid is None:
    removeWriteDirectory(path,compdir)
    makeWriteDirectory(path,compdir)
else:
    removeWriteDirectory(path,nameid)
    makeWriteDirectory(path,nameid)

#know make a list of all directories with the nameid within the title
raw_dirs = getDirs(path)

#now construct a dict based on name

dirdict = getDirectoryDict(nameid)
dirdict = sortDictArrays(dirdict)
print(dirdict)

if nameid is None:
    for key in dirdict:
        data = getData(path, dirdict, key)
        removeWriteDirectory(path + "/" + compdir ,key)
        dirname = makeWriteDirectory(path + compdir,key)
        dirpath = path + compdir + "/" + dirname

        os.system("cp " + path + str(dirdict[key][0]) + "/_simulation_parameters.txt " + dirpath)
        smin, smax = getSheerRange(dirdict, key)
        changeDF(dirpath + "/_simulation_parameters.txt", len(dirdict[key]))
        changeRange(dirpath + "/_simulation_parameters.txt", smin, smax)
        renameFiles(key,dirdict)
        handleHeatmaps(key, dirdict, dirpath)
        writeData(dirpath, data)
