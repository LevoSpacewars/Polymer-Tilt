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
        if '_' in filename:
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
    data = [[]]*7

    for element in input[key]:
        fn = str(path) + str(element) + "/data.txt"
        print(fn)
        file = open(fn,'r')
        index = 0
        print(file.readlines())
        for line in file.readlines():

            file_element = float(line.split(',')[-1].strip('\n'))
            data[index].append(file_element)
            index +=1
        file.close()
    
    print(data)
    exit()
    return data


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


if nameid is None:
    for key in dirdict:
        data = getData(path, dirdict, key)
        removeWriteDirectory(path + "/" + compdir ,key)
        dirname = makeWriteDirectory(path + compdir,key)
        dirpath = path + compdir + "/" + dirname
       
        os.system("cp " + path + str(dirdict[key][0]) + "/_simulation_parameters.txt " + dirpath)
        changeDF(dirpath + "/_simulation_parameters.txt", len(dirdict[key]))
        writeData(dirpath, data)

