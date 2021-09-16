from matplotlib import pyplot as plt
import glob
import numpy as np
import pandas as ps
from Simulations import PolymerSimulationParameters as PSP

BASEDIR = "compiledRuns/"


def getmatchingparameterruns(parameters):
    """
    parameters: dict(str : value)


    return: list of folders matching the wanted parameters

    see heatmapgeneration notebook for example. 
    if you need a list of what the parameters are then look in _simulation_parmaters.txt OR
    use the PolymerSimulationParameters.keys() function to see all the parameters

    """
    files = glob.glob(f"{BASEDIR}/*/_simulation_parameters.txt",recursive = True)
    
    folders = []
    
    for file in files:
        params = PSP()
        params.loadParameters(file)
        matches = True
        for key in parameters.keys():
            #print(params[key], parameters[key])
            matches *= (params[key] == parameters[key])
            
        
        if matches:
            folder = '/'.join(file.split('/')[0:-1])
            folders.append(folder)
    return folders

def getrunswithda(DA: str):
    """
    returns all the folders with the same DA value in the _simulation_parameters.txt
    """
    
    files = glob.glob(f"{BASEDIR}/*/_simulation_parameters.txt",recursive = True)
    folders = []
    
    for file in files:
        params = PSP()
        params.loadParameters(file)
        
        a = params.getDisorder()
        
        if a == DA:
            folder = '/'.join(file.split('/')[0:-1])
            folders.append(folder)
    return folders

def getheatmap(path, forcevalue):
    """
    Path: str
        path to the directory with heatmap
    forcevalue: str or float
        heatmaps contain the sheer values at the end. 
        remeber that those values are rounded to the first decimal
    
    returns: np.array or None
    """
    path = path.strip()
    if path[-1] == "/":
        path = path[:-1]
    try:
        print(f"{path}/heatmap_{forcevalue}.txt")
        return np.array(ps.read_csv(f"{path}/heatmap_{forcevalue}.txt" ,header=1))
    except:
        print("file not found")
        return None
        
def getdenmap(path, forcevalue):
    """
    Path: str
        path to the directory with denmap
    forcevalue: str or float
        denmap contain the sheer values at the end. 
        remeber that those values are rounded to the first decimal
    
    returns: np.array or None
    """
    path = path.strip()
    if path[-1] == "/":
        path = path[:-1]
    try:
        print(f"{path}/heatmap_{forcevalue}.txt")
        return np.array(ps.read_csv(f"{path}/ProfileDensity_{forcevalue}.txt" ,header=1))
    except:
        print("file not found")
        return None

def getoutput(path):
    """
    path: str
        path to the direcotry with the data.txt
    
    returns: (list output,list output_unc) or None
    """
    
    path = path.strip()
    if path[-1] == "/":
        path = path[:-1]
    try:
        data = np.array(ps.read_csv(path,header=0))
        return data[3][1:], data[4][1:]
    except:
        print(f"file not found or error with reading the file: {path}")
        return None
    
def getforceinterval(path):
    """
    returns a list of ordered sheers encountered over a particular run
    """
    path = path.strip()
    if path[-1] == "/":
        path = path[:-1]
    path += "/_simulation_parameters.txt"
    
    params = PSP()
    params.loadParameters(file)
    ab = params.getSheerForceRange()
    df = params.getDf()
    a = ab[0]
    b = ab[1]
    return np.linspace(a,b,df)

def getforceangles(path):
    """
    returns a list of ordered forceangles encountered over a particular run
    """
    path = path.strip()
    if path[-1] == "/":
        path = path[:-1]
    path += "/_simulation_parameters.txt"
    
    params = PSP()
    params.loadParameters(file)
    ab = params.getSheerForceRange()
    df = params.getDf()
    t = params.getPullForce()
    a = ab[0]
    b = ab[1]
    fr = []
    return np.linspace(a,b,df)/t

def getparameter(path,key):
    """
    return specific parameter value
    """
    path = path.strip()
    if path[-1] == "/":
        path = path[:-1]
    path += "/_simulation_parameters.txt"
    params = PSP()
    params.loadParameters(path)
    
    return params[key]

def getparameters(path,keys):
    """
    return specific parameter value
    """
    path = path.strip()
    if path[-1] == "/":
        path = path[:-1]
    path += "/_simulation_parameters.txt"
    params = PSP()
    params.loadParameters(path)
    
    return [params[key] for key in keys]
