
import math
from random import random
def get_amount_nodes(disorder_ratio):
    def func(x):
        return -0.00203498 +  0.0309991 * pow(x, 0.49721885)
    done = False
    ratio = 0
    num = 0
    while (ratio < disorder_ratio):
        num += 1
        ratio = func(num)

    return num

def get_amplitude_mod(disorder_ratio,nodes,amplitude):
    
    upper = 0.005
    lower = 0.005
    am = 0
    hb = 0
    index = 0
    print(disorder_ratio)
    
    while am < disorder_ratio or index < 10:
        am = 0
        
        for i in range(nodes):
            am = am + pow( random() * (hb * upper) + lower, 2.0)
        am = math.sqrt(am)/abs(amplitude)
        #print(am,hb)
        hb = hb + 0.01

        if am > disorder_ratio:
            index += 1
    print(am,hb)

    am = 0
        
    for i in range(nodes):
        am = am + pow( random() * (hb * upper) + lower, 2.0)
    am = math.sqrt(am)/abs(amplitude)
    print(am,hb)
    return hb
