
import math
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
    x = 0
    a = []
    hb = 0
    
    while(x <= disorder_ratio):
        x=0
        a=[]
        for i in range(nodes):
            am = (upper * hb - lower) + lower
            a.append(am)
        for element in a:
            x += element * element
        x = math.sqrt(x)/abs(amplitude)
        hb += 0.001
    return hb
