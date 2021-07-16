
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

def get_amplitude_mod(disorder_ratio,nodes):
    a = []
    import math as r
    hb = 1
    lower = 0.1
    upper = 0.2
    x = 0
    while(x < disorder_ratio):
        a=[]
        for i in range(nodes):
            am = r.random() * (upper - lower * hb) + lower
            a.append(am)
        for element in a:
            x += element * element
        x = math.sqrt(x)
        hb += 1
