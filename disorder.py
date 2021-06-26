
import math 
def get_amount_nodes( disorder_ratio):
    def func(x):
        return 0.00703853 +  0.0594096 * pow(x, 0.50259308)
    done = False
    ratio = 0
    num = 0
    while (ratio < disorder_ratio):
        num += 1
        ratio = func(num)
    
    return num