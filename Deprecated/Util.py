class ID:
    Name = "NULL"
    objectName = "NULL"
    IN = -1             #identification number

def add(x,y):
    if len(x) != len(y):
        print("lengths are not the same!",len(x),len(y))
        return -1
    s = []
    for i in range(len(x)):
        s.append(x[i] + y[i])
    return s
