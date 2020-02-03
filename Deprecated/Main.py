import Populator
def fu(i,N):
    if i > 0:
        return [i-1,i]
    return None
test = Populator.Polymer(1)
test.createPolymerChain(10)
print(test.getPos())
test.defineBonds(func = fu)
print(test.getBonds())
