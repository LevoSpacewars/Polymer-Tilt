from Analyzer import Analyze


#analyzer = Analyze()

print("plotting probability")
#analyzer.getPositionProbabilityData(rez=[50,50],Interval = 0,name = "initial")
print("finished \n")
print("plotting position")
#analyzer.plotPositionData(fps=500,Interval= 0)
print("finished")
for i in range(1,11):
    print("plotting force:" + str(i))
    analyzer = Analyze(names=["polymer_force_" + str(i)])
    analyzer.getPositionProbabilityData(rez=[50,50],Interval = 0.9 ,name = "force_" + str(i))
