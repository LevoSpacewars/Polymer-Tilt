from Analyzer import Analyze


analyzer = Analyze(names=["polymer_force_0"])

print("plotting probability")
analyzer.getPositionProbabilityData(rez=[100,100],Interval = 0.9,name = "p0")
print("finished \n")
print("plotting position")
#analyzer.plotPositionData(fps=500,Interval= 0)
print("finished")
# for i in range(1,11):
#     print("plotting force:" + str(i))
#     analyzer = Analyze(names=["polymer_force_" + str(i)])
#     analyzer.getPositionProbabilityData(rez=[50,50],Interval = 0.99 ,name = "force_" + str(i))
