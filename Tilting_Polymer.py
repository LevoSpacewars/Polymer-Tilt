import hoomd
import hoomd.md
import hoomd.md.update
import numpy
import math
import random
hoomd.context.initialize("")
random.seed()

def HardContact(r, rmin, ramx, l_0,K):
    V = K*(r-l_0)**2
    F = -2*K*(r-l_0)
    #print(F)
    return (V,F)


def Harmonicf(r, rmin, rmax,l_0,l_max,K):
    if r < l_0:

        V = K*(r-l_0)**2
        F = -2*K*(r-l_0)
        #print(F)
        return (V,F)
    # elif (1-(r-l_0)**2/l_max**2 ==0) or (1-(r-l_0)**2/l_max**2 < 0):
    #     V = K*(r-l_0)**2
    #     F = 2*K*(r-l_0)
    #     return (0,0)

    V = - K*l_max**2/2*math.log(1-(r-l_0)**2/l_max**2)
    F =  - (K*(r-l_0))/(1-(r-l_0)**2/l_max**2)
#    print(F)
    return (V,F)
def harmonicf(r, rmin, rmax, kappa, r0):
   V = 0.5 * kappa * (r-r0)**2;
   F = -kappa*(r-r0);
   return (V, F)
                                                    #creating and populating a snapshot
types = ['A','B','C']                                   # A: the defualt particle, B: the anchor particle
a = 1
length = 100
lines = 5
lbda = length
rez=0.3

K=10**3
l_0 = 0.3
l_max = 0.3
rmax = l_0 + l_max - 0.001
pull = 10
lpull =1
amplitude = 10

gamma = 0.5
kbT = 1
pos = []
bonds = []
id = []
mass = []


position_period = 10000;
log_period = position_period;
dt = 0.0001;
run_length = 10000000
backup = True

xbox = lines
ybox = 16*lbda
boundary = hoomd.data.boxdim(Lx = lines, Ly = 16*lbda,dimensions=2)
snapshot = hoomd.data.make_snapshot(N=lines*length, box=boundary, particle_types=types, bond_types=['polymer'])
# Position for each particle
for i in range(lines):
    x = random.uniform(0, lines) - lines/2
    y = 0
    for j in range(length):
        y = y + random.uniform(l_0, l_max)
        pos.append([x,y,0])
        print(pos[-1][0])
#defing bonds
for i in range(lines):
    for j in range(1,length):
        bonds.append([(j-1 + i*length),(j + i*length)])
print(bonds)

#defining ids
for i in range(lines):
    for j in range(length):
        if(j==0):
            id.append(0)
        elif (j==(length-1)):
            id.append(1)
        else:
            id.append(2)

#defining mass
for i in range(len(pos)):
    mass.append(1)

snapshot.particles.position[:] = pos
snapshot.bonds.resize(len(bonds))
snapshot.bonds.group[:] = bonds
snapshot.particles.typeid[:] = id
snapshot.particles.mass[:] = mass

hoomd.init.read_snapshot(snapshot)

#########################################################################

all = hoomd.group.all()
anchor = hoomd.group.type(name='anchor', type='A')
pulley = hoomd.group.type(name='pulley',type='B')
chain = hoomd.group.type(name='chain', type='C')

nl = hoomd.md.nlist.cell();

#########################################################################


table = hoomd.md.pair.table(width=100000, nlist=nl) #I think the radius cut is the ran
table.pair_coeff.set('A', 'A', func=HardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
table.pair_coeff.set('A', 'B', func=HardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
table.pair_coeff.set('B', 'B', func=HardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))

table.pair_coeff.set('B', 'C', func=HardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
table.pair_coeff.set('A', 'C', func=HardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
table.pair_coeff.set('C', 'C', func=HardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))

nl.reset_exclusions(exclusions = [])

harmonic = hoomd.md.bond.table(width = 10000);
harmonic.bond_coeff.set('polymer', func = Harmonicf,rmin=0, rmax=rmax,coeff = dict(l_0=l_0,l_max=l_max,K=K));
#harmonic.bond_coeff.set('polymer', func = harmonicf,rmin=0, rmax=100,coeff = dict(kappa=330,r0=0.84));
#fene = hoomd.md.bond.fene()

#fene.bond_coeff.set('polymer', k=10**3, r0=1, sigma=2, epsilon= 1.0)

#################################################################
#hoomd.md.constrain.oneD(group=chain, constraint_vector=[1,0,0])#
#hoomd.md.constrain.oneD(group=chain, constraint_vector=[0,1,0])#
#################################################################

############################################### -TODO-: IMPLEMENTING UNIFORM FORCE AND BOUNDARY FORCES

const = hoomd.md.force.constant(group = pulley, fvec=(0.0,pull,0.0))  # ?
periodic = hoomd.md.external.periodic()
periodic.force_coeff.set('A', A=amplitude, i=0, w=1, p=lines)
periodic.force_coeff.set('B', A=amplitude, i=0, w=1, p=lines)
periodic.force_coeff.set('C', A=amplitude, i=0, w=1, p=lines)

periodic.force_coeff.set('A', A=-10000000.0, i=1, w=1, p=5)

###############################################
N = len(all)

# force_strength = math.sqrt(2*gamma*kbT)
# activity = [ ( ((numpy.random.rand(1)[0]-0.5))*force_strength,
#                ((numpy.random.rand(1)[0]-0.5))*force_strength,
#                0)
#              for i in range(N)];
# hoomd.md.force.active(group=all,
#                       seed=43,
#                       f_lst=activity,
#                       rotation_diff=0.005,
#                       orientation_link=False);
###############################################

hoomd.md.integrate.mode_standard(dt=0.0001);
bs = hoomd.md.integrate.brownian(group=all, kT=1, seed=random.randint(0,999999));
bs.set_gamma('A', gamma=gamma)
bs.set_gamma('B', gamma=gamma)
bs.set_gamma('C', gamma=gamma)
#hoomd.variant.linear_interp([(0,0), (5000,1)]




#integrator.randomize_velocities(kT=0.8, seed=42)

################################################
hoomd.analyze.log(filename="energy.log",
                  quantities=['potential_energy', 'temperature'],
                  period=log_period,
                  overwrite=True);

hoomd.dump.gsd("polymer.gsd", period=position_period, group=all, overwrite=True);
hoomd.run(run_length)
for i in range(10):
    const = hoomd.md.force.constant(group = pulley, fvec=(lpull + i,pull,0.0))
    hoomd.dump.gsd("polymer_force_" + str(lpull + i) + ".gsd", period=position_period, group=all, overwrite=True);
    hoomd.run(run_length)
simulation_data_log = open("simulation_parameters_log.txt", 'a+')
simulation_data = open("simulation_parameters.txt", 'w')

from datetime import datetime
time = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

simulation_data_log.write("SIMULATION:" + time)
simulation_data_log.write("\nposition period:" + str(position_period) + "\nlog period:" + str(log_period) + "\ndt:" + str(dt) + "\nrun length:" + str(run_length) + "\nbox_dim:" + str(xbox) + "," + str(ybox) +  "\n\n\n")
simulation_data_log.close()

simulation_data.write("SIMULATION:" + time)
simulation_data.write("\nposition period:" + str(position_period) + "\nlog period:" + str(log_period) + "\ndt:" + str(dt) + "\nrun length:" + str(run_length) + "\nbox_dim:" + str(xbox) + "," + str(ybox) +  "\n\n\n")
simulation_data.close()

if (backup == True):
    import os
    backupname= "polymer_" + time + ".gsd"
    os.system("cp polymer.gsd GSD_backups/" + backupname)
