import hoomd
import hoomd.md
import hoomd.md.update
import numpy
import math
hoomd.context.initialize("")

                                                    #creating and populating a snapshot
boundary = hoomd.data.boxdim(Lx = 100, Ly = 100,dimensions=2)
types = ['A','B']                                   # A: the defualt particle, B: the anchor particle
snapshot = hoomd.data.make_snapshot(N=10, box=boundary, particle_types=types, bond_types=['polymer'])

snapshot.particles.position[:] = [[0, -4.5, 0], [0, -3.5, 0],
                                  [0, -2.5, 0], [0, -1.5, 0],
                                  [0, -0.5, 0], [0, 0.5, 0],
                                  [0, 1.5, 0], [0, 2.5, 0],
                                  [0, 3.5, 0], [0, 4.5, 0]];



snapshot.particles.typeid[1:10]=0                   #A
snapshot.particles.typeid[0] = 1                    #B

snapshot.bonds.resize(9)                            #since we have 10 particles, there are only 9 linear bonds (N-1)
bonds = [[0,1], [1, 2], [2,3],[3,4], [4,5], [5,6],[6,7], [7,8], [8,9]] # defining the molecular bonds between particle ids
snapshot.bonds.group[:] = bonds                     #assigning bond groups to the snapshot



hoomd.init.read_snapshot(snapshot)                  #loading snapshot to hoomd
hoomd.md.update.enforce2d()

###########################################
all = hoomd.group.all()
anchor = hoomd.group.type(name='anchor', type='B')
chain = hoomd.group.type(name='chain', type='A')

nl = hoomd.md.nlist.cell();

dpd = hoomd.md.pair.dpd(r_cut=2.0, nlist=nl, kT=0.0, seed=3456) #I think the radius cut is the ran
dpd.pair_coeff.set('A', 'A', A=20.0, gamma = 1.0)
dpd.pair_coeff.set('A', 'B', A=20.0, gamma = 1.0)
dpd.pair_coeff.set('B', 'B', A=20.0, gamma = 1.0)

nl.reset_exclusions(exclusions = [])


harmonic = hoomd.md.bond.harmonic();
harmonic.bond_coeff.set('polymer', k=10.0, r0=0);

hoomd.md.constrain.oneD(group=anchor, constraint_vector=[1,0,0])
#hoomd.md.constrain.oneD(group=chain, constraint_vector=[1,0,0])#
#hoomd.md.constrain.oneD(group=chain, constraint_vector=[0,1,0])#

############################################### -TODO-: IMPLEMENTING UNIFORM FORCE AND BOUNDARY FORCES
hoomd.md.force.constant(group=all,fy=0.0,fx=0.0,fz=0.0)
const = hoomd.md.force.constant(fvec=(0.0,1.0,0.0))  # ?
periodic = hoomd.md.external.periodic()
periodic.force_coeff.set('A', A=-1000.0, i=0, w=1, p=1)
periodic.force_coeff.set('B', A=-1000.0, i=0, w=1, p=1)


###############################################
N = len(all)
activity = [ ( ((numpy.random.rand(1)[0]-0.5)),
               ((numpy.random.rand(1)[0]-0.5)),
               0)
             for i in range(N)];
hoomd.md.force.active(group=all,
                      seed=43,
                      f_lst=activity,
                      rotation_diff=0.005,
                      orientation_link=False);
###############################################


#hoomd.md.integrate.mode_standard(dt=0.01);
#integrator = hoomd.md.integrate.nve(group=all);
hoomd.md.integrate.mode_standard(dt=0.001);
hoomd.md.integrate.brownian(group=all, kT=0.0, seed=651);



#integrator.randomize_velocities(kT=0.8, seed=42)

################################################
hoomd.analyze.log(filename="log-output.log",
                  quantities=['potential_energy', 'temperature'],
                  period=100,
                  overwrite=True);
hoomd.dump.gsd("polymer.gsd", period=1e4, group=all, overwrite=True);

hoomd.run(10e5);
