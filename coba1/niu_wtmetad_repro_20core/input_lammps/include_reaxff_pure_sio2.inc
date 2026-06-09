# ReaxFF setup for pure Si/O data files with atom types: 1=Si, 2=O.
# If your LAMMPS uses the older USER-REAXC syntax, replace pair_style reaxff with pair_style reax/c.
pair_style      reaxff NULL safezone 1.6 mincap 100
pair_coeff      * * ${POTENTIAL} Si O
fix             fqeq all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff maxiter 200
