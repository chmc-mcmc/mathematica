mcmc: hmc, rmhmc (riemannian manifold hamiltonian monte carlo) , chmc (complementary hamiltonian monte carlo).

It is based on the following ideas:
1. the system is consisted of (a few if CPU, many if UCDA) particles
2. these particles can exchange energy by colliding virtually
3. the energy of the whole system is conserved
4. use "metropolis" algorithm
5. tune all parameter automatically
6. use both trajectories of HMC and RMHMC
