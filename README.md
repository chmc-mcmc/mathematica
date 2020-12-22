mcmc: hmc, rhmc (riemannian manifold hamiltonian monte carlo) , chmc (complementary hamiltonian monte carlo).
主要基于如下想法:
1. 一个系统包含几个原子
2. 这些原子相互碰撞交换能量
3. 碰撞后的方向为等可能
4. 整个系统的总能量恒定
5. 接受概率与总能量成反比
6. 交替使用两种仿真轨迹，以兼顾所有主成分方向
7. 使用metropolis算法
