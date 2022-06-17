# Purpose
Dynamics based MCMC

# Codebase for
* "A Multi-Trajectory Monte Carlo Sampler"
* "Robust Inference Based On the Complementary Hamiltonian Monte Carlo"
* orthogonal multi-trajectory sampler
* method4/Sampler4

# Example
* example.nb
* sampler.wl

# Theory

Total energy:

$H(q,p)=K(p,q)+U(q)$.

Hamiltonian Dynamics:

$\dot{q}=H_{p}=K_{p}$,

$\dot{p}=-H_{q}=-U_{q}-K_{q}$.

Euler integration:

$q_{1}= q_{0}+ \delta\times K_{p}(p_{0},q_{0})$,

$p_{1}= p_{0} - \delta\times \left( U_{q}(q_1)+K_{q}(p_0,q_1) \right)$.

The Kinetic energy $K_r$:

$K_{r}(p,q)=K(p,q,r)=\frac{1}{2}q^{T}\cdot U_{qq}^{-r}\cdot q$.

For simplicity, $r$ is a parameter (unlike $K_p$ and $K_q$).

For $r$ can be a decimal, Eigen decomposition is carried out:

$U_{qq}=V\cdot\Lambda\cdot{}V^{T}$.

For non-positive (and positive) definite Hessian:

$U_{qq}^{-r}:=V\cdot(\left|\Lambda\right|^{-r} \odot \text{sign}(\Lambda)\cdot{}V^{T}$.

The kinetic energy function that permits correct sampling of multivariate normal distributions utilizes the $K_{0.5}$ .

If the system comprises several particles with collision effects, it is intuitively apparent that there is no difference in jumping probabilities between the simulation trajectory's two terminals. Therefore, for a multi-particle system with collision, if total energy conservation is enforced during collision, the mutual jump probabilities of any particle's position at either end of the trajectory are equal, and the following equation follows

$P(q_{0 } \to q_{1}) = P(q_{1 } \to q_{0})$.

Thus, the Metropolis algorithm's acceptance probability can be utilized to determine if a particle jumps to the new position:

$\alpha=e^{U(q_{1})-U(q_{0})}$.

# Usage

1. load the sampler.

```
SetDirectory[NotebookDirectory[]];
<< "Sampler`";
```
2. Ignore the $K_q$ (optionally) to speed up computation. 

```
Kq[p_, q_, r_] = 0;
```

3. define the potential energy, i.e. the negative logarithm of probability density function, e.g. a bivariate normal distribution:

$\Sigma=\left(
\begin{array}{cc}
 1 & 0.999999999999999 \\
 0.999999999999999 & 1 \\
\end{array}\right)$

```
rho = 1 - 1/10^15;
SIGMA = {{1, rho}, {rho, 1}};
U[x_, y_] = 1/2 Simplify[{x, y}.LinearSolve[SIGMA, {x, y}]];
```
4. Define the derivative functions. 

```
Uq[x_, y_] = D1[U[x, y], {x, y}];
Uqq[x_, y_] = D2[U[x, y], {x, y}];
Uqqq[x_, y_] = D3[U[x, y], {x, y}];

```

5. Run the sampler with $K_{0.5}$.
```
Dim = 2;
BURNIN=5000;
ITERATIONS=10000;
QS = hmc[U, Uq, Uqq, Uqqq, Dim, BURNIN, ITERATIONS, {.5}, {}];
```

6. Check the result.

```
QS1 = QS.MatrixPower[SIGMA, -.5];
StandardDeviation[QS1]
ListPlot[{QS, QS1}, PlotStyle -> Opacity[1], AspectRatio -> 1, PlotLegends -> {Samples, Transformed}]
```

> {0.99964, 1.0068}

![scatter plots](bn2.png "scatter plots")

# Additional

+ For difficult models, several kinetic energies can be adopted, e.g.:

```
QS = hmc[U, Uq, Uqq, Uqqq, Dim, BURNIN, ITERATIONS, {.45,.5,.55}, {}];
```

+ Run-time parameter can be set:

```
CHAINS=5;
STEPS=6;
```

+ Flat prior can be implemented by setting a reject region, e.g.:

```
outbnd[q_] := q[[-1]] <= 0;
```

+ Provide initial value.

```
qinit = RandomVariate[UniformDistribution[], {CHAINS, Dim}];
QS = hmc[U, Uq, Uqq, Uqqq, Dim, BURNIN, ITERATIONS, {.5}, qinit];
```
