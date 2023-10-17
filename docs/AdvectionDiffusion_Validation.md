# Validating Our Advection Diffusion Simulations
In this section I explain the validation of the advedtion diffusion simulation implemented using fenics against the 1D analytical solution of the advection diffusion equation presented in [Arzani, 2021]. To understand the 3D Advection Diffusion simulation implemented read the simulation pipelin [document](SimulationPipeline.md). Considering a steady-state advection diffusion simulations, the 1D analytical solution for the advection diffusion is presented below:

$v \frac{\nabla C}{\nabla x} = D \frac{\nabla^2 C}{\nabla^2 x} $

Where $v$ is the fluid  velocity in the x-axis, and $D$ is the diffusion coefficient. For $ 0<x<1 $ the analytical solution for the advection diffusion is:

$C(x) = \dfrac{c_2 - c_1}{e^\frac{v}{D}-1}e^{\frac{v}{D}x}-\dfrac{c_1e^\frac{v}{D} - c_2}{e^\frac{v}{D}-1} $

The implementation of the analytical solution is on the [AnalyticalSolution](../scripts/Validation/OasisAnalyticalSolution.py) script. In this script using frnics library, an IntervalMesh is defined with an specific number of elements along a min and max bounds. Then a Lagrange element type is defined for a function space to interpolate the equation 2 on the defined function space.
Then, using the [AdvectionDifusionValidation](../scripts/Validation/OasisAdvectionDiffusionValidation.py) script, a 3D mesh which is a result of advection difusion simulation is taken. Then, for $ 0<x<1 $ 100 slices was selected and the average value of the contrast concentration were stored. Then the average value of the first slice was assigned to $c_1$ and that of the last slide was assigned to $c_2$. Using the analytical solution of the advection diffusion then were calculated and the values from 1D solution were compared to those of the averaged method. The mean error between the two method were $1.34 \%$.
