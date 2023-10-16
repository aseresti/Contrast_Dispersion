# Simulation Pipeline

### Theory
Mainly the goal of this project was to calculate average blood velocity based on CT-MPI dataset using the Advection-Diffusion equation. The advection diffusion equation is defined as below where $C$ is the concentration of the contrast agent injected in the fluid, and $U$ is the velocity of the fluid:

$\partial C/\partial t + (U.\nabla )C = D\nabla ^2C$

Below is the one dimentional advection diffusion equation along the z-axis where $w$ is the velocity component of the fluid along the z-axis:

$dC/dt + w*dC/dz = D*d^2C/dz^2$

Based on the work presented by Eslami 2022, by neglecting the diffusion force, and assuming that the flow in the tube is unidirectional and it has a minimal variation across each cross section, the velocity of the fluid isestimated as below:

| $w = -\dfrac{\frac{dc}{dt}}{\frac{dc}{dz}} $ | eq. 3 |
| --- | --- |


Both of the $\frac{dc}{dt}$ and $\frac{dc}{dz}$ can be found using an instantanous derivative of the upslope of the time attenuation curve and the length attenuation curve.

### Simulations
To implement the contrast dispersion method on the contrast concentration data firstly acquires a geometery to be assigned. We used two geometries to test our pipeline on: i) A straight pipe geometry ii) A stenotic pipe geometery. Both pipes have a 81 mm length along the x-axis. The stenotic pipe has a $ 90\% $ stenosis at 9 mm after the inlet. The CFD simulations was implemented on the geometry using the Oasis solver. Next, using FEniCs library we implemented Advection Diffusion equation on the fluid velocity. Last step is to apply the Contrast Dispersion method on the contrast concentration data extracted from the Advection Diffusion Simulations. The details of the simulation can be found below:  

#### 1. Oasis CFD simulations

#### 2. Fenics Advection-Diffusion Simulations

#### 3. Contrast Dispersion Simulation
The simulation pipeline for the contrast dispersion method is implemented in [ContrastDispersion.py](../scripts/ContrastDispersion.py) script. This script takes the results folder containing the fenics simulation of the Advection-Diffusion simulations. It takes 10 concentration files linearly spaced in time to extract the upslope of the time-attenuation curve in the inlet cross-section of the pipe. Also, it takes the concentration at peak. To extract the length-attenuation curve is then extracted by taking 100 cross sections along the length of the pipe. Then, the average of the contrast concentration is taken to form the length-attenuation curve. Using the linear model package in the ScikitLearn library, I applied a linear model on both time-attenuation and length-attenuation curves. The slope of both models represent the istantanous derivatives of Contrast over temporal or spatial domain. Next step is to use eq. 3 to estimate the value of the average fluid velocity in the pipe.

#### 4. Data and Results

##### a) Straight Pipe Simulations: Test runtime

##### b) Stenotic Pipe Simulation: Changing the Reynolds Number in the inlet of the geometry

##### c) CT-MPI images results
