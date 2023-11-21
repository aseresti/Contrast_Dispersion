# Material and Methods
## Contrast Transport and Blood Flow Hemodynamics
The implementation of the TAFE method was discussed by detail in Eslami et al. [2015, 2022]. In short, two partial differential equations are the key components of implementation of this method: Navier-Stokes equations, simulating the blood flow assuming as a newtonian and incompressable fluid, and Advection-Diffusion equation, modeling the contrast agent transport in the blood. The governing equations of the Navier-Stokes equation is presented below:

$\dfrac{\partial U}{\partial t} + (U. \nabla)U + \dfrac{\nabla P}{\rho} = \nu \nabla^2U,\quad \nabla.U=0$

Where $U$ is the time-dependant velocity vector field, $P$ is the pressure field inside the geometry domain, $\rho$ is the blood flow density, and $\nu$ is the dynamic viscosity. A neuman type boundary condition was assigned to the inlet of the vessel and a non-slip boundary condition ($u=0$) was assigned to the walls. The advection-diffusion equation is defined as below where $C(\frac{mg}{ml})$ is the contrast agent concentration: 

$\dfrac{\partial C}{\partial t} + (U.\nabla )C = D\nabla ^2C$

Eslami et al [2022] in their work presented an equation for the quantitative flow rate estimate $Q_{TAFE}$ of the flow along the vessel. The equation was established based on assuming a unidirectional flow with minimal lateral flow distortion and the domainance advection in comparison to the diffusion. The method then suggest to use the cumulative volume at the axial locations down the vessel to find $Q_{TAFE}$. However, we used the same strategy to validate against the $v_{mean}$, the mean blood velocity along the centerline direction of the vessel.

$v_{mean} = \dfrac{\frac{\partial C}{\partial t}}{-\frac{\partial C}{\partial x}}$

The numerator of the above fraction, $\frac{\partial C}{\partial t}$ is the temporal derivative of the time-attenuation curve (TAC) with the $\frac{HU}{s}$ unit, and $\frac{\partial C}{\partial x}$ can be found with takeng the derivative of the contrast at the centerline of the vessel with respect to the centerline direction with the unit of $\frac{HU}{cm}$.

## Simulation Setup
\textit{CFD simulations} The direct numerical solution of the steady flow along a 3D non-axisymetric stenotic model were previously presented by M. Owais Khan et al. [2018]. Briefly, to characterize the turbulance transition in a stenosed artery we implemented the CFD simulations for Navier-Stokes equation on Oais flow solver [] for $Re=100-1000$. We used a second-order polynomial function for velocity and a first-order polynomial function for pressure. Fig 1 presents the axial and cross sectional view of the stenosed artery. The geometry was meshed in 6000 triangular elements. The simulation was run for 1 second with 10000 time steps. Fig 2 depicts the Q-criterion of the turbulance transition along the artery for $Re=100-1000$. 

\textit{Advection Diffusion Simulations} The advection diffusion equation was implemented using the FEniCs framework [] on the output of the CFD simulations. FEniCs is a library to implements PDEs based on the basis of the finite element methods on Python. To assign the inlet boundary conditions of the contrast agent we used the equation presented in Eslami et al. [2015]:

$ C_{inlet} =  C_{min} + 0.5 \times (C_{max}-C_{min}) \times [1-cos(\pi \times \dfrac{t-T_s}{T_d})]$

In this simulation $C_{min}$ and $C_{max}$ are the minimum and maximum of the concentration at the inlet. $T_s$ is the time that bolus arrives to the inlet, and $T_d$ is it takes for the bolus to arrive from the inlet to the outlet. In our simulations that we have a varying $Re$, we chose $T_d = 15-80s$. Fig 3 sheds more light on the contrast transition along the vessel.

\textit{Contrast Dispersion Simulations} The next step in our simulation to back-calculate the mean velocity according to waht presented in eq (3). In this simulation we fit a linear regression model to the time-attenuation curve at the inlet of the vessel, and the centerline-attenuation curve at the peak time of the advection diffusion output. The slope of the linear model fitted to each curves represents $\frac{\partial C}{\partial t}$ and $\frac{\partial C}{\partial x}$ respectively.