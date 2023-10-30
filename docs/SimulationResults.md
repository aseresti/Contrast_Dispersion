# Data and Results

In this ssection the results of our several tests on different phantoms are presented. The goal is to implement the contrast dispersion simulations on a pipe along which the contrast agent diffuses.

## a) Straight Pipe Simulations:
In this section we tried to implement the CFD and Advection Diffusion simulations on a straight pipe. The final goal is to calculate velocity using the Contrast Dispersion method presented by Eslami 2022. Our implementation contains of first generating the mesh of a straight pipe. The mesh is generated on SimVascular open source software. On SimVascular the Global Max Edge Size was assigned $0.15$ and the final mesh was created with $2,612,769$ number of elements. Figure 1 presents the mesh of the straight pipe used in our simulations. This mesh is along the x-axis and witha length of 81 cm bounding from $x \in [-9, 72] cm$. The diameter of the pipe is $d = 3 cm$. In the next step, we assign a velocity to the CFD simulations which are run using Oasis solver on Niagara. In the next step we run the advection diffusion simulation using fenics and running on a hard drive at CIMBL. Finally, we use our implementation of the contrast disperssion method to acquire the mean velocity of the fluid. The followings are the tests that we have run on the straight pipe.

![Alt text](<../Images/Screen Shot 2023-10-30 at 11.04.02 AM.png>)

### Test the run time
In this test, I changed the number of the cycles ran in the advection diffusion simulation from 5 cycles to 30 cycles. The duration of each cycle is 1 second and the number of time step in each cycle is 1000. The inflow boundary condition assigned to the pipe in each simulation changes with respect to the input number of the cycle using the following formula:
"(cmin+0.5*(cmax-cmin)*(1-cos(pi*((t-Ts)/(2*Td)))))",cmin=0.0,cmax=1.0,t=0.0,Ts=0.0,Td=self.Args.PeriodContrast,degree=2

$
c_{min} + 0.5 \times (c_{max}-c_{min}) \times [1-cos(\pi \times \dfrac{t-T_s}{2 \times T_d})]
$

In which the $c_{min} = 0.0$, $c_{max} = 0.0$, $T_{s} = 0.0$
### Test the centerline vs sectional average values
### Test the Reynold's number

## b) Stenotic Pipe Simulation: Changing the Reynolds Number in the inlet of the geometry

## c) CT-MPI images results