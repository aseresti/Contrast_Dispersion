# Data and Results

In this ssection the results of our several tests on different phantoms are presented. The goal is to implement the contrast dispersion simulations on a pipe along which the contrast agent diffuses.

## a) Straight Pipe Simulations:
In this section we tried to implement the CFD and Advection Diffusion simulations on a straight pipe. The final goal is to calculate velocity using the Contrast Dispersion method presented by Eslami 2022. Our implementation contains of first generating the mesh of a straight pipe. The mesh is generated on SimVascular open source software. On SimVascular the Global Max Edge Size was assigned $0.15$ and the final mesh was created with $2,612,769$ number of elements. Figure 1 presents the mesh of the straight pipe used in our simulations. This mesh is along the x-axis and witha length of 81 cm bounding from $x \in [-9, 72] cm$. The diameter of the pipe is $d = 3 cm$. In the next step, we assign a velocity to the CFD simulations which are run using Oasis solver on Niagara. In the next step we run the advection diffusion simulation using fenics and running on a hard drive at CIMBL. Finally, we use our implementation of the contrast disperssion method to acquire the mean velocity of the fluid. The followings are the tests that we have run on the straight pipe.

![Straight Pipe](<../assets/Screen Shot 2023-10-30 at 11.04.02 AM.png>)

### i) Test the run time
In this test, I changed the number of the cycles ran in the advection diffusion simulation from 5 cycles to 30 cycles. The duration of each cycle is 1 second and the number of time step in each cycle is 1000. The inflow boundary condition assigned to the pipe in each simulation changes with respect to the input number of the cycle using the following formula:

$
c_{min} + 0.5 \times (c_{max}-c_{min}) \times [1-cos(\pi \times \dfrac{t-T_s}{2 \times T_d})]
$

In which the $c_{min} = 0.0$, $c_{max} = 0.0$, $T_{s} = 0.0$, and $T_{d} = number\ of\ run\ cycles$ which in our simulations varies from 5 to 20. The output function of this equation starts rising at 0 from 0 to 0.5 by the end of the last time cycle. The figure below shows how the assigned inflow function:


Through doing several simulations with varying run-time from 5, 10, 20 to 30, in a pipe with a reynold's number of 500, I came to the results that the run-time of the Advection-Diffusion simulation should be selected with respect to the reynold's number assigned in the CFD simulations. The following shows the results for pipe with a reynold's number of 500 with 5, 10, 20, and 30 run-times in advection-diffusion simulation. The results contains the temporal and spatial attenuation curves using the values extracted along the centerline of the pipe from the centerline and the cross-sectional averages. Also, the velocity values were acquired using the centerline, and cross-sectional average values.


The results shows that the most accurate velocity is the one extracted from the results of the advection-diffusion simulations with 20 cycles. With 5 and 10 run-time cycles, it is obvious that the contrast agent has not yet reached the end of the pipe, while 30 cycles the contrast agent has reached the end of the pipe long before. As for the first two, the model fit on the spatial attenuation curve is not accurate. While, for the latter one, the calculation errors of the linear implementation of the time derivative flows into the further time steps as we let the simulation continue further. The first equation below shows the linear approximation of the function derivative using the taylor series. The second equation shows as we loop over the derivative of the function, the error propagates to the next iteration:

$ \dfrac{dC}{dt} = \dfrac{\Delta C}{\Delta t} + E $

$ \dfrac{dC^i}{dt} = \dfrac{\Delta C^i}{\Delta t} + \sum E_i = \dfrac{C^i-C^{i-1}}{\Delta t} + \sum E_i  $

As a result the velocity acquired from advection-diffusion simulation ran for 20 cycle, is the most accurate one. We need to consider that for the reynold's number of 500, the assigned velocity was about $6.66 \frac{cm}{s}$. Therefore, for the contrast bulos to reach the end of the pipe it takes $ \dfrac{81 cm}{6.66 \frac{cm}{s}} = 12.16 s $. Therefore, I tested running the advection-diffusion simulation for 13 time cycles.

### ii) Test the centerline vs sectional average values

The results of the velocity from centerline vs sectional average values shows that the velocity is more accurate when we use the centerline extracted values. As a result, we implement the pipelines to use centerline values instead of the sectional average values as was presented in Eslami 2022 work.

### iii) Test the Reynold's number

In a straight pipe, we ran the CFD simulation for a reynold's number of 500 and a reynold's number of 1500. The results of the $Re \# = 500 $ was presented above. The results for the $Re\# = 1500$ can be found below. We ran the advection-diffusion simulations for 10 time cycles. 

## b) Stenotic Pipe Simulation: Changing the Reynolds Number in the inlet of the geometry