# Implementing Crank-Nicolson Scheme on Our Advection Diffusion Pipeline
## Why do we need to implement Crank-Nicoloson method?
One of the problems that we have faced in finding the velocity out of the contrast transition along the vessel was the error regardong the linear implementation of the temporal derivative. As a result of the linear implementation some error propagates to the further time-steps. Therefore, we need a more accurate implementation of the temporal derivative.
## Advection-Diffusion equation with Crank-Nicoloson implementation
Crank-Nicoloson method is applied on PDEs based on the equation below:

${\color{blue} \textstyle \dfrac{P^{\color{red}k} - P^{\color{red}k-1}}{\Delta t}={\color{BlueGreen}\theta} Q^{\color{red}k} + (1-{\color{BlueGreen}\theta})Q^{\color{red}k-1}}$

${\color{blue} \textstyle \theta = 0.5}$

the advection diffusion equation is presented below:

${\color{blue} \textstyle \dfrac{\partial C}{\partial t} + (U.\nabla)C = D\nabla^2C}$

The time-discrete version of this equation is implemented as below where ${\color{blue} \textstyle \Delta t}$ is the time discritization parameter:

${\color{blue} \textstyle \dfrac{C^{\color{red}k} - C^{\color{red}k-1}}{\Delta t} = -(u.\nabla)C^{\color{red}k}+D\nabla^2C^{\color{red}k}}$

the right-side of this equation represents the $Q$ in ${\color{blue} \textstyle \theta}$-scheme. Therefore, the Crank-Nicolson (${\color{blue} \scriptstyle \theta = 0.5}$) implementation of the advection diffusion is:

${\color{blue} \textstyle \dfrac{C^{\color{red}k} - C^{\color{red}k-1}}{\Delta t} = {\color{BlueGreen}0.5}(-(u.\nabla)C^{\color{red}k}+D\nabla^2C^{\color{red}k}) + {\color{BlueGreen}0.5}((-(u.\nabla)C^{\color{red}k-1}+D\nabla^2C^{\color{red}k-1}))}$

## FEniCs implementation of the Advection Diffusion: Crank-Nicolson Scheme

The linear implementation of the temporal derivatave on FEniCs is as below:

```
F = ((u - u_n) / k)*v*dx + dot(w, grad(u))*v*dx + D*dot(grad(u), grad(v))*dx - source*v*dx
```

Plus, the second order implementation of the temporal derivative using the crank-nikolson scheme would be as what follows:

```
Theta = Constant(0.5)
F = ((u - u_n) / k)*v*dx + Theta*dot(w, grad(u))*v*dx + Theta*D*dot(grad(u), grad(v))*dx + Theta*dot(w, grad(u_n))*v*dx + Theta*D*dot(grad(u_n), grad(v))*dx - source*v*dx
```

To study more about the ${\color{blue} \textstyle \theta}$-scheme and its implementation on FEniCs look at [link](https://home.simula.no/~hpl/homepage/fenics-tutorial/release-1.0-nonabla/webm/timedep.html) and [link](https://en.wikipedia.org/wiki/Crank–Nicolson_method)
