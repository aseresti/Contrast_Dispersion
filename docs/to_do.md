# to do list of this project
The goal of this repository is to contain all of the scripts and their test uses in contrast dispersion.
to dos:
 - modify the scripts to a not case specific class
 - add test folder to test the process further
 - Laminar flow study
 - turbulant flow study + develope scripts
 - patients data study
 - develop back calculation models:
    - Parastou's
    - higher-order
    - fft
 - update readme with proper abstract information
 - add proper docs:
    - Background
    - Our Avection Diffusion pipeline validation against the theoretical solution
    - method (Latex?)
    - results (either md or excel)
    - discussion

Friday Sept 29 TODOs:
   - Straight Pipe Segmentation and Meshing: **Done**
   - Generate Oasis Mesh **Done**
   - Submit CFD for a Re of 500 (Check Re of Parastou's) **Submitted**
   - Advection Diffusion for the laminar flow 
   - fft func **Done**

Sat Sep 30 TODOs:
   - Advection Diffusion for the laminar flow **Done**
   - fft func **Done**
   - Parastou's method results for data, Results excel sheet **Done**


Wed Oct 04 TODOs:
   - find the issue with velocity **Done**
   - filter model **In progress** >> **Done**
   - full resolution contrast dynamics **NFN**
   - fix CFD for straight pipe **Done** 

Fri Oct06 TODOs:
   - rendering of the square mesh **Done**
   - apply 2d fft on the unified square mesh
   - 2d lowpass filter **Done**
   - study the effect of the turbulance on frequency domain data **Done** applying a filter doesn't have any specific effect on the results
   - Advection Diffusion of the straight pipe **Done**

Mon Oct09 TODOs:
   - Contrast Dispersion on a Taper vessel + using accumulative volume **Not For Now**
   - implement accumulative volume method **NFN**

Wed Oct11 TODOs:
   - Do we need to increase the time cycle of the Advection Diffusion? By how much? **Done**

Tue Oct24 TODOs:
   - The centerline results are slightly different from cross-sectional average.
   Need to adapt the code. **Done**

Mon Oct30 TODOs:
   - Run all of the advection diffusion simulations for 20 cycles with 1000 time steps in each cycle
   - assign a sphere with a dimater of the CT-MPI spatial resolution to the centerline of the lumen.
   - adapt the Advection Diffusion code to move the results to the hard drive at the end of the simulation.
   - test 20 cycle results for the re of 100 and re of 1500.
   - explain in an Md why the calculation error increases as we increase the number of time cycles.
   - Update Simulation Results Md **Done**

Wed Nov1 TODOs:
   - Update Simulation Results md with pictures
   - check the Re 500 for about 13 time cycles.
   - Update Stenotic pipe simulation results in simulation results md.