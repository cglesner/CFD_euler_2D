# CFD_euler_2D
Created by Colin Glesner \n
This code is a from scratch implementation of a computational fluid dynamics simulator.

Code Usage:

  - At present this code is hard-coded to use input grids and to generate output data
    formatted for use with the Tecplot visualization package. The next project for this code
    is to convert to VisIt or some other open-source visualization package. 
  - The executable can be built using make commands. Make Run_Simulation builds the executable.
  - At present the code is hard-coded to read input from a file named sim_config.txt, an 
    immediate priority is to change the code so that a command line argument can be passed to 
    the executable to specify the appropriate input file.
  - At present the code is only able to simulate one of three scenarios:
    * Method of manufactured solutions code verification scheme
    * Supersonic Inlet with angled ramp
    * Airfoil with varying angles of attack
    The lack of flexibility in the code lies in the presently very primitive input file,
    which only specifies which of the three hard-coded scenarios is desired. A more flexible
    input file with the ability to specify the type of boundary conditions at time of use
    will be a more challenging project, but one that would vastly increase the uses of the code.

This is a pretty basic code so far still, but I'm pretty proud of it, and hope to continue
working on improving it and expanding it, for my own education and just for fun. Some long
term goals for this code:

   - Convert to an object based implementation for maximum flexibility at time of use.
   - Add ability to perform axisymetric simulations
   - Add the ability to perform full 3d simulations
   - Implement openMPI or some other parallelization scheme to improve performance
   - Generally try to improve simulation speed
   - Experiment with line implicit and other more advanced numerical techniques to improve code
     performance.

Thanks for checking out my project! 
Colin Glesner
cglesner@vt.edu

