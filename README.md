# Direct Collocation Methods for Constrained Robotic Systems
This repository contains the Matlab implementation of the methods described in the paper "Direct Collocation Methods for Trajectory Optimization in Constrained Robotic Systems", published in IEEE Transactions on Robotics. The paper can be downloaded from https://bit.ly/3OsqTZS

![dae_methods](https://github.com/rbordalba/collocation_constrained_systems/blob/main/doc/images/dae_methods.png)

# Prerequisites #

- Install HSL package for IPOPT (includes linear solver MA27 used in IPOPT). Download from:
http://www.hsl.rl.ac.uk/ipopt/

- Install 'CasADi v.3.5.5': tool for nonlinear optimization and algorithmic differentiation (includes the interior point optimizer IPOPT). Download from:
https://web.casadi.org/get/  (Add folder to your Matlab path)

- Tested with Matlab 2020b.

# Other dependencies #

- 'Spatial_v2 library': spatial vector arithmetic and dynamics algorithms.  
https://royfeatherstone.org/spatial/v2/ (A modified version is included in this repository.)

- 'IPOPT': interior point optimizer for large-scale optimization. 
https://coin-or.github.io/Ipopt/ (Included within CasADi)

# Run #

## Five-bars example ##
To see an animation of the optimized trajectory for the different methods (Fig. 10 in the paper), execute the following script 
```
run_fivebars_snapshots
```

To see the plots of the kinematic and dynamic errors for the different methods (Fig. 11 and 12 in the paper), execute the following script 
```
run_fivebars_plots
```

To see the plot of the evolution of the cost per iteration for the different methods (Fig. 13 in the paper), execute the following script 
```
run_fivebars_convergence
```

To see computed action trajectory for different cost functions with the projection method  (Fig. 14 in the paper),  execute the following script 
```
run_fivebars_snapshots
```

To see summary results for the different methods,  execute the following script 
```
run_fivebars_statistics
```



## Collaborative arms example ##
To see the plots of the optimized trajectory with the projection and local coordinates methods, execute the following script 
```
run_collaborativeArms_plot
```

To see summary results for the different methods,  execute the following script 
```
run_collaborativeArms_statistics
```
