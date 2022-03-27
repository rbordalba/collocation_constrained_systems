# Direct Collocation Methods for Constrained Systems
This work contains implementations of different method for trajectory optimization for constrained systems. For more information please refer to the paper "Direct Collocation Methods for Trajectory Optimization in Constrained Robotic Systems".

![dae_methods](https://github.com/rbordalba/collocation_constrained_systems/blob/main/doc/images/dae_methods.png)

# Prerequisits #

- Install HSL package for IPOPT (includes linear solver MA27 used in IPOPT). Download from:
http://www.hsl.rl.ac.uk/ipopt/

- Install 'CasADi v.3.5.5': tool for nonlinear optimization and algorithmic differentiation (includes the interior point optimizer IPOPT). Download from:
https://web.casadi.org/get/  (Add folder to your Matlab path)

- Tested with Matlab 2020b.

# Other dependencies #

- 'Spatial_v2 library': spatial vector arithmetic and dynamics algorithms.  
https://royfeatherstone.org/spatial/v2/ (A modified version is included in the repository.)

- 'IPOPT': interior point optimizer for large-scale optimization. 
https://coin-or.github.io/Ipopt/ (Included within CasADi)

# Run #

## Five-bars example ##
To see an animation of the optimized trajectory for the different methods (Fig 10 in the paper), execute the following script 
```
run_fivebars_snapshots
```

To see the plots of the kinematic and dynamic error for the different methods (Fig 11 and 12 in the paper), execute the following script 
```
run_fivebars_plots
```

To see the plot of the evolution of the cost per iteration for the different methods (Fig 13 in the paper), execute the following script 
```
run_fivebars_convergence
```

To see solution action trajectory for different cost function with the projection method  (Fig 14 in the paper),  execute the following script 
```
run_fivebars_snapshots
```

To see summary results for the different methods,  execute the following script 
```
run_fivebars_statistics
```



## Collaborative arms example ##
To see the plots of the optimized trajectory with the projection and local method, execute the following script 
```
run_collaborativeArms_plot
```

To see summary results for the different methods,  execute the following script 
```
run_collaborativeArms_statistics
```
