# Optimal Control of the Mars Ascent Vehicle via Convex Optimization

The techniques used within this code were based on similar methods implemented in the PhD thesis by Xinfu Liu titled "Autonomous Trajectory Planning by Convex Optimization", Iowa State University, 2013.

The following algorithm implements a convex optimization approach to allow the Mars Ascent Vehicle (MAV) to autonomously plan and execute a fuel optimal ascent trajectory, from the surface of Mars to a stable circular orbit. The optimal guidance algorithm makes use of convex relaxations and sequential convex programming (also known as successive convexification) to solve a series of convex subproblems. This enables rapid progress towards an optimal solution that minimizes fuel use along the MAV's ascent trajectory.

Further techniques were also implemented from the following papers:
1) Successive Convexification for Fuel-Optimal Powered Landing with Aerodynamic Drag and Non-Convex Constraints. Szmuk M., Acikmese B., Berning A., 2016
2) Solving Nonconvex Optimal Control Problems by Convex Optimization. Liu X., Lu P., 2013



**Note that the following software is required to run the code contained within this repo:**
1) YALMIP optimization environment (https://yalmip.github.io/)
2) MOSEK optimization package (https://www.mosek.com/)
3) A norms function to allow for computation of multiple vector norms

## File descriptions:
* `mav_launch.m`: main algorithm used to implement and run autonomous guidance and control of the Mars Ascent Vehicle
* `controller.m`: gain scheduling PID feedback controller used to follow guidance reference trajectory
* `post_ascent_propagation.m`: dynamic equations used to propagate trajectory forward in time after final circularization burn is complete
* `drag.m`: atmospheric drag model of Martian atmosphere
* `plots.m`: produces convergence results and ascent trajectory output
* `control_subplots.m`: produces reference trajectory tracking from PID controller
* `control_animations.m`: produces reference trajectory tracking animations from PID controller
