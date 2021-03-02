# Hodgkin-Huxley
This project satfies project 1 of 2 for NE212 SPR 2021.
In this project I simulate the Hodgkin-Huxley differential equations, I used a self-written Runge-Kutta (RK4), to solve the set of diffyQ's as well as the ODE45() function of Matlab. 
The answers to the project questions are described below.
1. Plot the properties of the gating variables m,n, and h as a function of input voltage. Does this match yourunderstanding of how an action potential is generated?
(The plot for this question can be obtained by running solveHH.mat)
2. Simulate an action potential using the hodgkin huxley model. Get it to generate a spike!
(The plot for this question can be obtained by running solveHH.mat)
3. Simulate action potentials using the hodgkin huxley model you programmed and also use several different currents from 2μA to 10μA.
(The plot for this question can be obtained by running solveHH.mat)
4. Label all axes and include legends and make sure all labels are legible. Provide a title for all the graphs atthe top of the figure. (This will need some research)
(The plot for this question can be obtained by running solveHH.mat)
5.All functions need to be commented and code needs appropriate comments. Files should execute withminimal modification
The only modifications possibly necessary for my scripts is to alter the time interval, dt, in solveHH.m, or the Input current in HodgkinHuxley.m
