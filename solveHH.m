%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
This script will solve and simulate the Hodgkin-Huxley equations first 
described in:
HODGKIN AL, HUXLEY AF. A quantitative description of membrane current and 
its application to conduction and excitation in nerve. J Physiol. 1952;117
(4):500-544. doi:10.1113/jphysiol.1952.sp004764. 

For this project I will compare a self-written Runge-Kutta method of 
numerical integration with Matlab's popular ODE45() function for solving 
ordinary differential equations.

Author:Ryan Senne
Date: 2/16/2021
Class:NE212 SPR 2021
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare workspace
clear
close all

%define paramters to be used in runge-kutta and hh
restingPotential = 0; %resting potential of neuron, using OG HH parameter
dt = 0.001; 
%{ 
time step by which we are varying by (ms), be aware that
the highest order of magnitude  you should use is 10^-3 (microseconds, 
because base unit is ms here). This is due to an integration problem caused 
by an improper step size. You will get NaN values if you ignore this.
%}
time = [0:dt:50]; % array of time points (ms)

%create pre-allocated arrays of each of our variables, i don't really need
%these but it's good practice
v = zeros(1,max(time));
m = zeros(1,max(time));
n = zeros(1,max(time));
h = zeros(1,max(time));

%create initial values of m, n, and h, membrane potential was previously
%defined
v(1) = restingPotential;
m(1) = am(v(1))/(am(v(1))+bm(v(1)));
n(1) = an(v(1))/(an(v(1))+bn(v(1)));
h(1) = ah(v(1))/(ah(v(1))+bh(v(1)));

%{
This is the interesting part of the script, lets initialize a runge-kutta
method for numerical integration to solve this ODE, I will be using the
most familiar version of this known as "RK4" method. You can refer to 
https://en.wikipedia.org/wiki/Runge–Kutta_methods for a more detailed
explanation.
The output of the hodgkinHuxley equation gives a 4 by 1 column vector such
that dydt = [v, n, m, v] thus kn(1,1) = vn, kn(2,1) = nn, kn(3,1) = mn,
kn(4,1) = hn
%}
  
for i=1:length(time)-1 % initialize for loop to run through time points
    % fist step of the runge-kutta method, slope at beginning of interval
    k1 = dt*hodgkinHuxley(i,[v(i); n(i); m(i); h(i)]);
    % step 2 of runge-kutta, slope at midpoint of interval
    k2 = dt*hodgkinHuxley((i+(0.5*dt)),[v(i)+(0.5*k1(1,1)); 
        n(i)+(0.5*k1(2,1)); m(i)+(0.5*k1(3,1)); h(i)+(0.5*k1(4,1))]);
    % step 3 of runge-kutta, slope at midpoint of interval
    k3 = dt*hodgkinHuxley(i+(0.5*dt),[v(i)+(0.5*k2(1,1)); 
        n(i)+(0.5*k2(2,1)); m(i)+(0.5*k2(3,1)); h(i)+(0.5*k2(4,1))]);
    % step 4 of runge-kutta, slope at end of interval
    k4 = dt*hodgkinHuxley(i+dt,[v(i)+k3(1,1); n(i)+k3(2,1); m(i)+k3(3,1); 
        h(i)+k3(4,1)]);
    
    % define the next step for each variable v,n,m,h
    v(i+1)= v(i)+(1/6)*(k1(1,1)+2*k2(1,1)+2*k3(1,1)+k4(1,1));
    n(i+1)= n(i)+(1/6)*(k1(2,1)+2*k2(2,1)+2*k3(2,1)+k4(2,1));
    m(i+1)= m(i)+(1/6)*(k1(3,1)+2*k2(3,1)+2*k3(3,1)+k4(3,1));
    h(i+1)= h(i)+(1/6)*(k1(4,1)+2*k2(4,1)+2*k3(4,1)+k4(4,1));

end

%store values of each variable at each iteration
finalV = v;
finalN = n;
finalM = m;
finalH = h;

clear v n m h;

% now lets use matlab's built-in ODE45 method (boring but effective)
% set initial values by use of ODE45
v=restingPotential; 
m=am(v)/(am(v)+bm(v)); 
n=an(v)/(an(v)+bn(v)); 
h=ah(v)/(ah(v)+bh(v)); 

vars = [v;n;m;h]; % vector of variables values
timeSpan = [0,max(time)]; % vector to store time values

[t, V] = ode45(@hodgkinHuxley,timeSpan,vars);

ODE45v=V(:,1);
ODE45n=V(:,2);
ODE45m=V(:,3);
ODE45h=V(:,4);
 
clear v;

%let's create some plots so we can see how our two methods compared

%plot gating variables as a function of time
figure
plot(time,finalN,'b', time, finalM,'r', time,finalH,'g',t,ODE45n,'b--', ...
    t, ODE45m,'r--', t,ODE45h,'g--' );
title('Gating Variable: Hodgkin-Huxley')
xlabel('Time (ms)')
ylabel('Gating Variable')
legend('Runge-Kutta N', 'Runge-Kutta M', 'Runge-Kutta H', 'ODE45 N',...
    'ODE45 M', 'ODE45 H')

%plot voltage vs time
figure
plot(time, finalV, 'r', t, ODE45v, 'b')
ylabel('Voltage (mV)')
xlabel('Time (ms)')
title('Voltage vs. Time')
legend('Runge-Kutta Voltage', 'ODE45 Voltage')

%plot gating variables as a function of voltage
figure
plot(finalV, finalN, 'b', finalV, finalM, 'r', finalV, finalH, 'g',...
    ODE45v, ODE45n, 'b--',ODE45v, ODE45m,'r--', ODE45v, ODE45h, 'g--' );
title('Gating Variable vs. Voltage')
xlabel('Voltage (mV)')
ylabel('Gating Variable')
legend('Runge-Kutta N', 'Runge-Kutta M', 'Runge-Kutta H', 'ODE45 N',...
    'ODE45 M', 'ODE45 H')


