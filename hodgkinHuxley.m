% define a function that describes the titular Hodgkin-Huxley equation,
% this fucntion is a set of four differential equations
function dydt = hodgkinHuxley(~, varSpan)
% this will set each of the starting values equivalent to the function input
v = varSpan(1);
n = varSpan(2);
m = varSpan(3);
h = varSpan(4);

I = 10.0; % current(microAmps) for spikes
Cm = 0.09; % capacitance microfarads

% define constants for Hodgkin-Huxley
gNa = 120; % mS/cm^2 sodium conductance
gK = 36; % mS/cm^2 potassium conductance
gL = 0.3; % mS/cm^2 leak conductance
voltageNa = 115.0; % mV reversal potential sodium
voltageK = -12.0; % mV reversal potential potassium
voltageL = -10.36; % mV reversal potential leak

dvdt = [(((-gNa*(m^3)*h)*(v-voltageNa))-((gK*(n^4))*(v-voltageK))-(gL*(v-voltageL))+I)/Cm];
dndt = [an(v)*(1-n)-bn(v)*n];
dmdt = [am(v)*(1-m)-bm(v)*m];
dhdt = [ah(v)*(1-h)-bh(v)*h];
dydt = [dvdt; dndt; dmdt; dhdt];
end
    