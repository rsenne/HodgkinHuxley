% create beta function for variable h
function b1 = bh(v)
b1 = 1/(1+exp((30.0-v)/10.0));
end