%Spiking Neuron Simulation
%---------------------------------
%Stage1: simple model description
%C dV/dt = -gL*(V-VL) + I(t)
%solve using Runge kutta method
%---------------------------------

%Constants
gL=100;
VL=