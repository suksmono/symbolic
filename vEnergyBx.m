% spin energy of the Ising model
function E=vEnergyBx(bi, hi, Jij, vSpin)
%% Calculate the energy of the spin vector system
E = bi + hi'*vSpin+ vSpin'*(Jij*vSpin) ;
