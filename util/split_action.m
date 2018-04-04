function [ S0, S1, S2 ] = split_action( actions, alpha, gamma, sx, J, xkinks, Nt, Nx )
%SPLIT_ACTION split an array of actions into their contributions
%   action = S0*J + S1*gamma + S2*alpha
%   S0 = Nx*Nt*(2*xkinks - 1)
%   S1 = Nx*Nt*(2*sx - 1)
S0 = -Nx*Nt*(ones(size(xkinks)) - 2.0*xkinks);
S1 = -Nx*Nt*(ones(size(sx)) - 2.0*sx);
S2 = 1.0/alpha * (actions - gamma*S1 - J*S0);


end

