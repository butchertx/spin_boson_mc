function gamma = helicity(int, corr, type, angle, dim)
%%Calculate the helicity modulus given interactions, correlation function,
%%type of modulus, angle, and dimensions
%%int and corr should be 2-d matrices with the same number of elements along each dimension
%%angle should be a 2-component vector giving the angle in the x and y directions
%
%%type options: 'discrete' gives the energy gap for a discrete angle
%%'cont' gives the modulus defined as the second derivative of free energy wrt angle, at angle = 0

M = dim(1);
N = dim(2);
base_action = 0.5*M*N*(int'*corr);


if(size(int) ~= size(corr))
    disp('Error: Interactions and correlation sizes do not match')
    disp(size(int))
    disp(size(corr))
    return
end

if(type == 'disc')
    xangles = exp(1i*angle(1)/M * linspace(0,M - 1, M));
    yangles = exp(1i*angle(2)/N * linspace(0,N - 1, N));
    corr_tilde = kron(xangles,yangles).' .* corr;
    gamma = 0.5*M*N*(int'*corr_tilde) - base_action;
elseif(type == 'cont')
    if(angle(1) > 0)
        %twist along the x direction
        moment = 1/M*(0:1:(M-1));
        filler = ones(N,1);
        moment = kron(moment',filler);
    elseif(angle(2) > 0)
        %twist along the y direction
        moment = 1/N*(0:1:(N-1));
        filler = ones(M,1);
        moment = kron(filler, moment');
    else
        disp('Error: Must have a nonzero angle component to choose direction for modulus')
        return
    end
    gamma = -sum(int.*moment.*moment.*corr);
end

end
