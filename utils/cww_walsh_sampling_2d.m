% Samples the function f: [0,1)^2 -> R with walsh functions as the integral
% transform
%            Wf(n,m) = int_{0}^{1} int_{0}^{1} f(x,y)w_n(x)w_m(y) dx dy
%
% Here N-1 is the maximum frequency and Omega is a subset of {0, ..., N-1} 
%
% INPUT
% f     - Function f : [0,1)^2 -> R 
% N     - N-1 is the maximum frequency
% Omega - two columns with a subset of {0,...,N-1}, if no argument is provided Omega is chosen as 
%         Omega = [(0:N-1)', (0:N-1)'].
%
% OUTPUT
% The walsh samples specified in Omega.
%
function samples = cww_walsh_sampling_2d(f, N)

    Omega  = [(0:N-1)', (0:N-1)'];
    eps = 1e-14;
    samples = zeros(size(Omega,1), size(Omega,2));

    int_factor = 2^3;
    t = linspace(0, 1-eps, int_factor*N);
    [X,Y] = meshgrid(t,t);
    F = f(X,Y);
    samples = cww_fastwht_2d(F);

    samples = samples(1:N, 1:N);

end

