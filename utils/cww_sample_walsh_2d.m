% Computes the Walsh transform of a function f: [0,1)^2 -> R.
%
% Given a function f: [0,1)^2 -> R, this procedure computes the integral
%   Wf(n) = int_{0}^{1} int_{0}^{1} f(t_1,t_2)w_{n_1}(t_1)w_{n_2}(t_2) dt_1 dt_2,
% where w_n : [0,1) -> {-1,+1} is sequency ordered Walsh function, for 
% n_1, n_2 = 0,...,N-1 
%
% Arguments
% ---------
% f (function_handle): Function f : [0,1)^2 -> R. 
% N (int): N-1 is the maxiumum Walsh frequency beeing sampled.
% r (int): Numerical integration uses 2^r sampling points in each interval of
%         length 1/N (optional, default r=4).
%
% Return
% ------
% samples (mat): The N×N first Walsh samples.
%
function samples = cww_sample_walsh_2d(f, N, r)

    if nargin < 3
        r = 4;
    end

    int_factor = 2^r;
    Nf = int_factor*N;

    t = linspace(1/(2*Nf), 1-(1/(2*Nf)), Nf);
    
    [X,Y] = meshgrid(t,t);
    F = f(X,Y);
    samples = cww_fastwht_2d(F);

    samples = samples(1:N, 1:N);

end

