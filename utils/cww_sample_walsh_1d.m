% Computes the Walsh transform of a function f: [0,1) -> R.
%
% Given a function f: [0,1) -> R, this procedure computes the integral
%            Wf(n) = int_{0}^{1} f(t)w_n(t) dt,
% where w_n : [0,1) -> {-1,+1} is sequency ordered Walsh function, for 
% n=0,...,N-1 
%
% Arguments
% ---------
% f (function_handle): Function f : [0,1) -> R. 
% N (int): N-1 is the maxiumum Walsh frequency beeing sampled.
% r (int): Numerical integration uses 2^r sampling points in each interval of
%         length 1/N (optional, default r=4).
%
% Return
% ------
% samples (vector): The N first Walsh samples.
%
function samples = cww_sample_walsh_1d(f, N, r)
    
    if nargin < 3
        r = 4;
    end
    int_factor = 2^r;
    
    Nf = int_factor*N;

    t = linspace(1/(2*Nf), 1-(1/(2*Nf)), Nf)';

    func_values = f(t);
    samples = fastwht(func_values);    
    samples = samples(1:N);
end
