% Samples the function f: [0,1) -> R with walsh functions as the integral
% transform
%            Wf(n) = int_{0}^{1} f(x)w_n(x) dx
%
% Here N-1 is the maximum frequency.
%
% INPUT
% f     - Function f : [0,1) -> R 
% N     - N-1 is the maximum frequency
% r     - Numerical integration uses 2^r sampling points in each interval of
%         length 1/(2*N) 
%
% OUTPUT
% The N first walsh samples.
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
