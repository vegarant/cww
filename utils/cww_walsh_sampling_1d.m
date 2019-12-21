% Samples the function f: [0,1) -> R with walsh functions as the integral
% transform
%            Wf(n) = int_{0}^{1} f(x)w_n(x) dx
%
% Here N-1 is the maximum frequency and Omega is a subset of {0, ..., N-1} 
%
% INPUT
% f     - Function f : [0,1) -> R 
% N     - N-1 is the maximum frequency
% Omega - Subset of {0,...,N-1}, if no argument is provided Omega is chosen as 
%         Omega = {0, ..., N-1}.
%
% OUTPUT
% The walsh samples specified in Omega.
%
function samples = cww_walsh_sampling_1d(f, N, Omega)
    if (nargin < 3) 
        Omega  = 1:N;
    end

    int_factor = 2^4;
    Nf = int_factor*N;
    
    t = linspace(1/(2*Nf), 1-(1/(2*Nf)), Nf)';
    
    func_values = f(t);
    wfunc_values = fastwht(func_values);    
    samples = wfunc_values(Omega);

end
