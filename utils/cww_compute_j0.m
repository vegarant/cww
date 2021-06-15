% Computes the minimum wavelet decomposition level j0
%
% Given an orthonormal wavelet with `vm` vanishing moments, this function 
% computes the coarsest wavelet resolution level j0 for this wavelet. 
%
% Arguments
% ---------
% vm (int): Number of vanishing moments
%
% Return
% ------
% j0 (int): Minimum wavelet decomposition level
%
function j0 = cww_compute_j0(vm)
    j0 = 0;
    while (2^j0 <= 2*vm) % This is not optimal in all cases, but it is 
                         % necessary in order to make a unified interface with 
                         % the non-boundary wavelets. Do not change it! 
        j0 = j0 + 1;
    end
end


