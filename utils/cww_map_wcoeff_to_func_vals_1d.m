% Evaluates a linear combination of scaling functions at dyadic grid points given 
% their coefficients.
% 
% Evaluates the function 
%             f(t) = sc(0)*phi_{j,0}(t) + ... + sc(M-1)*phi_{j,M-1}(t)
% at the grid points `t = (0:N-1)/N`.
%
% Arguments
% ---------
% sc (vector): scaling function coefficents.
% log2N (int): `N=2^(log2N)`. Number of grid points.
% wname (str): Wavelet name.
% bd_mode (str): Boundary handling (either 'bd' or 'per').
%
% Return 
% ------
% x (vector): `f` evaluated at the grid points.  
function x = cww_map_wcoeff_to_func_vals_1d(sc, log2N, wname, bd_mode)

    M = length(sc);
    log2M = round(log2(M));

    if abs(2^log2M - M) > 1e-2
        disp('Error: cww_map_wcoeff_to_func_vals, length of sc, should be 2^m,');
        disp('       for an integer m.');
    end

    if isrow(sc)
        sc = sc';
    end

    T = cww_get_scaling_matrix(log2N, log2M, wname, bd_mode);
    x = T*sc;   

end











