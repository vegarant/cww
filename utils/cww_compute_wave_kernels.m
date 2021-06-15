% Computes the wavelet filter coefficents used by the wl library.
%
% A wavelet is given by its wavelet coefficients. This function extracts these
% coefficients and store them in wavelet kernels. These kernels have
% information about the boundary handling, dimension of the transform 
% and the minimum wavelet decomposition level. 
%
% This function aims to avoid repeated calls to the `wl_dwt_kernel`
% function since this function can be rather slow. Instead, this function allows 
% you to compute these kernels once and use them repeatedly without recomputing 
% them. Its intended usage is before calling the `cww_handle_1d_cs` and 
% `cww_handle_2d_cs` functions.
%
% Arguments
% ---------
% log2M (int): `ceil(log2(M))` where M^dims is the dimension of the input signal.
% dims (int): Dimension of the wavelet transform (dims = 1,2,3 are supported).
% wname (str): Wavelet name
% bd_mode (str): Boundary handling (either 'bd' or 'per').
% j0 (int): Minimum wavelet decomposition level.
%
% Return
% ------
% dwt_kernel (struct) : DWT kernel
% idwt_kernel (struct) : IDWT kernel
%
function [dwt_kernel, idwt_kernel] = cww_compute_wave_kernels(log2M, dims, wname, bd_mode, j0)

    nres = log2M - j0;
    is_per = cww_extract_is_per_from_bd_mode(bd_mode);

    M = 2^log2M;
    if dims == 1
        data_size = [M,1];
    elseif dims == 2
        data_size = [M, M];
    elseif dims == 3
        data_size = [M, M, M];
    else
        sprintf('Error: dims has to be 1, 2 or 3. dims=%d', dims);
    end

    if is_per
        dwt_kernel = wl_dwt_kernel(wname, data_size, dims, 'bd_mode', 'per', ...
                                     'm', nres, 'prefilter_mode', 'none');
        idwt_kernel = wl_idwt_kernel(wname, data_size, dims, 'bd_mode', 'per', ...
                                     'm', nres, 'prefilter_mode', 'none');
    else
        dwt_kernel = wl_dwt_kernel(wname, data_size, dims, 'bd_mode', 'bd', ...
                                     'm', nres, 'prefilter_mode', 'none');
        idwt_kernel = wl_idwt_kernel(wname, data_size, dims, 'bd_mode', 'bd', ...
                                     'm', nres, 'prefilter_mode', 'none');
    end

end
