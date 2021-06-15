% Computes a matrix-vector multiplication between a subsampled two-dimensional 
% Walsh sampling basis and a two-dimensional orthonormal wavelet reconstruction 
% basis.
%
% Let `U` denote the change-of-basis matrix between a two-dimensional Walsh basis 
% and a two-dimensional orthonormal wavelet basis. Let `P_{M^2}` 
% denote a projection matrix onto the first M^2 components of a sequence and 
% let `P_{omage}` be the projection extracting the vector elements enumerated 
% by the set `Omega \subset {1, ..., N^2}`.
%
% Let `A = P_{Omega} U P_{M^2}`. This function computes the matrix-vector 
% multiplication with `A` and `A'`.
%
% Arguments
% ---------
% x (vector): Input vector for the matrix-vector multiplication (vectorized matrix)
% mode (int): If mode == 1 or mode == 'notransp', we compute a matrix-vector 
%             multiplication with `A`, otherwise we compute matrix-vector with 
%             its transpose.
% idx (vector): This is the set `Omega` above. It is the indices of the samples 
%               we want to acquire (linear ordering). 
% log2N (int): `N=2^(log2N)`, where N^2 is an upper bound on the elements in `idx`.
% log2M (int): `M=2^(log2M)`, where M × M is the number of columns of `A`.
% dwt_kernel (struct): DWT kernel, used to compute the discrete wavelet transform.  
% idwt_kernel (struct): IDWT kernel, used to compute the inverse discrete wavelet 
%                       transform.  
% phi_walsh_pieces (struct): Walsh transform of the wavelet scaling function
%
% Return
% ------
% y (vector): The result of the matrix-vector multiplication with `A` or its 
%             transpose.
%
function y = cww_handle_2d_cs(x, mode, idx, log2N, log2M, dwt_kernel, idwt_kernel, phi_walsh_pieces)

    wname   = dwt_kernel.wave_name;
    bd_mode = dwt_kernel.bd_mode;
    nres    = dwt_kernel.m;

    vm = cww_extract_vm_from_wname(wname);
    is_per = cww_extract_is_per_from_bd_mode(bd_mode);

    M = 2^log2M;
    N = 2^log2N;
    j0 = log2M-nres;

    if isrow(idx);
        idx = idx';
    end
    if isrow(x)
        x = x';
    end

    if is_per
        if (or(mode == 1 , strcmpi('notransp', mode)))

            x = reshape(x, M, M);
            x = wl_idwt_impl_from_kernel(x, idwt_kernel);
            y = cww_handle_2d(x, mode, log2N, log2M, wname, bd_mode, j0, phi_walsh_pieces);
            y = y(idx);

        else % mode ~= 1 or 'transp' 

            z = zeros([N, N]);
            z(idx) = x;
            c = cww_handle_2d(z, mode, log2N, log2M, wname, bd_mode, j0, phi_walsh_pieces);
            y = wl_dwt_impl_from_kernel(reshape(c, M, M), dwt_kernel);
            y = y(:);

        end
    else
        if (or(mode == 1 , strcmpi('notransp', mode)))

            x = reshape(x, M, M);
            x = wl_idwt_impl_from_kernel(x, idwt_kernel);
            y = cww_handle_2d(x, mode, log2N, log2M, wname, bd_mode, j0, phi_walsh_pieces);
            y = y(idx);

        else % mode ~= 1 or 'transp' 

            z = zeros([N, N]);
            z(idx) = x;
            c = cww_handle_2d(z, mode, log2N, log2M, wname, bd_mode, j0, phi_walsh_pieces);
            y = wl_dwt_impl_from_kernel(reshape(c, M, M), dwt_kernel);
            y = y(:);

        end

end
