% Computes a matrix-vector multiplication with an N^2 × M^2 section of a 
% change-of-basis matrix between a two-dimensional Walsh basis and a 
% two-dimensional orthonormal wavelet basis.
%
% Let `U` denote the change-of-basis matrix between a two-dimensional Walsh basis 
% and a two-dimensional orthonormal wavelet basis. Let `P_{M^2}` 
% denote a projection matrix onto the first M^2 components of a sequence. This 
% function computes the matrix-vector multiplication with the matrix 
%                                  P_{N^2} U P_{M^2}  
% and its transpose. 
%
% Arguments
% ---------
% x (vector): Input vector for the matrix-vector multiplication (a vectorized matrix).
% mode (int): If mode == 1 or mode == 'notransp', we compute a matrix-vector 
%             multiplication with `P_N U P_M`, otherwise we compute matrix-vector 
%             with its transpose.
% log2N (int): `N=2^(log2N)`, where N × Nis the number of rows
% log2M (int): `M=2^(log2M)`, where M × M is the number of columns.
% j0 (int): Minimum wavelet decomposition level.
% phi_walsh_pieces (struct): Walsh transform of the wavelet scaling function
%
% Return
% ------
% y (vector): The result of the matrix-vector multiplication with 
%             `P_{N}^2 U P_{M^2}` or its transpose.
%
function y = cww_handle_2d(x, mode, log2N, log2M, wname, bd_mode, j0, phi_walsh_pieces)

    vm = cww_extract_vm_from_wname(wname);
    is_per = cww_extract_is_per_from_bd_mode(bd_mode);

    if (nargin < 7) 
        j0 = cww_compute_j0(vm);
    end
    if j0 < round(log2(2*vm))
        disp('illegal j0');   
    end

    if (nargin < 8)
        phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
    end

    N = 2^log2N;
    M = 2^log2M;

    if is_per
        cww_kernel = @(x, mode) cww_kernel_CWW_per(x, mode, wname, log2N, log2M, phi_walsh_pieces);
    else 
        cww_kernel = @(x, mode) cww_kernel_CWW_bd(x, mode, wname, log2N, log2M, phi_walsh_pieces);
    end

    if (or(mode == 1, strcmpi(mode, 'notransp')))

        sc = reshape(x,M,M);

        Y = zeros(N,N);

        %parfor i = 1:M
        for i = 1:M
            Y(:, i) = cww_kernel(sc(:, i), mode);
        end
        X = Y;
        
        %parfor i = 1:N
        for i = 1:N
            Y(i,:) = cww_kernel(X(i, 1:M)', mode);
        end
        Y = Y';
        y = reshape(Y,N*N,1);

    else

        % x is N*N vector
        x = reshape(x,N,N);
        Y_tmp = zeros(N, M);

        %parfor i = 1:N
        for i = 1:N
            Y_tmp(i, :) = cww_kernel(x(i,:)' , mode);
        end

        Y = zeros(M,M);
        %parfor i = 1:M
        for i = 1:M
            Y(:, i) = cww_kernel(Y_tmp(:,i) , mode);
        end
        Y = Y';

        y = reshape(Y,M*M,1);

    end

end

