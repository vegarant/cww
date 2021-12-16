% Computes a dense version of the matrix P_N U P_M, where U is the two dimensional 
% change-of-basis matrix between a Walsh basis and an orthonormal wavelet basis.
% The dense representation is computed using repeated calls to `cww_handle_2d`.
%
% Arguments
% ---------
% log2N (int): `N=2^(log2N)`, where N is the number of rows.
% log2M (int): `M=2^(log2M)`, where M is the number of columns.
% wname (str): Wavelet name.
% bd_mode (str): 'per' or 'bd'.  
% j0 (int): Minimum wavelet decomposition level.
%
% Return
% ------
% X (matrix): The N x M matrix  P_N U P_M. 
%
function X = cww_generate_full_matrix_2d(log2N, log2M, wname, bd_mode, j0)

    M = 2^log2M;
    N = 2^log2N;

    phi_walsh_pieces = cww_get_phi_walsh_pieces(log2N, log2M, wname, bd_mode, j0);

    G = @(x,mode) cww_handle_2d(x, mode, log2N, log2M, wname, bd_mode, j0, phi_walsh_pieces); 

    X = zeros(N*N,M*M);
    mode = 1;
    parfor i = 1:M*M
        ei = zeros(M*M,1);
        ei(i) = 1;
        col = G(ei, mode);
        X(:, i) = col;
    end


end

