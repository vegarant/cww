% Computes a matrix-vector multiplication with an N × M section of a 
% change-of-basis matrix between a Walsh basis and an orthonormal wavelet basis.
%
% Let `U` denote the change-of-basis matrix between a Walsh basis and an 
% orthonormal wavelet basis.  Let `P_M` denote a projection matrix onto the 
% first M components of a sequence. This function computes the matrix-vector 
% multiplication with the matrix 
%                                  P_N U P_M  
% and its transpose. 
%
% Arguments
% ---------
% x (vector): Input vector for the matrix-vector multiplication
% mode (int): If mode == 1 or mode == 'notransp', we compute a matrix-vector 
%             multiplication with `P_N U P_M`, otherwise we compute matrix-vector 
%             with its transpose.
% log2N (int): `N=2^(log2N)`, where N is the number of rows
% log2M (int): `M=2^(log2M)`, where M is the number of columns.
% j0 (int): Minimum wavelet decomposition level.
% phi_walsh_pieces (struct): Walsh transform of the wavelet scaling function
%
% Return
% ------
% y (vector): The result of the matrix-vector multiplication with `P_N U P_M` 
%             or its transpose.
%
function y = cww_handle_1d(x, mode, log2N, log2M, wname, bd_mode, j0, phi_walsh_pieces)

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

    if is_per
        if (or(mode == 1 , strcmpi('notransp', mode)))

            %x = IWT_CDJV_noP(x, j0, vm);
            y = cww_kernel_CWW_per(x, mode, wname, log2N, log2M, phi_walsh_pieces);

        else % mode ~= 1 or 'transp' 

            y = cww_kernel_CWW_per(x, mode, wname, log2N, log2M, phi_walsh_pieces);
            %y = FWT_CDJV_noP(c, j0, vm);

        end
    else
        if (or(mode == 1 , strcmpi('notransp', mode)))

            %x = IWT_CDJV_noP(x, j0, vm);
            y = cww_kernel_CWW_bd(x, mode, wname, log2N, log2M, phi_walsh_pieces);

        else % mode ~= 1 or 'transp' 

            y = cww_kernel_CWW_bd(x, mode, wname, log2N, log2M, phi_walsh_pieces);
            %y = FWT_CDJV_noP(c, j0, vm);

        end
    end
end
