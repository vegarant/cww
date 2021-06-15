% Computes the Walsh transform of the wavelet scaling function in each interval [l, l+1].
%
% Computes the Walsh transform of the wavelet scaling functions needed for the 
% fast computation of a matrix-vector multiplication with an N×M section of 
% the infinite-dimensional change of basis matrix between the Walsh basis and 
% a given wavelet basis.
% 
% The ratio N/M determines how many Walsh transform is computed. 
%
% Arguments
% ---------
% log2N (int): `ceil(log2(N))`, size of the output dimension.
% log2M (int): `ceil(log2(M))`, size of the input dimension.
% wname (str): Wavelet name
% bd_mode (str): Boundary handling (either 'per' or 'bd').
% j0 (int): Minimum wavelet decomposition level.
%
% Return
% ------
% pieces (struct): Walsh transform of the different wavelet pieces.
%
function pieces = cww_get_phi_walsh_pieces(log2N, log2M, wname, bd_mode, j0);

    r = log2N - log2M + 5;
    pieces = cww_get_phi_pieces(r, wname, bd_mode, j0);
    vm = cww_extract_vm_from_wname(wname);    
    is_per = cww_extract_is_per_from_bd_mode(bd_mode);

    if (is_per)
        for l = 1:2*vm -1
            pieces{l} = fastwht(pieces{l});
        end 
    else 

        for k = 1:vm
            for l = 1:vm+(k-1)
                pieces{k,l} = fastwht(pieces{k,l});
                pieces{vm+1+k,l} = fastwht(pieces{vm+1+k,l});
            end
        end

        for l = 1:2*vm-1
            pieces{vm+1,l} = fastwht(pieces{vm+1,l});
        end

    end
end


