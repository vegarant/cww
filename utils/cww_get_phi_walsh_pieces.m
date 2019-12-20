% Update code

% vm    - Number of vanishing moments
% j     - scaling function decomposition level i.e. 2^j = M
% log2N - 2^log2N = N, number of samples
% is_per - whether a periodic or boundary wavelet basis should be used
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


