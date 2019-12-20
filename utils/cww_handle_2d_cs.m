function y = cww_handle_2d_cs_internal(x, mode, idx, log2N, log2M, dwt_kernel, idwt_kernel, phi_walsh_pieces)

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
