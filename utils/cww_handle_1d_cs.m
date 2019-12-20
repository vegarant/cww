function y = cww_handle_1d_cs_internal(x, mode, idx, log2N, log2M, dwt_kernel, idwt_kernel, phi_walsh_pieces)
    
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

            %x = IWT_CDJV_noP(x, j0, vm);
            x = wl_idwt_impl_from_kernel(x, idwt_kernel);
            y = cww_kernel_CWW_per(x, mode, wname, log2N, log2M, phi_walsh_pieces);
            y = y(idx);
            
        else % mode ~= 1 or 'transp' 

            z = zeros([N,1]);
            z(idx) = x;
            c = cww_kernel_CWW_per(z, mode, wname, log2N, log2M, phi_walsh_pieces);
            y = wl_dwt_impl_from_kernel(c, dwt_kernel);

        end
    else
        if (or(mode == 1 , strcmpi('notransp', mode)))

            %x = IWT_CDJV_noP(x, j0, vm);
            x = wl_idwt_impl_from_kernel(x, idwt_kernel);
            y = cww_kernel_CWW_bd(x, mode, wname, log2N, log2M, phi_walsh_pieces);
            y = y(idx);

        else % mode ~= 1 or 'transp' 

            z = zeros([N,1]);
            z(idx) = x;
            c = cww_kernel_CWW_bd(z, mode, wname, log2N, log2M, phi_walsh_pieces);
            y = wl_dwt_impl_from_kernel(c, dwt_kernel);
            %y = FWT_CDJV_noP(c, j0, vm);

        end
    
    end

