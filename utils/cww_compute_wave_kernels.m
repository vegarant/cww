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
