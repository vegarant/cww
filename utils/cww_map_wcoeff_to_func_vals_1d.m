function x = cww_map_wcoeff_to_func_vals_1d(sc, log2N, wname, bd_mode)

    M = length(sc);
    log2M = round(log2(M));

    if abs(2^log2M - M) > 1e-2
        disp('Error: cww_map_wcoeff_to_func_vals, length of sc, should be 2^m,');
        disp('       for an integer m.');
    end

    if isrow(sc)
        sc = sc';
    end

    T = cww_get_scaling_matrix(log2N, log2M, wname, bd_mode);
    x = T*sc;   

end











