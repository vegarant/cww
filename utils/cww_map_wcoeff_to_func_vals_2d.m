function X = cww_map_wcoeff_to_func_vals_2d(sc, log2N, wname, bd_mode)

    [M1, M2] = size(sc);

    if M1 ~= M2
        disp('Error: cww_map_wcoeff_to_func_vals, sc should be a square matrix');
    end

    log2M = round(log2(M1));

    if abs(2^log2M - M1) > 1e-2
        disp('Error: cww_map_wcoeff_to_func_vals_2d, length of sc, should be 2^m,');
        disp('       for an integer m.');
    end

    T = cww_get_scaling_matrix(log2N, log2M, wname, bd_mode);

    M = 2^log2M;
    N = 2^log2N;

    X_tmp = T*sc;
    X = T*(X_tmp');

end

