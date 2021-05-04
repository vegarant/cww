function y = cww_findim_handle_1d_cs(x, mode, idx, N, dwt_kernel, idwt_kernel);
    
    if isrow(x)
        x = x.';
    end
    
    if (or(mode == 1 , strcmpi('notransp', mode)))
        x = wl_idwt_impl_from_kernel(x, idwt_kernel);
        z = fastwht(x)*sqrt(N);
        y = z(idx);

    else % Transpose

        z = zeros([N, 1]);
        z(idx) = x;
        z = fastwht(z)*sqrt(N);
        y = wl_dwt_impl_from_kernel(z, dwt_kernel);
    
    end
    
    if isrow(y)
        y = y.';
    end

end
