function X = cww_visualize_1d_samp_patt(N, idx)

    X = zeros([round(0.15*N), N], 'uint8');
    X(:,idx) = uint8(255);

end
