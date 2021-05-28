function [idx, str_id] = cww_sph1_power_law(N, nbr_samples, alpha, r_0)

    str_id = sprintf('power_law_1d_al_%g', alpha);

    if nargin == 4
        str_id = sprintf('power_law_1d_al_%g_r0_%d', alpha, r_0);
    elseif nargin == 3
        r_0 = 0;
        str_id = sprintf('power_law_1d_al_%g', alpha);
    elseif nargin == 2
        r_0 = 0;
        alpha = 1;
        str_id = 'power_law_1d';
    end

    if r_0 == 0;

        X = (2:N)';
        P = 1./((X).^(alpha));
        P = P/sum(P(:));

        idx = datasample(2:N, nbr_samples-1,'Replace', false, 'Weights', P);
        idx = sort(idx');
        idx = [1; idx]; % Always sample zero frequency

    else 

        N1 = 2^r_0;
        X = (1:(N-N1))';
        P = 1./((X).^(alpha));
        P = P/sum(P(:));
        idx = datasample(N1+1:N, nbr_samples-N1,'Replace', false, 'Weights', P);
        idx = [(1:N1)' ; idx'];

    end

end 

