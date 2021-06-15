% Computes a one-dimensional sampling pattern for Walsh functions using a 
% power-law sampling pattern.
% 
% This function pick `nbr_samples` elements from the set {0, ..., N-1} without 
% replacement. Here the numbers {0, ..., 2^{r_0}-1} are picked deterministically, 
% while we pick `nbr_samples - 2^{r_0}` numbers from the set {r^{r_0}, ..., N-1}
% according to a probability distribution. 
%
% Let { pi_{k} }_{k=2^{r_0}}^{N-1} be a probability distribution on 
% {2^{r_0},..., N-1} given by
%               pi_k = C/(max{1,k}^(alpha))
% where C is a normalizing constant. The `nbr_samples - 2^{r_0}`,  
% numbers from the set {2^{r_0}, ..., N-1} are picked according to this 
% distribution. 
%
% Arguments
% ---------
% N (int): Size of the set we sample from.
% nbr_samples (int): Number of samples.
% alpha (float): Exponent in the power law (optional, default alpha=1).
% r_0 (int): The 2^{r_0} first frequencies are sampled deterministically 
%             (optional, default r_0=0)
% Return
% ------
% idx (vector): The `nbr_samples` we picked from the set {0, ..., N-1}. 
% str_id (str): String identifying the sampling pattern we used.
% 
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

