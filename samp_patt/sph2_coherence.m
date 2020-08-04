% For a N×N image, it samples `nbr_samples` according to the coherence density. 
% The coherence structure is competed for the change of basis matrix between 
% Walsh functions and Dauberchies wavelets with more than two vanishing moments.
% Note that the densities will be different for the Haar wavelet. 
%
% INPUT:
% N - Size of sampling pattern (N×N).
% nbr_samples - Total number of samples  
%
% OUTPUT
% idx - linear indices of the chosen samples
function idx = sph2_coherence(N, nbr_samples)
    vm = 3;
    r = round(log2(N));
    J0 = computeJ0(vm);
    nres = r-J0;
    
    bv = [];
    for i = 1:nres
        bv =[bv, 2^(J0+i)];
    end
    bh = bv;
    
    density = compute_density(r, J0, vm);
    idx = sph2_multi(nbr_samples, density, bv, bh);

end

% Computes the densities for a given wavelet.
%
% INPUT:
% r - Number of levels. Sampling levels will be M = [2^(J_0+1), ..., 2^(J_0+r)].
% J_0 - Minimum wavelet decomposition level
% vm  - Number of vanishing moments
%
% OUTPUT:
% density - Density matrix.
function density = compute_density(r, J_0, vm)
    nres = r-J_0;
    N = 2^r;
    wname = sprintf('db%d', vm);
    bv = [];
    for i = 1:nres
        bv =[bv, 2^(J_0+i)];
    end
    bh = bv;

    im = phantom(N);
    
    [c,S] = wavedec2(im, nres, wname);

    M = [1];
    
    for i = 1:nres
        M = [M, 2^(2*(J_0+i))]; 
    end
    
    
    sparsity = [];
    
    for i = 1:nres
        idx = M(i):M(i+1);
        s = sum( abs( c( idx ) ) > 1e-12); 
        sparsity = [sparsity, s];
    end
    
    m_mat = zeros([nres, nres]);
    
    for k1=1:nres
        for k2=1:nres
            m = 0;
            if ( k1 > 1 & k2 > 1)
    
                for l = 1:nres
                    m = m + sparsity(l)*2^( -abs(l-k1) - abs(l-k2) );
                end
    
            elseif (k1 == 1 & k2 > 1)
    
                for l = 1:nres
                    m = m + 4*sparsity(l)*2^( -l - abs(l-k2) );
                end
    
            elseif (k1 > 1 & k2 == 1)
            
                for l = 1:nres
                    m = m + 4*sparsity(l)*2^( -l - abs(l-k1) );
                end
            
            elseif (k1 == 1 & k2 == 1)
    
                for l = 1:nres
                    m = m + 16*sparsity(l)*2^( -2*l);
                end
            
            else
                disp('Something in my logic is wrong')
            end
            m_mat(k1,k2) = m;
        end
    end
    
    a = bv(1:end) - [0, bv(1:end-1)];
    b = bh(1:end) - [0, bh(1:end-1)];
    max_level_size = a'*b;
    
    density = m_mat./max_level_size;
    density = density/sum(m_mat(:));


end
