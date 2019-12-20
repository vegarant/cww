%
% x     - Signal coefficients 
% mode  - 1: forward transform, 2: adjoint operation 
% vm    - Number of vanishing moments
% log2N - 2^log2N = N
% log2M - 2^log2M = M
%
function y = cww_kernel_CWW_per(x, mode, wname, log2N, log2M, phi_walsh_pieces)
    
    N = 2^log2N;
    M = 2^log2M;
    R = log2M;
    q = log2N - log2M;
    vm = cww_extract_vm_from_wname(wname);

    if (or(mode == 1 , strcmpi('notransp', mode)))
        
        c_sum = zeros(N,1);

        % The center wavelets
        for i = 1:2*vm-1
            l = i - vm;

            a = zeros(N,1);
            idx =  (vm-1+l:M-vm+l);
            a(1 + 2^q*idx) = x(1+idx-l);
            c = (N/sqrt(M))*fastwht(a);

            phi_walsh_l = phi_walsh_pieces{i}; 
            for j = 0:2^q-1
                idx = j*M + (1:M);
                c(idx) = phi_walsh_l(j+1).*c(idx);
            end
            c_sum = c_sum + c;
        end

        % Left 
        idx = 1:M;
        s = zeros([N,1]);
        for h = 0:2^q-1
            for k = 0:vm-2
                s1 = zeros([M,1]);
                for l = -vm+1:-k-1
                    s1 = s1 + had_mat_idx(N, M*h + idx, 2^q*(M+k+l)+1) * phi_walsh_pieces{l+vm}(1+h);
                end
                for l = -k:vm-1
                    s1 = s1 + had_mat_idx(N, M*h + idx, 2^q*(k+l)+1) * phi_walsh_pieces{l+vm}(1+h);
                end
                s(M*h + idx) = s(M*h + idx) + s1*x(k+1);
            end
        end
        c_sum_left  = s/sqrt(M);

        % Right
        c_sum_right = zeros(N,1);
        s = zeros([N,1]);
        for h = 0:2^q-1
            for k = M-vm+1:M-1
                s1 = zeros([M,1]);
                for l = -vm+1:M-k-1
                    s1 = s1 + had_mat_idx(N, M*h + idx, 2^q*(k+l)+1) * phi_walsh_pieces{l+vm}(1 + h);
                end
                for l = M-k:vm-1
                    s1 = s1 + had_mat_idx(N, M*h + idx, 2^q*(k+l-M)+1) * phi_walsh_pieces{l+vm}(1 + h);
                end
                s(M*h + idx) = s(M*h + idx) + s1*x(k+1);
            end
        end
        c_sum_right = s/sqrt(M);

        y = c_sum_left + c_sum + c_sum_right;

    else % mode ~= 1 or 'transp' 

        xl = cell(size(phi_walsh_pieces));
        for i = 1:2*vm-1
            phi_walsh_l = phi_walsh_pieces{i};
            x_copy = zeros(size(x));
            for n = 0:2^q-1 
                idx = M*n+(1:M);
                x_copy(idx) = x(idx)*phi_walsh_l(n+1);
            end
            xl{i} = (N/sqrt(M))*fastwht(x_copy);
        end

        % Center wavelets
        c = zeros(M,1);
        for l = -vm+1:vm-1
            i = l + vm;
            idx = vm:M-vm+1;
            c(idx) = c(idx) + xl{i}(1 + 2^q*(vm-1+l:M-vm+l));
        end

        % Left wavelets
        for m = 0:vm-2
            for l = -vm+1:-m-1
                i = l + vm;
                c(m+1) = c(m+1) + xl{i}( 1 + 2^q*(M+m+l) );
            end
            for l = -m:vm-1
                i = l + vm;
                c(m+1) = c(m+1) + xl{i}( 1 + 2^q*(m+l) );
            end
        end

        % Right wavelets
        for m = M-vm+1:M-1
            for l = -vm+1:M-m-1
                i = l + vm;
                c(m+1) = c(m+1) + xl{i}( 1 + 2^q*(m+l) );
            end
            for l = M-m:vm-1
                i = l + vm;
                c(m+1) = c(m+1) + xl{i}( 1 + 2^q*(m+l-M) );
            end
        end
        y = c;
    end
end
