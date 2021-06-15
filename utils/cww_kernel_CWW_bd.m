% Internal function used to compute the matrix-vector multiplication (for wavelets
% with vanishing moments preserving wavelets). 
%
% Let `U` denote the change-of-basis matrix between a Walsh basis and an 
% orthonormal wavelet basis with vanishing moments preserving wavelets at the 
% boundaries. Let `P_M` denote a projection matrix onto the first M
% components of a sequence. This function computes the matrix-vector 
% multiplication with the matrix 
%                                  P_N U P_M  
% and its transpose. 
%
% Arguments
% ---------
% x (vector): Input vector for the matrix-vector multiplication
% mode (int): If mode == 1 or mode == 'notransp', we compute a matrix-vector 
%             multiplication with `P_N U P_M`, otherwise we compute matrix-vector 
%             with its transpose.
% wname (str): Wavelet name
% log2N (int): `N=2^(log2N)`, where N is the number of rows
% log2M (int): `M=2^(log2M)`, where M is the number of columns.
% phi_walsh_pieces (struct): Walsh transform of the wavelet scaling function
%
% Return
% ------
% y (vector): The result of the matrix-vector multiplication with `P_N U P_M` 
%             or its transpose.
%
function y = cww_kernel_CWW_bd(x, mode, wname, log2N, log2M, phi_walsh_pieces)

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
            idx =  (vm+l:M-vm-1+l);
            a(1 + 2^q*idx) = x(1+idx-l);
            c = (N/sqrt(M))*fastwht(a);

            phi_walsh_l = phi_walsh_pieces{vm+1,i}; 
            for h = 0:2^q-1
                idx = h*M + (1:M);
                c(idx) = phi_walsh_l(h+1).*c(idx);
            end
            c_sum = c_sum + c;
        end

        % Left 
        idx = 1:M;
        s = zeros([N,1]);
        for h = 0:2^q-1
            for m = 0:vm-1
                s1 = zeros([M,1]);
                for l = 0:vm+m-1
                    s1 = s1 + had_mat_idx(N, M*h  + idx, (2^q)*l+1) * phi_walsh_pieces{m+1, l+1}(1+h);
                end
                s(M*h + idx) = s(M*h + idx) + s1*x(m+1);
            end
        end
        c_sum_left  = s/sqrt(M);

        % Right
        c_sum_right = zeros(N,1);
        s = zeros([N,1]);
        for h = 0:2^q-1
            for k = 0:vm-1
                s1 = zeros([M,1]);
                for l = 0:vm-1+k 
                    s1 = s1 + had_mat_idx(N, M*h + idx, 2^q*(M-vm-k+l)+1) * phi_walsh_pieces{vm+2+k, l+1}(1 + h);
                end
                s(M*h + idx) = s(M*h + idx) + s1*x(M-k);
            end
        end
        c_sum_right = s/sqrt(M);

        y = c_sum_left + c_sum + c_sum_right;

    else % mode ~= 1 or 'transp' 

        xl = cell(2*vm+1, 2*vm-1);
        % Boundary vectors
        idx = 1:M;
        for m = 0:vm-1
            for l = 0:vm-1+m
                x_copy_left  = zeros([N,1]);
                x_copy_right = zeros([N,1]);
                for h = 0:2^q-1
                    x_copy_left( M*h + idx) = x(M*h + idx)*phi_walsh_pieces{m+1,l+1}(1+h);
                    x_copy_right(M*h + idx) = x(M*h + idx)*phi_walsh_pieces{vm+m+2,l+1}(1+h);
                end
                xl{m+1, l+1}    = (N/sqrt(M))*fastwht(x_copy_left);
                xl{vm+m+2, l+1} = (N/sqrt(M))*fastwht(x_copy_right);
            end
        end

        for l = -vm+1:vm-1
            x_copy = zeros([N,1]);
            for h=0:2^q-1
                x_copy( M*h + idx) = x(M*h + idx)*phi_walsh_pieces{vm+1,l+vm}(1+h);
            end
            xl{vm+1, l+vm} = (N/sqrt(M))*fastwht(x_copy);
        end


        % Center wavelets
        c = zeros(M,1);
        idx = vm+1:M-vm;
        for l = -vm+1:vm-1
            c(idx) = c(idx) + xl{vm+1, l+vm}(1 + 2^q*(vm+l:M-vm-1+l));
        end

        % Left wavelets
        for m = 0:vm-1
            for l = 0:vm-1+m
                c(m+1) = c(m+1) + xl{m+1,l+1}( 1 + (2^q)*l );
            end
        end
        
        for m = 0:vm-1
            for l = 0:vm-1+m
                c(M-m) = c(M-m) + xl{vm+2+m,l+1}( 1 + (2^q)*(M-vm-m+l) );
            end
        end
        y = c;
    end

end


