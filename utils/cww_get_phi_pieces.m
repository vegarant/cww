% r - 2^r is the size of each pisce
% wname - Wavelet name
% bd_mode - 'per' or 'bd'. If 'per', pieces is a cell array with the 2*vm-1 unit pieces
% of phi. If 'bd', pieces is a 2*vm+1 times 2*vm-1, where the rows correpsonds
% to the different functions, and columns consists of the unit pieces. Notice
% that some of the entries will be empty due to the different supports.  
%
function pieces = cww_get_phi_pieces(r, wname, bd_mode, j0)

    dim = r + j0;
    N = 2^(r + j0);
    n = 2^r;

    vm = cww_extract_vm_from_wname(wname);
    is_per = cww_extract_is_per_from_bd_mode(bd_mode);

    nres = r;
    
    if (is_per)

        % Store current boundary extension mode
        %boundary_extension = dwtmode('status', 'nodisp');
        % Set the correct boundary extension for the wavelets
        %dwtmode('per', 'nodisp');


        %wname = sprintf('db%d', vm);
        %S = get_wavedec_s(dim, nres);

        phi = zeros([N,1]);
        phi(vm) = 2^((r)/2);    

        %phi = waverec(phi, S, wname);
        phi = wl_idwt_impl(phi, wname, 'm', nres, 'bd_mode', 'per', ...
                           'data_layout', 'resolution');
        pieces = cell(1,2*vm - 1);
        
        for k = 1:2*vm-1
            pieces{k} = phi( ((k-1)*2^r) + (1:2^r) );
        end

        % Restore dwtmode
%        dwtmode(boundary_extension, 'nodisp')

    else

        pieces = cell(2*vm + 1, 2*vm-1);

        for k = 1:vm

            ei_left = zeros([N,1]);
            ei_left(k) = 2^(r/2);
            wc_left = wl_idwt_impl(ei_left, wname, 'm', nres, 'bd_mode', 'bd', ...
                                   'prefilter_mode', 'bd_pre', 'data_layout', 'resolution');
            %wc_left = IWT_CDJV_noP(ei_left, j0, vm);

            ei_right = zeros([N,1]);
            ei_right(2^(j0)-k+1) = 2^(r/2);
            %wc_right = IWT_CDJV_noP(ei_right, j0, vm);
            wc_right = wl_idwt_impl(ei_right, wname, 'm', nres, 'bd_mode', 'bd', ...
                                   'prefilter_mode', 'bd_pre', 'data_layout', 'resolution');

            for i = 1:vm+k-1
                pieces{k,i} = wc_left( n*(i-1) + (1:n));
                pieces{k+vm+1,i} = wc_right( N-n*(vm+k-1) + n*(i-1) + (1:n));
            end 
        end

        % Ensure that there are interior wavelets
        if abs(round(log2(2*vm)) - j0) < 1e-5
            N = 2*N;
        end
        
        phi = zeros([N,1]);
        phi(vm) = 2^((r)/2);    

        %phi = waverec(phi, S, wname);
        phi = wl_idwt_impl(phi, wname, 'm', nres, 'bd_mode', 'per', ...
                           'data_layout', 'resolution');

        for i = 1:2*vm-1
            pieces{vm+1, i} = phi( ((i-1)*2^r) + (1:2^r) );
        end

    end

end




%% Example program used for debugging

% vm = 2;
% R = 5;
% q = 4;
% r = 4;
% is_per = 0;
% 
% 
% pieces = cww_get_phi_pieces(vm, r, is_per);
% 
% x_left = zeros((vm+k-1)*2^r,1);
% x_right = zeros((vm+k-1)*2^r,1);
% for l = 1:vm+(k-1)
%     x_left(idx + 2^r*(l-1)) =  pieces{k,l}; 
%     x_right(idx + 2^r*(l-1)) =  pieces{vm+1+k,l}; 
% end
% 
% subplot(1,2,1);
% plot(x_left)
% title('left');
% subplot(1,2,2);
% plot(x_right)
% title('right');


 
