
R = 7;
q_values = [1,2,3,4];
is_per = 0;


if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

all_wnames = {'db2', 'db3', 'db4', 'db5', 'db6', 'sym2', 'sym3', 'sym4', 'sym5', 'sym6'};
nbr_wavelets = length(all_wnames);
nbr_q_values = length(q_values);

for i = 1:nbr_wavelets
    wname = all_wnames{i};
    vm = cww_extract_vm_from_wname(wname);
    fprintf('%s ', wname);
    for j = 1:nbr_q_values
        q = q_values(j);
        log2N = R+q;
        log2M = R;
        j0 = cww_compute_j0(vm);
        
        phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
        G = @(x,mode) cww_handle_1d(x, mode, R+q, R, wname, bd_mode, j0, phi_walsh_pieces); 
        M = 2^log2M;
        N = 2^log2N;
        sing_vals = svds(G, [N,M], M);
        fprintf('& %6.3f ', 1/sing_vals(M));
    end
    fprintf('\\\\ \n');
end


fprintf('\n\n Two dimensions \n\n');

R = 5;
% Note: These computations takes quite some time
for i = 1:nbr_wavelets
    wname = all_wnames{i};
    vm = cww_extract_vm_from_wname(wname);
    fprintf('%s ', wname);
    for j = 1:nbr_q_values
        q = q_values(j);
        log2N = R+q;
        log2M = R;
        j0 = cww_compute_j0(vm);
        
        phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
        G = @(x,mode) cww_handle_2d(x, mode, R+q, R, wname, bd_mode, j0, phi_walsh_pieces); 
        M = 2^log2M;
        N = 2^log2N;
        sing_vals = svds(G, [N*N,M*M], M*M);
        fprintf('& %6.3f ', 1/sing_vals(M*M));
    end
    fprintf('\\\\ \n');
end




