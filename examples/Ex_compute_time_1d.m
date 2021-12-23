% This script computes the execution time for the fast transform and dense 
% matrix-vector multiplication for different input sizes. The result is stored 
% in .mat files. Use the script `Ex_create_time_plot_1d.m` to plot the result.

clear all;
dwtmode('per', 'nodisp')

dest = 'data';
if (exist(dest) ~= 7) 
    mkdir(dest);
end

vm = 4;
R_range = 4:16;
q = 1;
j0 = 3;     % cww_compute_j0(vm);
is_per = 0; % Set this to 0, to use boundary corrected wavelets
wname = sprintf('db%d', vm);
if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

mode = 1;

n = length(R_range);

time_forward_G = zeros(n,1);
time_forward_X = zeros(n,1);

for i = 1:n
    R = R_range(i);
    fprintf('R: %d\n', R);
    M = 2^R;
    N = 2^(R+q);
    X = randn(N,M);%cww_generate_full_matrix_1d(R+q, R, wname, bd_mode, j0);
    phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
    G = @(x,mode) cww_handle_1d(x, mode, R+q, R, wname, bd_mode, j0, phi_walsh_pieces); 
    x = randn(M,1);
    
    f_G = @() G(x, 1);
    f_X = @() X*x;
    
    fprintf('Timing the functions\n');
    t_G = timeit(f_G);
    t_X = timeit(f_X);
    
    time_forward_G(i) = t_G;
    time_forward_X(i) = t_X;

end

fname = sprintf('time_1d_forward_R_%d_%d_q_%d_db%d_%s_j0_%d.mat', R_range(1), R_range(end), q, vm, bd_mode, j0);

save(fullfile(dest, fname), 'time_forward_G', 'time_forward_X');

time_adjoint_G = zeros(n,1);
time_adjoint_X = zeros(n,1);

for i = 1:n
    R = R_range(i);
    fprintf('R: %d\n', R);
    M = 2^R;
    N = 2^(R+q);
    X = randn(M, N);%cww_generate_full_matrix_1d(R+q, R, wname, bd_mode, j0);
    phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
    G = @(x,mode) cww_handle_1d(x, mode, R+q, R, wname, bd_mode, j0, phi_walsh_pieces); 
    x = randn(N,1);
    
    f_G = @() G(x, 0);
    f_X = @() X*x;
    
    fprintf('Timing the functions\n');
    t_G = timeit(f_G);
    t_X = timeit(f_X);
    
    time_adjoint_G(i) = t_G;
    time_adjoint_X(i) = t_X;

end

fname = sprintf('time_1d_adjoint_R_%d_%d_q_%d_db%d_%s_j0_%d.mat', R_range(1), R_range(end), q, vm, bd_mode, j0);

save(fullfile(dest, fname), 'time_adjoint_G', 'time_adjoint_X');

