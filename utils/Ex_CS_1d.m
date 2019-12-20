clear all
R = 6;
q = 1;
q2 = 3;
vm = 4;
is_per = 0;
subsampling_rate = 0.25;
noise = 0.01;
dims = 1;

M = 2^R;
N = 2^(R+q);
N2 = 2^(R+q+q2);
wname = sprintf('db%d', vm);
j0 = cww_compute_j0(vm);
log2N = R+q;
log2M = R;
nbr_samples = round(N*subsampling_rate); 
if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
[dwt_kernel, idwt_kernel] = cww_compute_wave_kernels(log2M, dims, wname, bd_mode, j0);

f = @(x) 50*(x-0.5).*(x-0.05).*(x-0.95) + x;%cos(2*pi*x)  + 0.2 * cos(2*pi*8 *x)+0.5*x; 
%f = @(x) (x+1).*(x-0.5).*(x-0.25).*(x+3).*(x-0.6); %+ cos(2*pi*x).*(x <= 0.5); 

samples = walsh_sampling(f,N);

[idx, scales] = sph1_rect2(N, M, nbr_samples, j0);
idx = idx';

G = @(x, mode) cww_handle_1d_cs(x, mode, idx, log2N, log2M, ...
                                dwt_kernel, idwt_kernel, phi_walsh_pieces)

y = samples(idx);
opts_spgl1 = spgSetParms('verbosity', 1);
wc = spg_bpdn(G, y, noise, opts_spgl1); 

sc = wl_idwt_impl_from_kernel(wc, idwt_kernel);

A = cww_get_scaling_matrix(R+q+q2, R, wname, bd_mode);
t = linspace(0,1,N2);

x = A*sc; 

figure()
eps = 1e-14;
t1 = linspace(0,1-eps,N2)';
plot(t1,f(t1));
hold('on');
%t = linspace(0,1-eps,N)';
plot(t,x);
xlabel('t');
title(sprintf('vm: %d, N: %d, M: %d', vm, N, M));

% Truncated Walsh series
raw_samples = zeros([N2,1]);
ub = round(subsampling_rate*N);
raw_samples(1:ub) = samples(1:ub);
%s = 0;
%for n = 1:N
%    s = s + samples(n)*wal(n-1, t1);
%end


hold('on');
plot(t1,fastwht(raw_samples)*N2);


%%legend({'f(t)', 'GS rec.'});
%%legend({'f(t)', 'GS rec.', 'f_trunk'});
%legend({'f(t)', 'f_trunk'});

%fprintf('||f-f_{N,M}||_{inf}: %g\n', norm(f(t)-x, inf));
%fprintf('||f-f_trunk||_{inf}: %g\n', norm(f(t)-s', inf));








%X = zeros([round(0.1*N), N]);
%X(:, idx) = 1;
%imagesc(X); colormap('gray');





