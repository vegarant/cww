clear all
R = 6;
q = 2;
q2 = 2; % Degree of high resolution image
vm = 2;
is_per = 1;
dims = 2;
subsampling_rate = 0.75;
noise = 0.01;

% Sampling patterns
r_factor = 2;
p_norm = inf;
nbr_levels = 10;
a = 1;
r0 = 1;

M = 2^R;
N = 2^(R+q);
N2 = 2^(R+q+q2);
wname = sprintf('db%d', vm);
j0 = cww_compute_j0(vm);
log2N = R+q;
log2M = R;
nbr_samples = round(N*N*subsampling_rate); 

%f = @(x,y) 5*(x-1).*(x-0.5).*(x-0.25).*cos(4*pi*y) + 5;
f = @(x,y) 10*cos(1.5*pi*x).*sin(3*pi*y); % + 10*(x>= 0.5);


if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
[dwt_kernel, idwt_kernel] = cww_compute_wave_kernels(log2M, dims, wname, bd_mode, j0);
[idx, str_id] = sph2_gsquare(N, nbr_samples, a, r0, nbr_levels);
%[idx, str_id] = sph2_2level(N, nbr_samples, p_norm, r_factor)


Nf = (2^3)*N;
im = phantom(Nf);
samples = cww_fastwht_2d(im);
samples = samples(1:N,1:N);

phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
G = @(x,mode) cww_handle_2d_cs(x, mode, idx, R+q, R, dwt_kernel, idwt_kernel, phi_walsh_pieces); 

y = samples(idx);
opts_spgl1 = spgSetParms('verbosity', 1);
wc = spg_bpdn(G, y, noise, opts_spgl1); 

wc = reshape(wc, M, M);
sc = wl_idwt_impl_from_kernel(wc, idwt_kernel);

W = cww_map_2d_wcoeff_to_func_vals(sc, R+q+q2, wname, bd_mode);

F_approx = zeros(N2,N2);
F_approx(1:N, 1:N) = samples;
F_approx = cww_fastwht_2d(F_approx)*N2*N2;

figure();
subplot(131); imagesc(im); colormap('gray'); title('Original');
subplot(132); imagesc(W); colormap('gray'); title('Wavelet approximation');
subplot(133); imagesc(F_approx); colormap('gray'); title('Walsh approximation');





















