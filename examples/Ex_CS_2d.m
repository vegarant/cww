% This script compares compressive sensing and truncated walsh reconstruction in 
% two dimensions, using a finite number of samples.

clear all

load('cww_defaults.mat') % load font size, line width, etc.

R = 8;
q = 1;
q2 = 3; % Degree of high resolution image
vm = 4;
is_per = 0;
dims = 2;
disp_plot = 'on';
dest = 'plots';

% Sampling patterns
K = 2^(R+q-2);
nbr_samples = K*K
nbr_levels = 50;
a = 2;
r0 = 1;

M = 2^R
N = 2^(R+q)
N2 = 2^(R+q+q2);
wname = sprintf('db%d', vm);
j0 = cww_compute_j0(vm);
log2N = R+q;
log2M = R;

h = @(x,y) 1 - (x > 0.25).*(x < 0.28).*(y > 0.6).*(y<0.8);
g = @(x,y) 1 - (x > 0.31).*(x < 0.34).*(y > 0.6).*(y<0.7);
box = @(x,y) (x > 0.2).*(x < 0.5).*(y > 0.6).*(y < 0.8);
f = @(x,y) 0.5*(cos(1.5*pi*x).*sin(3*pi*y) + 1).*h(x,y).*g(x,y).*(1-box(x,y)) + 0.5*box(x,y).*(cos(4*pi*(x-0.2)).*sin(4*pi*(y-0.6))+1).*h(x,y).*g(x,y);


if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

[dwt_kernel, idwt_kernel] = cww_compute_wave_kernels(log2M, dims, wname, bd_mode, j0);
[idx, str_id] = cil_sph2_gcircle(N, nbr_samples, a, r0, nbr_levels);

%figure('visible', disp_plot);
Z_CS = zeros([N,N], 'uint8');
Z_CS(idx) = uint8(255);
Z_TW = zeros([N,N], 'uint8');
Z_TW(1:K, 1:K) = uint8(255);
%imagesc(Z);

samples = cww_sample_walsh_2d(f, N);

phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
A = @(x,mode) cww_handle_2d_cs(x, mode, idx, R+q, R, dwt_kernel, idwt_kernel, phi_walsh_pieces); 

y = samples(idx);
opts_spgl1 = spgSetParms('verbosity', 1, 'iterations', 2000);
wc = spg_bpdn(A, y, 0.001, opts_spgl1); 

wc = reshape(wc, M, M);
sc = wl_idwt_impl_from_kernel(wc, idwt_kernel);

W = cww_map_wcoeff_to_func_vals_2d(sc, R+q+q2, wname, bd_mode);

t2 = linspace(0,1,N2);
[X2,Y2] = meshgrid(t2,t2);
F2 = f(X2,Y2);

F_approx = zeros(N2,N2);
F_approx(1:K, 1:K) = samples(1:K,1:K);
F_approx = cww_fastwht_2d(F_approx)*N2*N2;

mi = min(F2(:));
ma = max(F2(:));

F_plot = (F2 - mi)/(ma-mi);
F_approx_plot = (F_approx - mi)/(ma-mi);
W_plot = (W - mi)/(ma-mi);

F_approx_plot(F_approx_plot > 1) = 1;
F_approx_plot(F_approx_plot < 0) = 0;

W_plot(W_plot > 1) = 1;
W_plot(W_plot < 0) = 0;



figure();
subplot(131); imagesc(F2); colormap('gray'); title('Original');
subplot(132); imagesc(W_plot); colormap('gray'); title('Wavelet approximation');
subplot(133); imagesc(F_approx_plot); colormap('gray'); title('Walsh approximation');



fname_core_f = sprintf('acomp_CS_2d_N_%d_true_func', N);
fname_core_CS_rec = sprintf('acomp_CS_2d_N_%d_M_%d_K_%d_cs_rec', N, M, K);
fname_core_TW_rec = sprintf('acomp_CS_2d_N_%d_M_%d_K_%d_tw_rec', N, M, K);
fname_CS_samp_patt = sprintf('acomp_CS_2d_CS_N_%d_K_%d_samp_patt', N,K);
fname_TW_samp_patt = sprintf('acomp_CS_2d_TW_N_%d_K_%d_samp_patt', N,K);

imwrite(im2uint8(F2), fullfile(dest, [fname_core_f, '.', cww_dflt.image_format]));
imwrite(im2uint8(W_plot), fullfile(dest, [fname_core_CS_rec, '.', cww_dflt.image_format]));
imwrite(im2uint8(F_approx_plot), fullfile(dest, [fname_core_TW_rec, '.', cww_dflt.image_format]));
imwrite(Z_CS, fullfile(dest, [fname_CS_samp_patt, '.', cww_dflt.image_format]));
imwrite(Z_TW, fullfile(dest, [fname_TW_samp_patt, '.', cww_dflt.image_format]));

bd = 20;
lb_x = round(0.55*N2);
ub_x = round(0.85*N2);
lb_y = round(0.15*N2);
ub_y = round(0.45*N2);

idx_y = lb_y:ub_y;
idx_x = lb_x:ub_x;


fname_core_f_red = sprintf('acomp_CS_2d_N_%d_true_func_red', N);
F2_red = im2uint8(F2);
F2_red(F2_red == 0) = uint8(1);
F2_red(lb_x-bd:ub_x+bd, lb_y-bd:lb_y) = uint8(0);
F2_red(lb_x-bd:ub_x+bd, ub_y:ub_y+bd) = uint8(0);
F2_red(lb_x-bd:lb_x, lb_y-bd:ub_y+bd) = uint8(0);
F2_red(ub_x:ub_x+bd, lb_y-bd:ub_y+bd) = uint8(0);
cmap = gray(256);
cmap(1,:) = [1,0,0];

imwrite(F2_red, cmap, fullfile(dest, [fname_core_f_red, '.', cww_dflt.image_format]));


F2_crop = F2(idx_x, idx_y);
W_plot_crop = W_plot(idx_x, idx_y);
F_approx_plot_crop = F_approx_plot(idx_x, idx_y);


imwrite(im2uint8(F2_crop), fullfile(dest, [fname_core_f, '_crop.', cww_dflt.image_format]));
imwrite(im2uint8(W_plot_crop), fullfile(dest, [fname_core_CS_rec, '_crop.', cww_dflt.image_format]));
imwrite(im2uint8(F_approx_plot_crop), fullfile(dest, [fname_core_TW_rec, '_crop.', cww_dflt.image_format]));























