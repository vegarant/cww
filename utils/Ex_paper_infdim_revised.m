% This scripts samples a Daubechies scaling function with walsh samples, 
% and preform reconstruction with direct inversion, in addition to a
% finite and infinite dimensional model with wavelets. 
%
clear all;
fsize = 18;
lwidth = 2;
figure_visible = 'off';
dwtmode('per', 'nodisp');


R = 6;
q = 1;

M = 2^R;
N = 2^(R+q);
is_per = 1;
vm = 4;
j0 = 4; %computeJ0(vm); 

fact = 8;
nbr_samples = 16
%[idx, scales] = sph1_rect1(N, M, nbr_samples, j0);
idx = 1:nbr_samples;%simple_patt(N,nbr_samples);
idx_full = simple_patt(N, N/2);

idx1 = 1:16; % Inf dim
idx2 = 1:16;
idx3 = 1:32;
idx4 = 1:32; % Dir inv

if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
[dwt_kernel, idwt_kernel] = cww_compute_wave_kernels(log2M, dims, wname, bd_mode, j0);

w_idx = vm+1
nres = R+q-j0;
phi = zeros([2^(j0+nres),1]);
phi(w_idx) = 2^(nres/2);
f_values_samp = wl_idwt_impl(phi, wname, 'm', nres, 'bd_mode', bd_mode);

samples = fastwht(f_values_samp)*sqrt(N);

G = @(x, mode) cww_handle_1d_cs(x, mode, idx1, log2N, log2M, ...
                                dwt_kernel, idwt_kernel, phi_walsh_pieces)

r2 = 5;
r3 = 5;
% Create inf-dim matrix
A1 = create_A1d_mat(vm, j0, R+q, R, idx, is_per);

% Create fin-dim matrix
U1 = create_U_matrix_discrete(2^r2, vm, r2-j0, is_per);
U2 = create_U_matrix_discrete(2^r3, vm, r3-j0, is_per);
A2 = U1(idx2,:);
A3 = U2(idx3,:);

Nlarge = fact*N;
a = round(log2(Nlarge));

x = zeros([Nlarge,1]);
x(5) = 2^((a-j0+1)/2);
2^(j0-2)+1
S = get_wavedec_s(a, a+1-j0);
wname = sprintf('db%d', vm);
f = waverec(x, S, wname)';

samples = fastwht( f );
samples_full = samples(1:N)';

y1 = samples_full(idx1);
y2 = samples_full(idx2)*sqrt(2^r2);
y3 = samples_full(idx3)*sqrt(2^r3);
y4 = samples_full(idx4); % dir inv


opts = spgSetParms('verbosity', 1);
z1 = spg_bpdn(A1, y1, 0.0001, opts);
z2 = spg_bpdn(A2, y2, 1e-8, opts);
z3 = spg_bpdn(A3, y3, 1e-8, opts);

tmp2 = A2*z2;
tmp3 = A3*z3;
n2 = norm(y2-tmp2);
n3 = norm(y3-tmp3);


% Infinite recovery
if (is_per)
    S = get_wavedec_s(R, R-j0);
    sc1 = waverec(z1, S, sprintf('db%d', vm));
    X1 = get_scaling_matrix_per(vm, 4*N, M);
else 
    sc1 = IWT_CDJV_noP(z1, j0, vm);
    X1  = get_scaling_matrix_bd(vm, R, R+q);
end

rec_inf = X1*sc1; 

% Finite recovery
if (is_per)
    S2 = get_wavedec_s(r2, r2-j0);
    S3 = get_wavedec_s(r3, r3-j0);
    rec_fin      = waverec(z2, S2, sprintf('db%d', vm));
    rec_fin_full = waverec(z3, S3, sprintf('db%d', vm));
else 
    rec_fin      = IWT_CDJV(z2,j0,vm)
    rec_fin_full = IWT_CDJV(z3,j0,vm)
end

% Direct inversion
z4 = zeros([Nlarge, 1]);
z4(idx4) = y4;
rec_dir_inv = fastwht(z4)*Nlarge;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   Start plotting                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps = 1e-8;
t1 = linspace(0,1-eps, Nlarge);

% Direct inversion
fig = figure('Visible', figure_visible);

plot(t1,f, 'linewidth', lwidth);

hold('on');
plot(t1, rec_dir_inv, 'linewidth', lwidth);
legend({'f(t)', 'Recovery'},'Location', 'southwest', 'FontSize', fsize);
set(gca,'fontsize',fsize)

saveas(fig, 'rec_dir_inv1.png')

% Finite dim reconstruction
fig = figure('Visible', figure_visible);
plot(t1,f, 'linewidth', lwidth);

hold('on');
plot(linspace(0, 1-eps, length(rec_fin)), rec_fin, 'linewidth', lwidth);
set(gca,'fontsize',fsize)

legend({'f(t)', 'Recovery'},'Location', 'southwest', 'FontSize', fsize);

saveas(fig, 'rec_fin1.png')

% Finite full dim reconstruction
fig = figure('Visible', figure_visible);
plot(t1,f, 'linewidth', lwidth);

hold('on');
plot(linspace(0, 1-eps, length(rec_fin_full)), rec_fin_full, 'linewidth', lwidth);
set(gca,'fontsize', fsize)

legend({'f(t)', 'Recovery'},'Location', 'southwest', 'FontSize', fsize);

saveas(fig, 'rec_fin_full1.png')

% Finite dim reconstruction
fig = figure('Visible', figure_visible);
plot(t1,f, 'linewidth', lwidth);

hold('on');
plot(linspace(0, 1-eps, length(rec_inf)), rec_inf, 'linewidth', lwidth);
set(gca,'fontsize',fsize)

legend({'f(t)', 'Recovery'},'Location', 'southwest', 'FontSize', fsize);

saveas(fig, 'rec_inf1.png')


