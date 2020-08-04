clear all

load('cww_defaults.mat') % load font size, line width, etc.

R = 4;
q = 1;
q2 = 3;
vm = 2;
is_per = 0;
subsampling_rate = 0.25;
noise = 0.01;
dims = 1;
dest = 'plots';
disp_plot = 'on';
do_save = 0;

if (exist(dest) ~= 7) 
    mkdir(dest);
end


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

f = @(x) 50*(x-0.5).*(x-0.05).*(x-0.95) + x;
f = @(x) 2*x.*(x <= 0.5) + (2 - 2*x).*(x > 0.5);
%f = @(x) cos(2*pi*x)  + 0.2 * cos(2*pi*8 *x)+0.5*x; 
%f = @(x) (x+1).*(x-0.5).*(x-0.25).*(x+3).*(x-0.6); %+ cos(2*pi*x).*(x <= 0.5); 
r = 5;

w_idx = vm+1
nres = R+q+r;
phi = zeros([2^(j0+nres),1]);
phi(w_idx) = 2^(nres/2);
f_values_samp = wl_idwt_impl(phi, wname, 'm', nres, 'bd_mode', 'per');

samples = fastwht(f_values_samp);
samples = samples(1:N);

nres = R+q+q2-j0;
phi = zeros([2^(j0+nres),1]);
phi(w_idx) = 2^(nres/2);
f_values_plot = wl_idwt_impl(phi, wname, 'm', nres, 'bd_mode', 'per');


[idx, scales] = sph1_rect2(N, M, nbr_samples, j0);
idx = idx';

G = @(x, mode) cww_handle_1d_cs(x, mode, idx, log2N, log2M, ...
                                dwt_kernel, idwt_kernel, phi_walsh_pieces);

y = samples(idx);
opts_spgl1 = spgSetParms('verbosity', 1);
wc = spg_bpdn(G, y, noise, opts_spgl1); 

sc = wl_idwt_impl_from_kernel(wc, idwt_kernel);

A = cww_get_scaling_matrix(R+q+q2, R, wname, bd_mode);
t = linspace(0,1,N2);

x = A*sc; 

eps = 1e-14;
t1 = linspace(0,1-eps,N2)';

fig = figure('visible', disp_plot);

plot(t1, f_values_plot, 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
hold('on');

plot(t1,x, 'color', cww_dflt.yellow, 'linewidth', cww_dflt.line_width);

%raw_samples = zeros([N2,1]);
%raw_samples(1:nbr_samples) = samples(1:nbr_samples);
%walsh_approx = fastwht(raw_samples)*N2;
%plot(t1, walsh_approx, 'color', cww_dflt.red, 'linewidth', cww_dflt.line_width);

legend({'$f(t)$', 'CS-infinite'}, 'location', 'northeast', 'fontsize', cww_dflt.font_size, ...
        'Interpreter','latex');
set(gca, 'FontSize', cww_dflt.font_size);

fname = sprintf('CS_1d_wavelet_N_%d_M_%d_%s_%s', N, M, wname, bd_mode);
if do_save
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
end

[dwt_kernel, idwt_kernel] = cww_compute_wave_kernels(log2N, dims, wname, bd_mode, j0);

G_fin = @(x, mode) cww_findim_handle_1d_cs(x, mode, idx, N, ...
                                       dwt_kernel, idwt_kernel);

wc = spg_bpdn(G_fin, y*sqrt(N), noise, opts_spgl1); 
sc = wl_idwt_impl_from_kernel(wc, idwt_kernel);

t = linspace(0,1, N);

fig = figure('visible', disp_plot);

plot(t1, f_values_plot, 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
hold('on');

t = linspace(0,1, N);
plot(t, sc, 'color', cww_dflt.magenta, 'linewidth', cww_dflt.line_width);
legend({'$f(t)$', 'CS-finite'}, 'location', 'northeast', 'fontsize', cww_dflt.font_size, ...
        'Interpreter','latex');
set(gca, 'FontSize', cww_dflt.font_size);

fname = sprintf('CS_1d_wavelet_fin_wcrime_N_%d_M_%d_%s_%s', N, M, wname, bd_mode);

if do_save
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
end




