clear all

load('cww_defaults.mat') % load font size, line width, etc.

R = 7;
q = 2;
q2 = 6;
vm = 2;
is_per = 0;
s_factor = 4;
subsampling_rate = 1/(2^s_factor);
noise = 1e-8;
dims = 1;
dest = 'plots';
disp_plot = 'on';
do_save = 0;
plot_walsh = 1;
location = 'north';


if (exist(dest) ~= 7) 
    mkdir(dest);
end

M = 2^R;
N = 2^(R+q);
N2 = 2^(R+q+q2);
wname = sprintf('db%d', vm);
j0 = 2; % cww_compute_j0(vm);
log2N = R+q;
log2M = R;
nbr_samples = round(N*subsampling_rate) 
R_GS = round(log2(nbr_samples))-q;
N_GS = 2^(R_GS + q);

if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

f = @(x) 2*x.*(x <= 0.5) + (2 - 2*x).*(x > 0.5);
%f = @(x) 50*(x-0.5).*(x-0.05).*(x-0.95) + 2*x;
%f = @(x) cos(2*pi*x)  + 0.2 * cos(2*pi*6 *x);
%f = @(x) (x+1).*(x-0.5).*(x-0.25).*(x+3).*(x-0.6); %+ cos(2*pi*x).*(x <= 0.5); 

samples = cww_walsh_sampling_1d(f,N);
samples_GS = cww_walsh_sampling_1d(f,N_GS);

[idx, scales] = sph1_rect2(N, M, nbr_samples, j0);
idx = idx';

phi_walsh_pieces_CS = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
phi_walsh_pieces_GS = cww_get_phi_walsh_pieces(R_GS+q, R_GS, wname, bd_mode, j0);
[dwt_kernel, idwt_kernel] = cww_compute_wave_kernels(log2M, dims, wname, bd_mode, j0);


A_CS = cww_get_scaling_matrix(R+q+q2, R, wname, bd_mode);
A_GS = cww_get_scaling_matrix(R+q+q2, R_GS, wname, bd_mode);
t = linspace(0,1,N2);


G_GS = @(x,mode) cww_handle_1d(x, mode, R_GS+q, R_GS, wname, bd_mode, j0, phi_walsh_pieces_GS); 
G_CS = @(x, mode) cww_handle_1d_cs(x, mode, idx, log2N, log2M, ...
                                dwt_kernel, idwt_kernel, phi_walsh_pieces_CS);

y = samples(idx);
opts_spgl1 = spgSetParms('verbosity', 1);
wc = spg_bpdn(G_CS, y, noise, opts_spgl1); 

residual = norm(y-G_CS(wc,1));
fprintf('Computed CS solution with error: %g\n', residual )

sc_CS = wl_idwt_impl_from_kernel(wc, idwt_kernel);


sc_GS = lsqr(G_GS,samples_GS, [], 600);

x_GS = A_GS*sc_GS; 
x_CS = A_CS*sc_CS; 


fig = figure('visible', disp_plot);

eps = 1e-14;
t1 = linspace(0,1-eps,N2)';
f_data = f(t1);

if plot_walsh
    raw_samples = zeros([N2,1]);
    raw_samples(1:N) = samples;
    walsh_approx = fastwht(raw_samples)*N2;
    plot(t1, abs(f_data - walsh_approx), 'color', cww_dflt.red, 'linewidth', cww_dflt.line_width);
    hold('on');
end

plot(t1,abs(f_data-x_GS), 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
hold('on');
plot(t1,abs(f_data-x_CS), 'color', cww_dflt.yellow, 'linewidth', cww_dflt.line_width);


if plot_walsh
    legend({'Error Walsh', 'Error GS', 'Error CS'}, 'location', location, 'fontsize', cww_dflt.font_size)
else
    legend({'Error GS', 'Error CS'}, 'location', location, 'fontsize', cww_dflt.font_size)
end
set(gca, 'FontSize', cww_dflt.font_size);
set(gca, 'YScale', 'log')
if do_save
    fname = sprintf('Ex_error_GS_vs_CS_1d_nbr_s_%d_s_fac_%d_%s_%s', nbr_samples, s_factor, wname, bd_mode);
    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
    
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
end


