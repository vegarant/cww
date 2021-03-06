clear all

load('cww_defaults.mat') % load font size, line width, etc.

R = 5;
q = 1;
q2 = 6;
vm = 4;
is_per = 0;
subsampling_rate = 0.5;
noise = 1e-8;
dims = 1;
dest = 'plots';
disp_plot = 'on';
do_save = 0;
plot_walsh = 0;
location = 'northwest';

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
nbr_samples = round(N*subsampling_rate) 
if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
[dwt_kernel, idwt_kernel] = cww_compute_wave_kernels(log2M, dims, wname, bd_mode, j0);

%f = @(x) 2*x.*(x <= 0.5) + (2 - 2*x).*(x > 0.5);
f = @(x) 50*(x-0.5).*(x-0.05).*(x-0.95) + 2*x;
%f = @(x) cos(2*pi*x)  + 0.2 * cos(2*pi*8 *x)+0.5*x; 
%f = @(x) (x+1).*(x-0.5).*(x-0.25).*(x+3).*(x-0.6); %+ cos(2*pi*x).*(x <= 0.5); 

samples = cww_walsh_sampling_1d(f,N);

[idx, scales] = sph1_rect2(N, M, nbr_samples, j0);
idx = idx';

G = @(x, mode) cww_handle_1d_cs(x, mode, idx, log2N, log2M, ...
                                dwt_kernel, idwt_kernel, phi_walsh_pieces);

y = samples(idx);
opts_spgl1 = spgSetParms('verbosity', 1);
wc = spg_bpdn(G, y, noise, opts_spgl1); 

residual = norm(y-G(wc,1));
fprintf('Computed solution with error: %g\n', residual )

sc = wl_idwt_impl_from_kernel(wc, idwt_kernel);

A = cww_get_scaling_matrix(R+q+q2, R, wname, bd_mode);
t = linspace(0,1,N2);

x = A*sc; 

fig = figure('visible', disp_plot);

eps = 1e-14;
t1 = linspace(0,1-eps,N2)';
f_data = f(t1);

plot(t1,abs(f_data-x), 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
hold('on');


if plot_walsh
    raw_samples = zeros([N2,1]);
    raw_samples(1:nbr_samples) = samples(1:nbr_samples);
    walsh_approx = fastwht(raw_samples)*N2;
    plot(t1, abs(f_data - walsh_approx), 'color', cww_dflt.red, 'linewidth', cww_dflt.line_width);
end

if plot_walsh
    legend({'Error CS', 'Error Walsh'}, 'location', location, 'fontsize', cww_dflt.font_size)
else
    legend({'Error CS'}, 'location', location, 'fontsize', cww_dflt.font_size)
end
set(gca, 'FontSize', cww_dflt.font_size);


if do_save
    fname = sprintf('CS_1d_N_%d_M_%d_%s_%s', N, M, wname, bd_mode);
    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
end








