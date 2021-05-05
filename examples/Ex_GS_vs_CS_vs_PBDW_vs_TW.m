clear all

load('cww_defaults.mat') % load font size, line width, etc.

R = 4;
R_cs = 7;
q = 1;
q_cs = 1;
q2 = 7;
vm = 4;
is_per = 0;             % Set this to 0, to use boundary corrected wavelets
dest = 'plots';
disp_plot = 'off';
do_save = 1;
plot_walsh = 0;
location = 'north';
dims = 1;

if (exist(dest) ~= 7) 
    mkdir(dest);
end

M = 2^R;
N = 2^(R+q);
N_cs = 2^(R_cs+q_cs);
N2 = 2^(R+q+q2);
wname = sprintf('db%d', vm);
j0 = 3; %cww_compute_j0(vm);
log2N = R+q;
log2N_cs = R_cs + q_cs;
log2N_pbdw = R_cs + q + q_cs;
log2M = R;
nbr_samples = N;

if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

% Add the `+ x` part to make the function smooth but non-periodict. 
f = @(x) cos(2*pi*x).*(x <= 0.5) + (x > 0.5).*sin(6*pi*x) + 0.5*x; 

samples = cww_sample_walsh_1d(f,N_cs);
alpha = 2
idx = cww_sph1_power_law(N_cs, nbr_samples, alpha, R-1);

phi_walsh_pieces_gs = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);
phi_walsh_pieces_cs = cww_get_phi_walsh_pieces(log2N_cs, R_cs, wname, bd_mode, j0);
phi_walsh_pieces_pbdw = cww_get_phi_walsh_pieces(log2N_pbdw, R, wname, bd_mode, j0);
[dwt_kernel, idwt_kernel] = cww_compute_wave_kernels(R_cs, dims, wname, bd_mode, j0);

G = @(x,mode) cww_handle_1d(x, mode, R+q, R, wname, bd_mode, j0, phi_walsh_pieces_gs); 
A = @(x, mode) cww_handle_1d_cs(x, mode, idx, log2N_cs, R_cs, ...
                                dwt_kernel, idwt_kernel, phi_walsh_pieces_cs);
G_PBDW = @(x,mode) cww_handle_1d(x, mode, log2N_pbdw, R, wname, bd_mode, j0, phi_walsh_pieces_pbdw); 

sc_gs = lsqr(G, samples(1:N), [], 600);
wal_coeff_pbdw = G_PBDW(sc_gs, 1);
wal_coeff_pbdw(1:N) = samples(1:N);

y = samples(idx);
opts_spgl1 = spgSetParms('verbosity', 1, 'iterations',10000);
wc = spg_bp(A, y, opts_spgl1); 
sc_cs = wl_idwt_impl_from_kernel(wc, idwt_kernel);


x_gs = cww_map_wcoeff_to_func_vals_1d(sc_gs, R+q+q2, wname, bd_mode);
x_cs = cww_map_wcoeff_to_func_vals_1d(sc_cs, R+q+q2, wname, bd_mode);
t = linspace(0,1,N2);


fig = figure('visible', disp_plot);

eps = 1e-14;
t1 = linspace(0,1-eps,N2)';
ft1 = f(t1);
ymax = max(max(ft1(:), max(x_gs(:))));
ymin = min(min(ft1(:), min(x_gs(:))));

% -------------- GS ---------------

%fig = figure('visible', disp_plot);
%
%plot(t1,f(t1), 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
%hold('on');
%
%plot(t1,x_gs, 'color', cww_dflt.orange, 'linewidth', cww_dflt.line_width);
%
%set(gca, 'FontSize', cww_dflt.font_size);
%r = 0.1;
%axis([0,1, ymin - r*sign(ymin)*ymin, ymax + r*sign(ymax)*ymax])
%legend({'f(t)', sprintf('GS - %s.', bd_mode)}, 'location', location, 'fontsize', cww_dflt.font_size)
%
%if do_save
%    fname = sprintf('compare_met_GS_1d_N_%d_M_%d_%s_%s', N, M, wname, bd_mode);
%    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
%    
%    saveas(fig, fullfile(dest, fname), 'png');
%    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
%    
%    X_pattern = cww_visualize_1d_samp_patt(N_cs, 1:N); 
%    fig = figure('visible', disp_plot); 
%    imagesc(X_pattern); colormap('gray'); 
%    set(gca, 'ytick', []); 
%    set(gca, 'FontSize',cww_dflt.font_size);
%    %set(gca,'LooseInset',get(gca,'TightInset'));
%    %daspect([1, 15, 1])
%    fig.PaperUnits = 'inches';
%    fig.PaperPosition = [0, 0, 8, 1.2];
%    fig.PaperSize = [8, 1.2];
%    saveas(fig, [fullfile(dest, fname), '_samp'], cww_dflt.plot_format);
%end
%
%% -------------- PBDW ---------------
%
%fig = figure('visible', disp_plot);
%
%plot(t1,f(t1), 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
%hold('on');
%
%raw_samples = zeros([N2,1]);
%raw_samples(1:2^log2N_pbdw) = wal_coeff_pbdw;
%walsh_approx_pbdw = fastwht(raw_samples)*N2;
%
%plot(t1, walsh_approx_pbdw, 'color', cww_dflt.orange, 'linewidth', cww_dflt.line_width);
%
%set(gca, 'FontSize', cww_dflt.font_size);
%r = 0.1;
%axis([0,1, ymin - r*sign(ymin)*ymin, ymax + r*sign(ymax)*ymax])
%legend({'f(t)', 'PBDW'}, 'location', location, 'fontsize', cww_dflt.font_size)
%
%if do_save
%    fname = sprintf('compare_met_PBDW_1d_N_%d_M_%d_%s_%s', N, M, wname, bd_mode);
%    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
%    
%    saveas(fig, fullfile(dest, fname), 'png');
%    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
%end


% -------------- CS ---------------

fig = figure('visible', disp_plot);

plot(t1,f(t1), 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
hold('on');

plot(t1,x_cs, 'color', cww_dflt.orange, 'linewidth', cww_dflt.line_width);

set(gca, 'FontSize', cww_dflt.font_size);
r = 0.1;
ymax = max(max(ft1(:), max(x_cs(:))));
ymin = min(min(ft1(:), min(x_cs(:))));
axis([0,1, ymin - r*sign(ymin)*ymin, ymax + r*sign(ymax)*ymax])
legend({'f(t)', sprintf('CS - %s.', bd_mode)}, 'location', location, 'fontsize', cww_dflt.font_size)

if do_save
    fname = sprintf('compare_met_CS_1d_N_%d_M_%d_%s_%s', N, M, wname, bd_mode);
    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
    
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);

    X_pattern = cww_visualize_1d_samp_patt(N_cs, idx); 
    fig = figure('visible', disp_plot); 
    imagesc(X_pattern); colormap('gray'); 
    set(gca, 'ytick', []); 
    set(gca, 'FontSize',cww_dflt.font_size);
    %set(gca,'LooseInset',get(gca,'TightInset'));
    %daspect([1, 15, 1])
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0, 0, 8, 1.2];
    fig.PaperSize = [8, 1.2];
    saveas(fig, [fullfile(dest, fname), '_samp'], cww_dflt.plot_format);
    
end


%% -------------- Trunkcated Walsh (TW) ---------------
%
%fig = figure('visible', disp_plot);
%
%plot(t1,f(t1), 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
%hold('on');
%
%raw_samples = zeros([N2,1]);
%raw_samples(1:N) = samples(1:N);
%walsh_approx = fastwht(raw_samples)*N2;
%plot(t1, walsh_approx, 'color', cww_dflt.orange, 'linewidth', cww_dflt.line_width);
%
%legend({'f(t)', 'TW'}, 'location', location, 'fontsize', cww_dflt.font_size)
%
%set(gca, 'FontSize', cww_dflt.font_size);
%r = 0.1;
%ymax = max(max(ft1(:), max(walsh_approx(:))));
%ymin = min(min(ft1(:), min(walsh_approx(:))));
%axis([0,1, ymin - r*sign(ymin)*ymin, ymax + r*sign(ymax)*ymax])
%
%if do_save
%    fname = sprintf('compare_met_TW_1d_N_%d_M_%d_%s_%s', N, M, wname, bd_mode);
%    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
%    
%    saveas(fig, fullfile(dest, fname), 'png');
%    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
%end










