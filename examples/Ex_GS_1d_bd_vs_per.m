% This script uses generalized sampling to reconstruction two functions
% using wavelets with different boundary handling. Both functions are smooth, 
% but one of them is non-periodic on the interval [0,1]. We consider a wavelet 
% reconstruction bases with periodic and boundary corrected wavelet function.  
% 
% To generate all figures, you need to change some of the parameters 

clear all

load('cww_defaults.mat') % load font size, line width, etc.

R = 4;
q = 1;
q2 = 5;
vm = 4;
func_name = 'f';
is_per = 1;             % Set this to 0, to use boundary corrected wavelets
dest = 'plots';
disp_plot = 'off';
do_save = 1;
plot_walsh = 0;
location = 'north';

if (exist(dest) ~= 7) 
    mkdir(dest);
end

M = 2^R;
N = 2^(R+q);
N2 = 2^(R+q+q2);
wname = sprintf('db%d', vm);
j0 = 3; %cww_compute_j0(vm);
log2N = R+q;
log2M = R;

if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

% Add the `+ x` part to make the function smooth but non-periodict. 
f = @(x) cos(2*pi*x);%  + x;

samples = cww_sample_walsh_1d(f,N);

phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);

G = @(x,mode) cww_handle_1d(x, mode, R+q, R, wname, bd_mode, j0, phi_walsh_pieces); 

sc = lsqr(G,samples, [], 600);

x = cww_map_wcoeff_to_func_vals_1d(sc, R+q+q2, wname, bd_mode);
t = linspace(0,1,N2);


eps = 1e-14;
t1 = linspace(0,1-eps,N2)';
ft1 = f(t1);
ymax = max(max(ft1(:)), max(x(:)));
ymin = min(min(ft1(:)), min(x(:)));

fig = figure('visible', disp_plot);
plot(t1,f(t1), 'color', 'k', 'linewidth', cww_dflt.line_width);
legend(sprintf('%s(t)', func_name), 'location', location, 'fontsize', cww_dflt.font_size)
set(gca, 'FontSize', cww_dflt.font_size);
r = 0.1;
axis([0,1, ymin - r*sign(ymin)*ymin, ymax + r*sign(ymax)*ymax])

if do_save
    fname = sprintf('bd_vs_per_GS_1d_%s', func_name);
    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
    
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
end

fig = figure('visible', disp_plot);
plot(t1,x, 'color', 'k', 'linewidth', cww_dflt.line_width);
legend(sprintf('GS - %s.', bd_mode), 'location', location, 'fontsize', cww_dflt.font_size)
set(gca, 'FontSize', cww_dflt.font_size);
r = 0.1;
axis([0,1, ymin - r*sign(ymin)*ymin, ymax + r*sign(ymax)*ymax])

if do_save
    fname = sprintf('bd_vs_per_GS_1d_%s_N_%d_M_%d_%s_%s', func_name, N, M, wname, bd_mode);
    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
    
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
end






