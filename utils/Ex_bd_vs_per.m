clear all

load('cww_defaults.mat') % load font size, line width, etc.

R = 4;
q = 1;
q2 = 5;
vm = 4;

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
j0 = cww_compute_j0(vm);
log2N = R+q;
log2M = R;

bd_mode_per = 'per';
bd_mode_bd = 'bd';

%f = @(x) 2*x.*(x <= 0.5) + (2 - 2*x).*(x > 0.5);
%f = @(x) 50*(x-0.5).*(x-0.05).*(x-0.95) + 2*x;
f = @(x) cos(2*pi*x)  + x; %0.2 * cos(2*pi*6 *x);
%f = @(x) (x+1).*(x-0.5).*(x-0.25).*(x+3).*(x-0.6); %+ cos(2*pi*x).*(x <= 0.5); 

samples = cww_sample_walsh_1d(f,N);

phi_walsh_pieces_per = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode_per, j0);
phi_walsh_pieces_bd = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode_bd, j0);

G_per = @(x,mode) cww_handle_1d(x, mode, R+q, R, wname, bd_mode_per, j0, phi_walsh_pieces_per); 
G_bd = @(x,mode) cww_handle_1d(x, mode, R+q, R, wname, bd_mode_bd, j0, phi_walsh_pieces_bd); 

A_per = cww_get_scaling_matrix(R+q+q2, R, wname, bd_mode_per);
A_bd = cww_get_scaling_matrix(R+q+q2, R, wname, bd_mode_bd);

t = linspace(0,1,N2);

sc_per = lsqr(G_per, samples, [], 600);
sc_bd = lsqr(G_bd, samples, [], 600);

x_per = A_per*sc_per; 
x_bd = A_bd*sc_bd; 

% Periodic reconstruction
fig = figure('visible', disp_plot);

eps = 1e-14;
t1 = linspace(0,1-eps,N2)';
plot(t1,f(t1), 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
hold('on');

plot(t1,x_per, 'color', cww_dflt.orange, 'linewidth', cww_dflt.line_width);

if plot_walsh
    raw_samples = zeros([N2,1]);
    raw_samples(1:N) = samples;
    walsh_approx = fastwht(raw_samples)*N2;
    plot(t1, walsh_approx, 'color', cww_dflt.red, 'linewidth', cww_dflt.line_width);
end

if plot_walsh
    legend({'g(t)', 'GS - per.', 'Walsh'}, 'location', location, 'fontsize', cww_dflt.font_size)
else 
    legend({'g(t)', 'GS - per.'}, 'location', location, 'fontsize', cww_dflt.font_size)
end

set(gca, 'FontSize', cww_dflt.font_size);

if do_save
    fname = sprintf('GS_bd_vs_per_1d_N_%d_M_%d_%s_%s', N, M, wname, bd_mode_per);
    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
    
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
end

% boundary wavelet reconstruction
fig = figure('visible', disp_plot);

eps = 1e-14;
t1 = linspace(0,1-eps,N2)';
plot(t1,f(t1), 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
hold('on');

plot(t1,x_bd, 'color', cww_dflt.orange, 'linewidth', cww_dflt.line_width);

if plot_walsh
    raw_samples = zeros([N2,1]);
    raw_samples(1:N) = samples;
    walsh_approx = fastwht(raw_samples)*N2;
    plot(t1, walsh_approx, 'color', cww_dflt.red, 'linewidth', cww_dflt.line_width);
end

if plot_walsh
    legend({'g(t)', 'GS - bd.', 'Walsh'}, 'location', location, 'fontsize', cww_dflt.font_size)
else 
    legend({'g(t)', 'GS - bd.'}, 'location', location, 'fontsize', cww_dflt.font_size)
end

set(gca, 'FontSize', cww_dflt.font_size);

if do_save
    fname = sprintf('GS_bd_vs_per_1d_N_%d_M_%d_%s_%s', N, M, wname, bd_mode_bd);
    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
    
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
end






