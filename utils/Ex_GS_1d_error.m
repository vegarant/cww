clear all

load('cww_defaults.mat') % load font size, line width, etc.

R = 4;
q = 1;
q2 = 6;
vm = 4;
is_per = 0;
dest = 'plots';
disp_plot = 'on';
do_save = 0;
plot_walsh = 1;
location = 'northwest';


if (exist(dest) ~= 7) 
    mkdir(dest);
end

M = 2^R;
N = 2^(R+q);
N2 = 2^(R+q+q2);
wname = sprintf('db%d', vm);
j0 =  cww_compute_j0(vm);
log2N = R+q;
log2M = R;

if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

%f = @(x) 2*x.*(x <= 0.5) + (2 - 2*x).*(x > 0.5);
f = @(x) 50*(x-0.5).*(x-0.05).*(x-0.95) + 2*x;
%f = @(x) cos(2*pi*x)  + 0.2 * cos(2*pi*6 *x);
%f = @(x) (x+1).*(x-0.5).*(x-0.25).*(x+3).*(x-0.6); %+ cos(2*pi*x).*(x <= 0.5); 

samples = cww_sample_walsh_1d(f,N);

phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);

G = @(x,mode) cww_handle_1d(x, mode, R+q, R, wname, bd_mode, j0, phi_walsh_pieces); 

A = cww_get_scaling_matrix(R+q+q2, R, wname, bd_mode);
t = linspace(0,1,N2);

sc = lsqr(G,samples, [], 600);

x = A*sc; 


fig = figure('visible', disp_plot);

eps = 1e-14;
t1 = linspace(0,1-eps,N2)';
f_data = f(t1);
plot(t1,abs(f_data-x), 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
hold('on');

if plot_walsh
    raw_samples = zeros([N2,1]);
    raw_samples(1:N) = samples;
    walsh_approx = fastwht(raw_samples)*N2;
    plot(t1, abs(f_data - walsh_approx), 'color', cww_dflt.red, 'linewidth', cww_dflt.line_width);
end

if plot_walsh
    legend({'Error GS', 'Error Walsh'}, 'location', location, 'fontsize', cww_dflt.font_size)
else
    legend({'Error GS'}, 'location', location, 'fontsize', cww_dflt.font_size)
end
set(gca, 'FontSize', cww_dflt.font_size);

if do_save
    fname = sprintf('GS_1d_N_%d_M_%d_%s_%s_error', N, M, wname, bd_mode);
    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
    
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
end

