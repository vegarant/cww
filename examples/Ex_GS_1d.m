% This script reconstructs two hat functions using generalised sampling from 
% Walsh samples

clear all

load('cww_defaults.mat') % load font size, line width, etc.

R = 5;
q = 1;
q2 = 5;
vm = 2;
is_per = 0;
dest = 'plots';
disp_plot = 'off';
do_save = 1;
plot_walsh = 1;
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

if is_per
    bd_mode = 'per';
else
    bd_mode = 'bd';
end

box1 = @(x) (x <= 0.5);
box2 = @(x) (x >=0.5);
f = @(x) 4*x.*(x <= 0.25).*box1(x) + (2 - 4*x).*(x > 0.25).*box1(x) + 4*(x-0.5).*(x <= 0.75).*box2(x) + (2 - 4*(x-0.5)).*(x > 0.75).*box2(x);

samples = cww_sample_walsh_1d(f,N);

phi_walsh_pieces = cww_get_phi_walsh_pieces(R+q, R, wname, bd_mode, j0);

G = @(x,mode) cww_handle_1d(x, mode, R+q, R, wname, bd_mode, j0, phi_walsh_pieces); 

A = cww_get_scaling_matrix(R+q+q2, R, wname, bd_mode);
t = linspace(0,1,N2);

sc = lsqr(G,samples, [], 600);

x = A*sc; 





eps = 1e-14;
t1 = linspace(0,1-eps,N2)';
ft1 = f(t1);
ymax = max(max(ft1(:), max(x(:))));
ymin = min(min(ft1(:), min(x(:))));

fig = figure('visible', disp_plot);
plot(t1,f(t1), 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
hold('on');

plot(t1,x, 'color', cww_dflt.orange, 'linewidth', cww_dflt.line_width);


legend({'f(t)', 'GS'}, 'location', location, 'fontsize', cww_dflt.font_size)

set(gca, 'FontSize', cww_dflt.font_size);
r = 0.1;
axis([0,1, ymin - r*sign(ymin)*ymin, ymax + r*sign(ymax)*ymax])

if do_save
    fname = sprintf('GS_1d_N_%d_M_%d_%s_%s', N, M, wname, bd_mode);
    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
    
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
end



if plot_walsh
    fig = figure('visible', disp_plot);
    
    raw_samples = zeros([N2,1]);
    raw_samples(1:N) = samples;
    walsh_approx = fastwht(raw_samples)*N2;
    
    plot(t1,f(t1), 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
    hold('on');

    plot(t1, walsh_approx, 'color', cww_dflt.orange, 'linewidth', cww_dflt.line_width);


    legend({'f(t)', 'Walsh'}, 'location', location, 'fontsize', cww_dflt.font_size)

    set(gca, 'FontSize', cww_dflt.font_size);
    axis([0,1, ymin - r*sign(ymin)*ymin, ymax + r*sign(ymax)*ymax])
    if do_save
        fname = sprintf('GS_TW_1d_N_%d', N);
        fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
        
        saveas(fig, fullfile(dest, fname), 'png');
        saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
    end
end



