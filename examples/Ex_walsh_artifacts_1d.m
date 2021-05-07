clear all

load('cww_defaults.mat') % load font size, line width, etc.

R = 4;
q = 1;
q2 = 6;
dest = 'plots';
disp_plot = 'off';
do_save = 1;

if (exist(dest) ~= 7) 
    mkdir(dest);
end

N = 2^(R+q);
N2 = 2^(R+q+q2);
log2N = R+q;

%f = @(x) 50*(x-0.6).*(x-0.5).*(x-0.5).*(x-0.05).*(x-0.95) + x;
f = @(x) 2*x.*(x <= 0.5) + (2 - 2*x).*(x > 0.5);

samples = cww_sample_walsh_1d(f,N);

fig = figure('visible', disp_plot);

eps = 1e-14;
t1 = linspace(0,1-eps,N2)';
plot(t1,f(t1), 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
hold('on');

raw_samples = zeros([N2,1]);
raw_samples(1:N) = samples;
walsh_approx = fastwht(raw_samples)*N2;
plot(t1, walsh_approx, 'color', cww_dflt.orange, 'linewidth', cww_dflt.line_width);

legend({'f(t)', 'Walsh'}, 'location', 'northwest', 'fontsize', cww_dflt.font_size)
set(gca, 'FontSize', cww_dflt.font_size);

fname = sprintf('walsh_artifacts_N_%d', N);

if do_save
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
end

