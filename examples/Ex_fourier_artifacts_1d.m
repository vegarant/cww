% This script demonstrate how a truncated Fourier series causes artifacts around 
% a discontinuity
%
% Change the value of R, to increse the number of Fourier coefficents. 


clear all

load('cww_defaults.mat') % load font size, line width, etc.

R = 5; % 8
q2 = 6;
dest = 'plots';
disp_plot = 'off';
do_save = 1;

if (exist(dest) ~= 7) 
    mkdir(dest);
end

M = 2^R;
N2 = 2^(R+q2);

f = @(x) (x <= 0.5);

t = linspace(0,1, N2)';

F = f(t);
samples = fft(F);
Mh = M/2;

samples_zero = zeros([N2,1]);
samples_zero(1:Mh+1) = samples(1:Mh+1);
samples_zero(N2-Mh+1:N2) = samples(N2-Mh+1:N2);

fourier_approx = real(ifft(samples_zero));

ymax = max(fourier_approx(:));
ymin = min(fourier_approx(:));
mag = ymax-ymin;

fig = figure('visible', disp_plot);

eps = 1e-14;
t1 = linspace(0,1-eps,N2)';
plot(t1,f(t1), 'color', cww_dflt.blue, 'linewidth', cww_dflt.line_width);
hold('on');


plot(t1, fourier_approx, 'color', cww_dflt.orange, 'linewidth', cww_dflt.line_width);

legend({'f(t)', 'Fourier'}, 'location', 'northeast', 'fontsize', cww_dflt.font_size)
set(gca, 'FontSize', cww_dflt.font_size);
ylim([ymin-0.05*mag, ymax+0.05*mag])

fname = sprintf('fourier_artifacts_N_%d', M);

if do_save
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);
end

