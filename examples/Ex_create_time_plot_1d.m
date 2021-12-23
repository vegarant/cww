% This script creates a plot of the execution time used by the fast transform
% and dense matrix-vector multiplication. Run the script `Ex_compute_time_1d.m`
% prior to running this script, to create the data. 

load('cww_defaults.mat') % load font size, line width, etc.

dest = 'plots';
if (exist(dest) ~= 7) 
    mkdir(dest);
end

src = 'data';
disp_plot = 'off';
vm = 4; 
is_per = 1;
j0 = 3;
q = 1;
R_min = 4;
R_max = 16;
% Only include the p+1 last values
p = 7;

if is_per
    bd_mode = 'per';
else 
    bd_mode = 'bd';
end

gray_color = [0.79, 0.79, 0.79];

types = {'forward', 'adjoint'};
for i = 1:length(types)

    type = types{i};
    fname = sprintf('time_%s_R_%d_%d_q_%d_db%d_%s_j0_%d.mat', type, R_min, R_max, q, vm, bd_mode, j0);
    data = load(fullfile(src, fname));
    if strcmpi(type, 'forward') 
        time_G = data.time_forward_G;
        time_X = data.time_forward_X;
    else
        time_G = data.time_adjoint_G;
        time_X = data.time_adjoint_X;
    end

    R_range = R_min:R_max;
    n = length(R_range);
    
    
    time_G1 = time_G(end-p:end);
    time_X1 = time_X(end-p:end);
    
    ideal_G2 = ones(1, p+1)*time_G1(1).*2.^(0:p);
    ideal_X2 = ones(1, p+1)*time_X1(1).*4.^(0:p);
    R_range1 = R_range(end-p:end);
    
    
    R = R_range(end-p);
    fig = figure('visible', disp_plot);
    
    semilogy(0:p, ideal_G2, '-.', 'color', gray_color, 'linewidth', cww_dflt.line_width);
    hold('on');
    semilogy(0:p, ideal_X2, '--', 'color', gray_color, 'linewidth', cww_dflt.line_width);
    hold('on');
    semilogy(0:p, time_G1, 'x', 'color', 'k', 'linewidth', cww_dflt.line_width);
    hold('on');
    semilogy(0:p, time_X1, '*', 'color', 'k', 'linewidth', cww_dflt.line_width);
    
    legend({'2^{k}Constant','4^{k}Constant', 'Fast transform', 'Dense matrix'}, ...
            'location', 'northwest', ...
            'fontsize', cww_dflt.font_size-5);
    xlabel('k', 'fontsize', cww_dflt.font_size)
    %ylabel('time (seconds)', 'fontsize', cww_dflt.font_size)
    set(gca, 'FontSize', cww_dflt.font_size);
    
    fname = sprintf('time_1d_%s_R_%d_q_%d_db%d_%s', type, R, q, vm, bd_mode);
    fprintf('saving as: %s.%s\n', fullfile(dest, fname), cww_dflt.plot_format(1:3));
    saveas(fig, fullfile(dest, fname), 'png');
    saveas(fig, fullfile(dest, fname), cww_dflt.plot_format);

end



