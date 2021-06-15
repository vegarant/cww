% This file set variuous parameters to get consistency between example files
% The parameters found in this file, will be considered as default parameters

clear('all'); close('all');

% Create destination for the data
if (exist('../var') ~= 7) 
    mkdir('../var');
end

% cslib defaults
cww_dflt.font_size = 18; 
cww_dflt.marker_size = 8;
cww_dflt.color =[0,    0.4470,    0.7410];% [0,0,0]; % Black
cww_dflt.marker = 'o';
cww_dflt.marker_edge_color = [0,    0.4470,    0.7410]; % Black
cww_dflt.marker_face_color = [0,    0.4470,    0.7410]; % Black
cww_dflt.cbar_line_width = 1.5;
cww_dflt.line_width = 2;
cww_dflt.math_font = 'Latin Modern Math';
cww_dflt.image_format = 'png';
cww_dflt.plot_format = 'epsc';
cww_dflt.blue    = [0,    0.4470,    0.7410];
cww_dflt.yellow  = [1,1,0]; 
cww_dflt.green   = [102, 255, 51]./255;
cww_dflt.brown   = [153, 102, 51]./255;
cww_dflt.orange  = [255,140,0]./255;
cww_dflt.red     = [1, 0, 0];
cww_dflt.magenta = [1, 0, 1];
cww_dflt.cyan    = [0, 1, 1]; 
cww_dflt.black   = [0, 0, 0];

cww_dflt.line_color = 'k';
cww_dflt.spgl1_verbosity = 1;
cww_dflt.verbose = 1;
cww_dflt.data_path = '/home/vegant/db-cs/cslib_data';

cww_dflt.cmap_matrix = jet(256);

save('../var/cww_defaults.mat');



