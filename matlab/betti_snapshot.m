% This is a MATLAB function that reads in files containing the
% coordinate data of an agent-based model for zebrafish in the form of a matrix 
% and plots barcodes and bettigraphs for betti 0 and 1 at a specified time.
%
% Inputs: path to model file, time of interest (0 to 46). If the input is empty we run simulation 808 at time 46 by default.
% 
% Outputs: 5 plots.
% Author: Daniel Tolosa, Purdue University, 2022. Based on work by Melissa R. McGuirl, Brown University, 2019.

%----------------------------------- MAKE IT SO USER INPUTS PATH AND TIME ---------------------------------


function [handle] = betti_snapshot(path_to_file, time_pt)

if ~exist('path_to_file','var')
    % if parameter does not exist, default it to something
    path_to_file = 'WT_default/Out_WT_default_808';
end

if ~exist('time_pt','var')
    % if parameter does not exist, default it to something
    time_pt = 46;
end

load(path_to_file);

% MAKE SURE ALL THE NEEDED MODULES/FUNCTIONS ARE IN MATLAB PATH
% Run it from the directory that contains rip_on_ML.py

% we loaded all the cell coordinates of all cell types at each time, from 1 to 46

% select only the Mel cells that live inside the boundary

cutoff = 0.1*boundaryY(time_pt) ;

cells_mel = cellsM(find(cellsM(1:numMel(time_pt), 2,time_pt) > cutoff &  cellsM(1:numMel(time_pt), 2,time_pt) < boundaryY(time_pt) - cutoff), :, time_pt);

% get periodic distance matrix D_mel

[D_mel, ~, ~] = getPeriodicDistMats(cells_mel, boundaryX(time_pt));

% use ripser(python), from inside Matlab on D_mel, calling the function rip_on_ML. It is a little messy because of the file formats accross languages.

D_mel_py = py.numpy.array(D_mel);
barcodes_temp = py.rip_on_ML.main(D_mel_py);
barcodes_temp_0 = barcodes_temp(1);
barcodes_temp_1 = barcodes_temp(2);
barcodes_betti0 = double(barcodes_temp_0{1});
barcodes_betti1 = double(barcodes_temp_1{1});

% Below we plot stuff. We do a scatterplot at the given time, barcodes of betti 0, betti 1, then bettigraphs of betti 0 and 1.

handle = figure;
hold on;

scatter(cells_mel(:,1),cells_mel(:,2))

plot_barcodes_Dan(barcodes_betti0, 0)

plot_barcodes_Dan(barcodes_betti1, 1)

plot_bettigraph_Dan(barcodes_betti0, 0)

plot_bettigraph_Dan(barcodes_betti1, 1)

hold off;

end


