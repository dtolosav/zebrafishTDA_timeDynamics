% This is a MATLAB function whichs reads in files containing the
% coordinate data of an agent-based model for zebrafish in the form of a matrix and counts stripes of Xantophores
% as a function of time.
%
% Inputs:number of simulation; 
% 
% Outputs: list of number of stripes for times 2 through 46.
% Author: Daniel Tolosa, Purdue University, 2022. Based on work by Melissa R. McGuirl, Brown University, 2019.

%----------------------------------- CHOOSE WHETHER USER INPUTS SIM NUMBER OR NOT  ---------------------------------

% MAKE SURE ALL THE NEEDED MODULES/FUNCTIONS ARE IN MATLAB PATH

% Uncomment the line below to get simulation number as user input
function [handle] = num_stripes_xan(sim_number)
if ~exist('sim_number','var')
    % if parameter does not exist, default it to something
    sim_number = 808;
end
file = sprintf('WT_default/Out_WT_default_%d', sim_number);

% Edit the line below to choose run the script specifying the path to file
%function [handle] = plot_surfer()
%file='WT_default/Out_WT_default_808';

% --------------------------------------------------------------------------------------------------

load(file);

% we loaded all the cell coordinates of all cell types at each time, from 2 to 46

% select only the Xan cells that live inside the boundary at each time from 2 to 46

cutoff = 0.1*boundaryY ;
cells_xan = [];

for i = 2 : 1 : 46
	cells_xan{i} = cellsXc(find(cellsXc(1:numXanc(i), 2,i) > cutoff(i) &  cellsXc(1:numXanc(i), 2, i) < boundaryY(i) - cutoff(i)), :, i);
end

% Get periodic distance matrix D_mel . 
% use ripser(python), from inside Matlab on D_mel, calling the function rip_on_ML. 
% It is a little messy because of the file formats accross languages.
%barcodes_betti0 = cell(46,1);
barcodes_betti1 = cell(46,1);
D_xan=cell(46,1);
for i = 2:1 : 46 
	[D_xan{i}, ~, ~] = getPeriodicDistMats(cells_xan{i}, boundaryX(i));
	D_xan_py = py.numpy.array(D_xan{i});
	barcodes_temp_xan = py.rip_on_ML.main(D_xan_py);
	%barcodes_temp_xan_0 = barcodes_temp_xan(1);
	barcodes_temp_xan_1 = barcodes_temp_xan(2);
	%barcodes_betti0{i} = double(barcodes_temp_0{1});
	barcodes_xan_betti1{i} = double(barcodes_temp_xan_1{1});
end

num_increments = 500;
eps_min=0;
%thresholds_b0=zeros(46,1);
%thresholds_b1=zeros(46,1);
%for i =2:1:46
    %temp_b0 = barcodes_betti0{i};
    %endpoints_betti0=temp_b0(:,2);
    %temp_b1 = barcodes_betti1{i};
    %endpoints_betti1=temp_b1(:,2);
    %thresholds_b0(i) = max(endpoints_betti0(~isinf(endpoints_betti0)))*(7/6);
    %thresholds_b1(i) = max(endpoints_betti1(~isinf(endpoints_betti1)))*(7/6);
%end
%Threshold_b0 = max(thresholds_b0);
%Threshold_b1 = max(thresholds_b1);
%increment= zeros(2);
%increment(1) = Threshold_b0/num_increments;
%increment(2) = Threshold_b1/num_increments;
%eps_values_b0 = eps_min : increment(1) : Threshold_b0;
%eps_values_b1 = eps_min : increment(2) : Threshold_b1;
%time_values = 1:1:46 ;

num_stripes = zeros(46,1);
for j = 2: 1 : 46
    num_stripes(j) = length(find(barcodes_xan_betti1{j}(:,1) < 90 & barcodes_xan_betti1{j}(:,2) > 210));
end

num_stripes

%___________________________ SCATTERPLOTS CLUSTERED BY STRIPES ____________

sample_times =2:1:45;
Clusters_Xan = [];
for i =2:1 : 45
    h(i)=figure;
    caption1= ' Xan Pattern at time ';
    caption2= '. Number of stripes =';
    temp_cells=cells_xan{sample_times(i)};
    Clusters_Xan{i} = clusterdata(temp_cells, 'Maxclust' , max(num_stripes(i), 1));
    scatter(temp_cells(:,1),temp_cells(:,2) , 50, Clusters_Xan{i} , 'filled');

    title(sprintf("%s %d %s %d", caption1 ,sample_times(i), caption2, num_stripes(i)));
    axis equal;
    %figname=strcat("Counting Stripes\WT_808\Mel\90-210\Stripe_Count_WT_808_Mel_90-210_time",string(sample_times(i)));
    %savefig(h(i),figname);
end
%savefig(h,"Counting Stripes\WT_808\Xan\90-210\Allfigs_colored");