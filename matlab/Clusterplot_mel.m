% This is a MATLAB function whichs reads in files containing the
% coordinate data of an agent-based model for zebrafish in the form of a matrix and counts stripes 
% as a function of time.
%
% Inputs:number of simulation; 
% 
% Outputs: list of number of stripes for times 1 through 46.
% Author: Daniel Tolosa, Purdue University, 2022. Based on work by Melissa R. McGuirl, Brown University, 2019.

%----------------------------------- CHOOSE WHETHER USER INPUTS SIM NUMBER OR NOT  ---------------------------------

% MAKE SURE ALL THE NEEDED MODULES/FUNCTIONS ARE IN MATLAB PATH

% Uncomment the line below to get simulation number as user input
function [handle] = num_stripes_mel(sim_number)
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

% we loaded all the cell coordinates of all cell types at each time, from 1 to 46

% select only the Mel cells that live inside the boundary at each time from 1 to 46

cutoff = 0.1*boundaryY ;
cells_mel = [];

for i = 1 : 46
	cells_mel{i} = cellsM(find(cellsM(1:numMel(i), 2,i) > cutoff(i) &  cellsM(1:numMel(i), 2, i) < boundaryY(i) - cutoff(i)), :, i);
end

% Get periodic distance matrix D_mel . 
% use ripser(python), from inside Matlab on D_mel, calling the function rip_on_ML. 
% It is a little messy because of the file formats accross languages.
%barcodes_betti0 = cell(46,1);
barcodes_betti1 = cell(46,1);
D_mel=cell(46,1);
for i = 1 : 46 
	[D_mel{i}, ~, ~] = getPeriodicDistMats(cells_mel{i}, boundaryX(i));
	D_mel_py = py.numpy.array(D_mel{i});
	barcodes_temp = py.rip_on_ML.main(D_mel_py);
	%barcodes_temp_0 = barcodes_temp(1);
	barcodes_temp_1 = barcodes_temp(2);
	%barcodes_betti0{i} = double(barcodes_temp_0{1});
	barcodes_betti1{i} = double(barcodes_temp_1{1});
end

num_increments = 500;
eps_min=0;
%thresholds_b0=zeros(46,1);
thresholds_b1=zeros(46,1);
for i =1:46
    %temp_b0 = barcodes_betti0{i};
    %endpoints_betti0=temp_b0(:,2);
    temp_b1 = barcodes_betti1{i};
    endpoints_betti1=temp_b1(:,2);
    %thresholds_b0(i) = max(endpoints_betti0(~isinf(endpoints_betti0)))*(7/6);
    thresholds_b1(i) = max(endpoints_betti1(~isinf(endpoints_betti1)))*(7/6);
end
%Threshold_b0 = max(thresholds_b0);
Threshold_b1 = max(thresholds_b1);
increment= zeros(2);
%increment(1) = Threshold_b0/num_increments;
increment(2) = Threshold_b1/num_increments;
%eps_values_b0 = eps_min : increment(1) : Threshold_b0;
eps_values_b1 = eps_min : increment(2) : Threshold_b1;
time_values = 1:1:46 ;

num_stripes = zeros(46,1);
for j = 1 : 46
    num_stripes(j) = length(find(barcodes_betti1{j}(:,1) < 150 & barcodes_betti1{j}(:,2) > 210));
end

num_stripes

%___________________________ SCATTERPLOTS CLUSTERED BY STRIPES ____________

sample_times =1:1:46;
Clusters = [];
for i =1 : 46
    t(i)=figure;
    caption1= 'Pattern at time ';
    caption2= '. Number of stripes =';
    temp_cells=cells_mel{sample_times(i)};
    Clusters{i} = clusterdata(temp_cells, 'Maxclust' , max(num_stripes(i), 1));
    scatter(temp_cells(:,1),temp_cells(:,2) , 50, Clusters{i} , 'filled');

    title(sprintf("%s %d %s %d", caption1 ,sample_times(i), caption2, num_stripes(i)));
    axis equal;
    %figname=strcat("Counting Stripes\WT_808\Mel\90-210\Stripe_Count_WT_808_Mel_90-210_time",string(sample_times(i)));
    %savefig(t(i),figname);
end
%savefig(t,"Counting Stripes\WT_808\Mel\90-210\Allfigs");