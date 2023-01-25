% This is a MATLAB function whichs reads in files containing the
% coordinate data of an agent-based model for zebrafish in the form of a matrix and returns two CROCKER plots 
% (as seen on Topaz et al) 
% for betti 0 and 1 as functions of epsilon and time.
%
% Inputs:number of simulation; 
% 
% Outputs: 2 CROCKER plots.
% Author: Daniel Tolosa, Purdue University, 2022. Based on work by Melissa R. McGuirl, Brown University, 2019.

%----------------------------------- CHOOSE WHETHER USER INPUTS SIM NUMBER OR NOT  ---------------------------------

% MAKE SURE ALL THE NEEDED MODULES/FUNCTIONS ARE IN MATLAB PATH

% Uncomment the line below to get simulation number as user input
function [handle] = plot_crocker(sim_number)
file = sprintf('WT_default/Out_WT_default_%d', sim_number);

% Edit the line below to choose the simulation number
%function [handle] = plot_crocker()
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
barcodes_betti0 = cell(46,1);
barcodes_betti1 = cell(46,1);
D_mel=cell(46,1);
for i = 1 : 46 
	[D_mel{i}, ~, ~] = getPeriodicDistMats(cells_mel{i}, boundaryX(i));
	D_mel_py = py.numpy.array(D_mel{i});
	barcodes_temp = py.rip_on_ML.main(D_mel_py);
	barcodes_temp_0 = barcodes_temp(1);
	barcodes_temp_1 = barcodes_temp(2);
	barcodes_betti0{i} = double(barcodes_temp_0{1});
	barcodes_betti1{i} = double(barcodes_temp_1{1});
end

num_increments = 500;
eps_min=0;
thresholds_b0=zeros(46,1);
thresholds_b1=zeros(46,1);
for i =1:46
    temp_b0 = barcodes_betti0{i};
    endpoints_betti0=temp_b0(:,2);
    temp_b1 = barcodes_betti1{i};
    endpoints_betti1=temp_b1(:,2);
    thresholds_b0(i) = max(endpoints_betti0(~isinf(endpoints_betti0)))*(7/6);
    thresholds_b1(i) = max(endpoints_betti1(~isinf(endpoints_betti1)))*(7/6);
end
Threshold_b0 = max(thresholds_b0);
Threshold_b1 = max(thresholds_b1);
increment= zeros(2);
increment(1) = Threshold_b0/num_increments;
increment(2) = Threshold_b1/num_increments;
eps_values_b0 = eps_min : increment(1) : Threshold_b0;
eps_values_b1 = eps_min : increment(2) : Threshold_b1;
time_values = 1:1:46 ;

% Now we create two functions betti0(epsilon,time) , betti1(epsilon,time)

betti0 = zeros(46, num_increments+1);
betti1 = zeros(46, num_increments+1);

% j is time

for j = 1:46
    for i =1:length(eps_values_b0)
        betti0(j,i) =nnz(barcodes_betti0{j}(:,1) <= eps_min + increment(1)*i & barcodes_betti0{j}(:,2) >= eps_min + increment(1)*i );
    end
    for i =1:length(eps_values_b1)
        betti1(j,i) =nnz(barcodes_betti1{j}(:,1) <= eps_min + increment(2)*i & barcodes_betti1{j}(:,2) >= eps_min + increment(2)*i );
    end
end

% surf alone is just the plot

% surf(time_values, eps_values , betti);

% surfc includes a contour plot below, but honestly it doesn't do much
% diference.

% ----------------- NEED TO REWRITE THIS PART BELOW, IT'S THE ACTUAL PLOT -----------------------


handle = figure;
hold on;


contour(eps_values_b0 ,time_values, betti0, [0,1,2,3,4]);
LineWidth = 3 ;
caption = 'Betti(time,eps) of dimension ';
title(sprintf('%s %d', caption, 0));
xlabel('epsilon');
ylabel('time');
view([90 -90]) %// instead of normal view, which is view([0 90])

% comment out the line below if the colomap CROCKER is not available on workspace
colormap(gca, CROCKER)


figure;

contour(eps_values_b1 ,time_values, betti1, [0,1,2,3,4]);
LineWidth = 3 ;
caption = 'Betti(time,eps) of dimension ';
title(sprintf('%s %d', caption, 1));
xlabel('epsilon');
ylabel('time');

% comment out the line below if the colomap CROCKER is not available on workspace
colormap(gca, CROCKER)

view([90 -90]) %// instead of normal view, which is view([0 90])

% Maybe a scatterplot here of 4 or 5 times would make sense to compare with the surf data, one the same figure
% scatter(cells_mel(:,1),cells_mel(:,2))

hold off;

end