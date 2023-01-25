% THIS CODE PRODUCES THE SURFER PLOT FROM THE DATA IN TXT FORMAT, PROVIDED ON DROPBOX, SO IT ONLY RUNS ON THAT DATASET. 
% THERE IS ANOTHER SCRIPT THAT RUNS ON ANY SIMULATION FROM THE FIGSHARE BATCH.

function [handle] = plot_surfer_from_txt_sim(dimension)

% Load barcodes from odd times from 1 to 51. 

for i=1:26 % Will need to change this when considering more or less time values
    filenames(i) = sprintf("mel_distances_1_time%d_dim%d",i*2-1,dimension);
end

% M{} contains the data of birth-death of homology classes, as we vary
% epsilon and go through the time frames.

for i = 1 : length(filenames)
    M{i} = load(filenames(i));
end

num_increments = 500;
eps_min=0;
thresholds=[];
for i =1:length(filenames)
    endpoints=M{i}(:,2);
    thresholds(i) = max(endpoints(~isinf(endpoints)))*(7/6);
end
Threshold = max(thresholds);
increment = Threshold/num_increments;
eps_values = eps_min : increment : Threshold;
time_values = 0:2:51 ;

% Now we create a function betti(epsilon,time)

betti = zeros(length(filenames), num_increments+1);

% j is time

for j = 1:length(filenames)
    for i =1:length(eps_values)
        betti(j,i) =nnz(M{j}(:,1) <= eps_min + increment*i & M{j}(:,2) >= eps_min + increment*i );
    end
end
handle = figure;
hold on;

% surf alone is just the plot

% surf(time_values, eps_values , betti);

% surfc includes a contour plot below, but honestly it doesn't do much
% diference.

surfc(eps_values ,time_values, betti);
caption = 'Betti(time,eps) of dimension ';
title(sprintf('%s %d', caption, dimension));
xlabel('epsilon');
ylabel('time');
    
hold off;

end    