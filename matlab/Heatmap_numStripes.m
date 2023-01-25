% This is a Matlab script that reads the barcodes for 100/200/500 simulations across times 1 to 46 and produces a heatmap
% illustrating the distribution of the stripe count in the patterns on zebrafish skin. This script uses the function
% num_stripes_mel_from_txt that is located in the folder "Matlab".
% Input =  None
% Output = heatmap in .fig format
% 
% Author: Daniel Tolosa, Purdue University, 2022.

addpath("Matlab/");

%function [handle] = Heatmap_numStripes()

ttimes = 1:46;
ttimes = transpose(ttimes);

% Use the function num_stripes_mel_from_txt to get the stripe counts for each simulation, here we do simulations 1 through 10

all_together = [];
stripe_count_by_sim = cell(0);

for i = 1 : 100
	stripe_count_by_sim{i} = [num_stripes_mel_from_txt(i) ttimes];
	all_together = cat (1 , all_together , stripe_count_by_sim{i} )
end

% Now let's make a table

stripes = all_together(:,1);
times = all_together(:,2);

T = table( times , stripes);

% Now a heatmap from the table

FIG = figure('visible', 'off');
heatmap(T , 'times' , 'stripes');
savefig(FIG,"Heatmap_numStripes_WT_Mel_100sims");