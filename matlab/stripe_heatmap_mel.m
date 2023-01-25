function [handle] = stripe_heatmap_mel()

ttimes = 1:46;
ttimes = transpose(ttimes);

% Use the function num_stripes_mel to get the stripe counts for each simulation, here we do simulations 1 through 10

all_together = [];
stripe_count_by_sim = cell(0);

for i = 1 : 3
	stripe_count_by_sim{i} = [num_stripes_mel(i) ttimes];
	all_together = cat (1 , all_together , stripe_count_by_sim{i} )
end

% Now let's make a table

stripes = all_together(:,1);
times = all_together(:,2);

T = table( times , stripes);

% Now a heatmap from the table

FIG = figure('visible', 'off');
heatmap(T , 'times' , 'stripes');
savefig(FIG,"Test_heatmap");