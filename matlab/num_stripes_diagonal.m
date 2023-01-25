function[num_stripes] = num_stripes_diagonal(barcodes, nnd , times, plotOn , plotName)

% if nearest-neighbor-distances are not provided then default it to 200, which is the persistence threshold
% used in Melissa's code, i.e. the one used in the paper [Volkening et al.]
if ~exist('nnd','var')
		nnd = 200 * ones(46,1);
	end

% find the distance of each point to the y=x line, we will use this to cluster the barcode data

barcodes_shifted = sin(atan(barcodes(:,2)./barcodes(:,1))- pi/4).*sqrt(barcodes(:,1).^2 + barcodes(:,2).^2);

% now evaluate using the Davies-Bouldin criterion to find the optimal cluster number, 
% we consider up to 8 clusters

eva = evalclusters(barcodes_shifted,"kmeans","DaviesBouldin","KList",[1:8]);
clus = eva.OptimalY;

% set num_stripes to NumObservations, then we will substract most until we get the actual stripe count.

num_bars = eva.NumObservations;
num_stripes = num_bars;

% This handles the case where there is only one cluster or only one bar

if isequaln(eva.OptimalK ,NaN)
	num_clusters = 1;
	if mean(barcodes_shifted) < nnd*3
		num_stripes = 0;
	end
	centroid = mean(barcodes_shifted);
	fprintf("centroid of only cluster at time %d = %d \n" , times , centroid);
	if plotOn == 1
		f=figure;
		scatter(barcodes(:,1),barcodes(:,2));
		xlabel("birth");
		ylabel("death");
		axlim = max(barcodes,[],"all")+5;
		ylim([0 axlim]);
		xlim([0 axlim]);
		titulo = sprintf("time=%d , numStripes=%d, numCl=%d",times, num_stripes ,num_clusters);
		title(titulo);
		set(gca,'fontsize',24);
		hold on;
		x= linspace(0,axlim);
		y= linspace(0,axlim);
		plot(x,y,"Color", "#ff0000");
	end
else
	num_clusters = eva.OptimalK;
	

	% no we get rid of the large clusters, this part is not great in the sense that we are 
	% writing code assuming we know the expected number of stripes, or at least an upper bound or 
	% an expected behavior of the data
	c_sizes = zeros(num_clusters,1);
	
	for i = 1 : num_clusters

%  This is the old version, where we get rid of large clusters
%		c_sizes(i)= length(find(clus == i));		
%		if c_sizes(i) > 4
%			num_stripes = num_stripes - c_sizes(i);
%		end

%  This is the new version, where we get rid of clusters whose centroid is too close to the diagonal
		this_cluster = barcodes_shifted(find(clus == i));
		centroid = mean(this_cluster);
		fprintf("centroid of cluster %d at time %d = %d \n" , i , times , centroid);
		if centroid < nnd*3
			num_stripes = num_stripes - length(this_cluster);
		end
	end
	if ~exist('times','var')
		    times = -1;
	end
	if plotOn == 1
		f=figure;
		gscatter(barcodes(:,1),barcodes(:,2),clus,[],[],30,"on","birth","death");
		axlim = max(barcodes,[],"all")+5;
		ylim([0 axlim]);
		xlim([0 axlim]);
		titulo = sprintf("time=%d , numStripes=%d, numCl=%d",times,num_stripes,num_clusters);
		title(titulo);
		set(gca,'fontsize',24);
		hold on;
		x= linspace(0,axlim);
		y= linspace(0,axlim);
		plot(x,y,"Color", "#ff0000");
		
	end
end
if exist("plotName","var")
	saveas(f, plotName);
	%clf(f);
end