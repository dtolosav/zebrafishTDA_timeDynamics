function[num_stripes] = num_stripes_clust(barcodes, times, plotOn , plotName)

% get rid of the loop with lastest death time, which is related to domain size

[~,max_loc] = max(barcodes(:,2)); 
barcodes(max_loc,:)=[];

% find the distance of each point to the y=x line, we will use this to cluster the barcode data

barcodes_shifted = sin(atan(barcodes(:,2)./barcodes(:,1))- pi/4).*sqrt(barcodes(:,1).^2 + barcodes(:,2).^2);

% now evaluate using the Davies-Bouldin criterion to find the optimal cluster number, 
% we consider up to 8 clusters

eva = evalclusters(barcodes_shifted,"kmeans","DaviesBouldin","KList",[1:8]);
clus = eva.OptimalY;

% set num_stripes to NumObservations, then we will substract most until we get the actual stripe count.

num_stripes = eva.NumObservations +1;

% This handles the case where there is only one cluster or only one bar

if isequaln(eva.OptimalK ,NaN)
	num_clusters = 1;
	if mean(barcodes_shifted) < 90
		num_stripes = 0;
	end
	centroid = mean(barcodes_shifted);
	fprintf("centroid of only cluster at time %d = %d \n" , times , centroid);
	if plotOn == 1
		f=figure;
		scatter(barcodes(:,1),barcodes(:,2)-barcodes(:,1));
		xlabel("birth");
		ylabel("persistence");
		axlim = max(barcodes,[],"all")+5;
		ylim([0 axlim]);
		xlim([0 axlim]);
		titulo = sprintf("time=%d , numStripes=%d, numCl=%d",times, num_stripes ,num_clusters);
		title(titulo);
		set(gca,'fontsize',24);
		hold on;
		
	end
else
	num_clusters = eva.OptimalK;
	

	% no we get rid of the large clusters, this part is not great in the sense that we are 
	% writing code assuming we know the expected number of stripes, or at least an upper bound or 
	% an expected behavior of the data
	c_sizes = zeros(num_clusters,1);
	
	if num_clusters == 2
		cluster1 = barcodes_shifted(find(clus == 1));
		cluster2 = barcodes_shifted(find(clus == 2));
		centroid1 = mean(cluster1);
		centroid2 = mean(cluster2);
		if max(centroid1,centroid2) < 150
			num_stripes = 1;
		elseif centroid1 > centroid2
			num_stripes=num_stripes - length(cluster2);
		else
			num_stripes = num_stripes - length(cluster1);
		end
	else
		for i = 1 : num_clusters

	%  We get rid of clusters whose centroid is too close to the diagonal
			this_cluster = barcodes_shifted(find(clus == i));
			centroid = mean(this_cluster);
			fprintf("centroid of cluster %d at time %d = %d \n" , i , times , centroid);
			if centroid < 200
				num_stripes = num_stripes - length(this_cluster);
			end
		end
		if ~exist('times','var')
			    times = -1;
		end
	end
	if plotOn == 1
		f=figure;
		gscatter(barcodes(:,1),barcodes(:,2)-barcodes(:,1),clus,[],[],30,"on","birth","persistence");
		axlim = max(barcodes,[],"all")+5;
		ylim([0 axlim]);
		xlim([0 axlim]);
		titulo = sprintf("time=%d , numStripes=%d, numCl=%d",times,num_stripes,num_clusters);
		title(titulo);
		set(gca,'fontsize',24);
		hold on;
		
	end
end
if exist("plotName","var")
	saveas(f, plotName);
	%clf(f);
end