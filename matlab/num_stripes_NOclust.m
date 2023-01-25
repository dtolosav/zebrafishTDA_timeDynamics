function[num_stripes] = num_stripes_NOclust(barcodes, times, plotOn , plotName)

% find the distance of each point to the y=x line, we will use this to cluster the barcode data

barcodes_shifted = sin(atan(barcodes(:,2)./barcodes(:,1))- pi/4).*sqrt(barcodes(:,1).^2 + barcodes(:,2).^2);
num_stripes = sum(barcodes_shifted > 200);
		if ~exist('times','var')
			    times = -1;
		end
	if plotOn == 1
		f=figure;
		scatter(barcodes(:,1),barcodes(:,2)-barcodes(:,1));
		axlim = max(barcodes,[],"all")+5;
		ylim([0 axlim]);
		xlim([0 axlim]);
		titulo = sprintf("time=%d , numStripes=%d",times,num_stripes);
		title(titulo);
		set(gca,'fontsize',24);
		hold on;
		
	end
if exist("plotName","var")
	saveas(f, plotName);
	%clf(f);
end