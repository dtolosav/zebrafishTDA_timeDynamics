%{
%_________________________________________________________________________________________________________________________________
%_________________________________________________________________________________________________________________________________

,ggggggggggg,                                          ,gggg,                                                                   
dP"""88""""""Y8,            I8      I8                ,88"""Y8b,                                   ,dPYb,                        
Yb,  88      `8b            I8      I8               d8"     `Y8                                   IP'`Yb                        
 `"  88      ,8P         8888888888888888  gg       d8'   8b  d8                                   I8  8I                        
     88aaaad8P"             I8      I8     ""      ,8I    "Y88P'                                   I8  8bgg,                     
     88""""Y8ba   ,ggg,     I8      I8     gg      I8'            ,gggggg,    ,ggggg,      ,gggg,  I8 dP" "8   ,ggg,    ,gggggg, 
     88      `8b i8" "8i    I8      I8     88      d8             dP""""8I   dP"  "Y8ggg  dP"  "Yb I8d8bggP"  i8" "8i   dP""""8I 
     88      ,8P I8, ,8I   ,I8,    ,I8,    88      Y8,           ,8'    8I  i8'    ,8I   i8'       I8P' "Yb,  I8, ,8I  ,8'    8I 
     88_____,d8' `YbadP'  ,d88b,  ,d88b, _,88,_    `Yba,,_____, ,dP     Y8,,d8,   ,d8'  ,d8,_    _,d8    `Yb, `YbadP' ,dP     Y8,
    88888888P"  888P"Y88888P""Y8888P""Y888P""Y8      `"Y8888888 8P      `Y8P"Y8888P"    P""Y8888PP88P      Y8888P"Y8888P      `Y8
                                                                                                 
%}
%_________________________________________________________________________________________________________________________________
%_________________________________________________________________________________________________________________________________
%
% Betti_Crocker is a matlab function that produces a (better) version of a CROCKER plot as in [Topaz et Al].
%
function [handle] = betti_crocker(sim_number,cell_type , dimension, Save)
if ~exist('sim_number','var')
    % if parameter does not exist, default it to something
    sim_number = 1;
end
if ~exist('cell_type','var')
    cell_type = "XanC";
end
if ~exist('dimension','var')
	dimension = 1; 
end
%
% load barcodes for all times for the selected simulation and cell type.
barcodes=cell(46,1);
% SOME FILES ARE MISSING FOR THE EARLIEST TIMES, BECAUSE THERE ARE NO BARS ON THOSE TIMES
for i = 1:46
	%input_file = sprintf('WT_barcodes/sim%d/%s/PD_%ssim%dtime%d_dim1', sim_number,cell_type , cell_type,sim_number, i);
	try
		input_file = sprintf('barcodes/sim%d/%s/PD_%ssim%dtime%d_dim%d', sim_number,cell_type , cell_type,sim_number, i, dimension);
		[barcodes{i}, ~, ~]= importdata(input_file);
		if length(barcodes{i}) == 0
			barcodes{i} = [0 0];
		end
	catch
		try
			input_file = sprintf('barcodes/sim%d/%s/BC_%ssim%dtime%d_dim%d', sim_number,cell_type , cell_type,sim_number, i, dimension);
			[barcodes{i}, ~, ~]= importdata(input_file);
			if length(barcodes{i}) == 0
				barcodes{i} = [0 0];
			end
		catch
			input_file = sprintf('barcodes/sim%d/%s/PD_%ssim%dtime%d_dim%d', sim_number,cell_type , cell_type,sim_number, i, dimension);
			disp(['Could not find file: ',input_file]);
			barcodes{i}= [0 0];
		end
	end
end
% Now we create a function of two variables: betti(epsilon,time)
num_increments = 500;
eps_min=0;
thresholds=zeros(46,1);
for i =1:46
    temp = barcodes{i};
    endpoints = temp(:,2);
    thresholds(i) = max(endpoints(~isinf(endpoints)));
end
Threshold = max(thresholds);
increment = Threshold/num_increments;
eps_values = eps_min : increment : Threshold;
tix = 100: 100: 1200;
eps_values = unique(sort([eps_values tix]));
time_values = 1:46 ;
betti = zeros(46, num_increments+1);
color_mask = zeros(46, num_increments+1);
% j is time point
for j = 1:46
    for i = 1: length(eps_values)
        betti(j,i) =nnz(barcodes{j}(:,1) <= eps_min + increment*i & barcodes{j}(:,2) >= eps_min + increment*i );
        %{
if betti(j,i) = 0
        	color_mask(j,i) = "grey";
        elseif betti(j,i) = 1
        	color_mask(j,i) = "blue";
        elseif betti(j,i) = 2
        	color_mask(j,i) = "yellow";
        elseif betti(j,i) = 3
        	color_mask(j,i) = "green";
        elseif betti(j,i) = 4
        	color_mask(j,i) = "dark green";
        elseif betti(j,i) > 4
        	color_mask(j,i) = "red";
        end
%}

    end
end

betti_flipped = flipud(betti);
h = heatmap(eps_values ,time_values, betti_flipped, 'ColorLimits',[0 5],'GridVisible','off','CellLabelColor','none','Colormap',turbo);
caption = sprintf('Betti %d as function of eps and time', dimension);
title(caption);
evals= [];
for i = 1 : length(eps_values)
	if ismember(eps_values(i),tix) | i == 1
		evals=[evals string(eps_values(i))];
	else
		evals=[evals ""];
	end
end
tvals= [];
for i = 1 : length(time_values)
	if mod(i,5) == 0 | i == 1
		tvals=[tvals string(time_values(i))];
	else
		tvals=[tvals ""];
	end
end
h.XDisplayData = eps_values;
h.XDisplayLabels = evals;
h.YDisplayLabels = flip(tvals);
axs = struct(gca); %ignore warning that this should be avoided
cb = axs.Colorbar;
cb.Ticks = [0 1 2 3 4 5];
cb.TickLabels = {'0','1','2',"3","4"," > 4"};
cb.Label.String = 'Homology classes';
axs.XAxis.TickLabelRotation = 0;   % horizontal
xlabel('epsilon');
ylabel('time');
if exist("Save","var")
	if Save == 1
		outfile = sprintf("Figures/Crocker/WT/sim%d/%s/crocker_WT_sim%d_%s_dim%d",sim_number,cell_type ,sim_number, cell_type, dimension);
		saveas(h, outfile);
	elseif Save ~= 0 & isstring(Save) == 1
		outfile = sprintf(Save+"crocker_WT_sim%d_%s_dim%d",sim_number,cell_type ,sim_number, cell_type, dimension);
		saveas(h, outfile);
	end
end
end