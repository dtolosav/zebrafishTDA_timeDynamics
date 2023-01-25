% This is Matlab function for connecting and visualizing paths of birth-death pairs 
% through time in persistence diagrams. It is an implementation of the algorithm described in 
% [V.Ciocanel et al., 2020] on section 2.3.2.
%
% Input: ( simulation number, cell type Mel/XanC/XanSn/IriL/IriD  , dimension for betti 0/1, path to save plot , timeVpersistence switch)
% 		timeVpersistence switches between plotting birth-death as in [Cio.] or plotting timeVpersistence. 
%		If you want the latter then input 1 for timeVpersistece, the default is b-d.
%
% Output: 1 figure containing all paths. The 3 most persitent will be colored Green,Yellow and Red,
% 		  respectively. For dimension 1 we produce a birth-death diagram, for dimension 0 we produce a
%		  time-persitence diagram, since birth is always 0 for connected components.
% Parameters: There is one built-in parameter, namely dist_thresh, which determines the maximum distance 
% between barcodes for them to be matched. It is set to be 100 nm. This number is approximately twice the  
% average cell-to-cell distance.
%
% Written by Daniel Tolosa, Purdue University, 2022.
% __________________________________________________________________________________________
function [handle] = ciocanel(sim_number,cell_type , dimension, Save, timeVpersistence)
if ~exist('sim_number','var')
    % if parameter does not exist, default it to something
    sim_number = 1;
end
if ~exist('timeVpersistence','var')
    % if parameter does not exist, default it to something
    timeVpersistence = 0;
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
% order barcodes by persistence (largest to shortest) at every time step
%barcodes{1} = [0 0];
persistences = cell(46,1);
for i = 1:46
	persistences{i}=barcodes{i}(:,2)-barcodes{i}(:,1);
	[ persistences{i} , mask ] = sort(persistences{i},'descend');
	barcodes{i} = barcodes{i}(mask,:);
	barcodes{i}(barcodes{i} == Inf) = 1500+i;
end
% for each pair of succesive times we compute the pairwise distances between barcodes
label_counter = 1;
dist_thresh = 100 ; % this value should be informed by the problem.
%
% We add a coordinate for labels on the barcodes. These will determine the paths.
for i =1 :46
	barcodes{i} = [barcodes{i} zeros(length(barcodes{i}(:,1)) , 1 )];
end
%
% FORMING PATHS:
%
for i = 1 : 45
	% pdist2 is stupid so we need to make all matrices of the same size. We are adding empty bars, i.e. bars that look like (-200, -200)
	l_i =length(barcodes{i});
	l_ip1 =length(barcodes{i+1});
	if l_i > l_ip1
		barcodes{i+1} = [barcodes{i+1} ; (-200)*ones( l_i-l_ip1 , 3)];
	else 
		barcodes{i} = [barcodes{i} ; (-200)*ones( -l_i+l_ip1 , 3)];
	end
	[D,I] = pdist2(barcodes{i+1}(:,1:2),barcodes{i}(:,1:2), 'chebychev' , 'Smallest' , 2 );
% I is a row vector whose jth entry is the index of barcodes{i+1} whose L_infty distance is minimum to the jth barcode in barcodes{i}
% we barcodes by pairs only, giving priority to the most persistent.
	j=1;
	if D(1,j) ~= 0 & D(1,j) < dist_thresh
		if barcodes{i}(j,3) == 0
			barcodes{i}(j,3) = label_counter ;
			barcodes{i+1}(I(1,j),3) = label_counter ;
			label_counter = label_counter + 1 ;
		else 
			barcodes{i+1}(I(1,j),3) = barcodes{i}(j,3) ;
		end
	end
	if length(I(1,:)) > 1
		j=2;
		if D(1,j) ~= 0 & D(1,j) < dist_thresh & ~ismember(I(1,j) , I(1,1:j-1))
			if barcodes{i}(j,3) == 0
				barcodes{i}(j,3) = label_counter ;
				barcodes{i+1}(I(1,j),3) = label_counter ;
				label_counter = label_counter + 1 ;
			else 
				barcodes{i+1}(I(1,j),3) = barcodes{i}(j,3) ;
			end
		elseif D(1,j) ~= 0 & D(2,j) ~= 0 & D(2,j) < dist_thresh & ismember(I(1,j) , I(1,1:j-1)) & ~ismember(I(2,j),I(1,1:j-1))
			I(1,j) = I(2,j);
			if barcodes{i}(j,3) == 0
					barcodes{i}(j,3) = label_counter ;
					barcodes{i+1}(I(1,j),3) = label_counter ;
					label_counter = label_counter + 1 ;
			else 
				barcodes{i+1}(I(1,j),3) = barcodes{i}(j,3) ;
			end
		end
		for j = 3 : length(I(1,:))
			if D(1,j) ~= 0 & D(1,j) < dist_thresh & ~ismember(I(1,j) , I(1,1:j-1))
				if barcodes{i}(j,3) == 0
					barcodes{i}(j,3) = label_counter ;
					barcodes{i+1}(I(1,j),3) = label_counter ;
					label_counter = label_counter + 1 ;
				else 
					barcodes{i+1}(I(1,j),3) = barcodes{i}(j,3) ;
				end
			elseif D(1,j) ~= 0 & D(2,j) ~= 0 & D(2,j) < dist_thresh & ismember(I(1,j) , I(1,1:j-1)) & ~ismember(I(1,j),I(1,1:j-2)) & ~ismember(I(2,j),I(1,1:j-1))
				I(1,j) = I(2,j);
				if barcodes{i}(j,3) == 0
						barcodes{i}(j,3) = label_counter ;
						barcodes{i+1}(I(1,j),3) = label_counter ;
						label_counter = label_counter + 1 ;
				else 
					barcodes{i+1}(I(1,j),3) = barcodes{i}(j,3) ;
				end
			end
		end
	end
end
%
% Now we just need to plot this monster.
% We put all the barcodes for all times together in a matrix, the rows of this matrix will have
% entries: ( birth, death, path label, time point, persistence )
num_bars_all = 0;
num_bars = zeros(46,1);
for i = 1 :46
	barcodes{i} = barcodes{i}(~ismember(barcodes{i},[-200 -200 -200],'rows'),:); % get rid of the extra barcodes of the form [-200 -200 -200]
	num_bars(i)= length(barcodes{i}(:,1));
	num_bars_all = num_bars_all + num_bars(i);
	barcodes{i} = [barcodes{i} i*ones(num_bars(i),1) barcodes{i}(:,2)-barcodes{i}(:,1) ];
end
M= zeros(num_bars_all,5);
% populating the matrix M
bar_counter=0;
M( (1 : num_bars(1)) , : ) = barcodes{1};
for i = 1 : 45
	bar_counter = bar_counter + num_bars(i);
	M( (bar_counter+1 : bar_counter+num_bars(i+1) ) , : ) = barcodes{i+1};
end
% consider a matrix Mp consisting on the barcodes with label ~= 0
Mp= M(find(M(:,3) ~= 0),:);
[~ , ind ] = max(Mp(:,5));
mpp_label = Mp(ind,3); % mpp = most persistent path: label of the path that contains the most persitent bar across all time points.
mpb_time = Mp(ind, 4); % mpb = most persistent bar: time point at which the most persistent bar is found
%
% Now we find bars in the largest path, i.e. find all bars with label mpp_label, we store them in 
% "largest_path", and sort it by time point
largest_path = Mp(find ( Mp(:,3) == mpp_label   ) ,:);
[~ , mask2 ]=sort(largest_path(:,4));
largest_path = largest_path( mask2 ,:);
%
% now let's find the second largest path
% We compute a new matrix M2, which is M withouth the largest path
M2 = Mp(find(Mp(:,3) ~= mpp_label) , : );
[~ , ind ] = max(M2(:,5));
second_mpp_label = M2(ind,3);
second_mpb_time = M2(ind, 4);
% find bars in second largest path, i.e. find all bars with label second_mpp_label
second_largest_path = M2(find ( M2(:,3) == second_mpp_label   ) ,:);
% sort second_largest_path by time point
[~ , mask2 ]=sort(second_largest_path(:,4));
second_largest_path = second_largest_path( mask2 ,:);
% now let's find the THIRD largest path
%
M3 = M2(find(M2(:,3) ~= second_mpp_label) , : );
[~ , ind ] = max(M3(:,5));
third_mpp_label = M3(ind,3);
third_mpb_time = M3(ind, 4);
% find bars in THIRD largest path, i.e. find all bars with label third_mpp_label, and sort them by time point
third_largest_path = M3(find ( M3(:,3) == third_mpp_label   ) ,:);
[~ , mask2 ]=sort(third_largest_path(:,4));
third_largest_path = third_largest_path( mask2 ,:);
%
%-------------------THERE ARE SOME DEPRECATED PLOTS BELOW, THAT MAY COME IN HANDY LATER, SO I WON´T DELETE THE CODE YET -------------
% now we are ready to plot the largest paths on a b-d figure
%{
Paths_bd=figure;
hold on;
plot(largest_path(:,1),largest_path(:,2), "LineWidth", 3);
plot(second_largest_path(:,1),second_largest_path(:,2), "LineWidth", 3);
plot(third_largest_path(:,1),third_largest_path(:,2),"color",[1 0 1], "LineWidth", 3);
title("Most significant topological features in time",'FontSize', 24);
xlim([0 1200]);
ylim([0 1200]);
xlabel("birth",'FontSize', 20);
ylabel("death",'FontSize', 20);
axis square;
%}
% now plot the largest paths on a b-persistence figure
%{
Paths_bp=figure;
hold on;
plot(largest_path(:,1),largest_path(:,5), "LineWidth", 3);
plot(second_largest_path(:,1),second_largest_path(:,5), "LineWidth", 3);
plot(third_largest_path(:,1),third_largest_path(:,5),"color",[1 0 1], "LineWidth", 3);
title("Most significant topological features in time",'FontSize', 24);
xlim([0 1200]);
ylim([0 1200]);
xlabel("birth",'FontSize', 20);
ylabel("persistence",'FontSize', 20);
axis square;
%}
% now we plot all the bars of all times with persistence > 60
%M_many = M(find(M(:,5) > 60) , : );
%BIGplot = figure;
%scatter(M_many(:,1),M_many(:,2),4,'filled');
%{
scatter(M(:,1),M(:,2),6,'filled');
title("All bars, all times","FontSize",24);
xlim([0 1200]);
ylim([0 1200]);
xlabel("birth",'FontSize', 20);
ylabel("death",'FontSize', 20);
axis square;
%}
%-----------------------END OF DEPRECATED CODE-------------------------------------------------------
%
% Plot all paths and color the 1st, 2nd and 3rd most persistent paths green, yellow and red, respectively.
if dimension == 0 | timeVpersistence == 1
	AllPathsPlot_bd = figure;
	hold on;
	num_labels = max(Mp( : , 3 ));
	for i = 1 :num_labels
		ipath = Mp(find(Mp( : , 3) == i),:);
		if i == mpp_label
			first_plt = plot(largest_path(:,4),largest_path(:,5),"color","#008040"  ,  "LineWidth", 3, "DisplayName", "Path that achieves highest persistence");
		elseif i == second_mpp_label
			second_plt = plot(second_largest_path(:,4),second_largest_path(:,5), "color","#ffcd58","LineWidth", 3,"DisplayName", "Path with 2nd highest persistence");
		elseif i == third_mpp_label
			third_plt = plot(third_largest_path(:,4),third_largest_path(:,5),"color","#80000b", "LineWidth", 3, "DisplayName", "Path with 3rd highest persistence");
		else
			plot(ipath(:,4),ipath(:,5), "LineWidth", 3);
		end
	end
	M_non_paths = M(find( M(:,3) == 0 & M(:,5) > 50),:);
	non_paths = scatter( M_non_paths(:,4), M_non_paths(:,5),30, "blue" , "DisplayName", "unmatched bars with persistence > 50");
	legend([first_plt second_plt third_plt non_paths]);
	titulo= sprintf("Paths of topological features in time (Betti %d)", dimension);
	title(titulo,'FontSize', 24);
	xlim([0 47]);
	ylim([0 1500]);
	xlabel("time in days",'FontSize', 20);
	ylabel("persistence",'FontSize', 20);
	axis square;
else
	AllPathsPlot_bd = figure;
	hold on;
	num_labels = max(Mp( : , 3 ));
	for i = 1 :num_labels
		ipath = Mp(find(Mp( : , 3) == i),:);
		if i == mpp_label
			first_plt = plot(largest_path(:,1),largest_path(:,2),"color","#008040"  ,  "LineWidth", 3, "DisplayName", "Path that achieves highest persistence");
		elseif i == second_mpp_label
			second_plt = plot(second_largest_path(:,1),second_largest_path(:,2), "color","#ffcd58","LineWidth", 3,"DisplayName", "Path with 2nd highest persistence");
		elseif i == third_mpp_label
			third_plt = plot(third_largest_path(:,1),third_largest_path(:,2),"color","#80000b", "LineWidth", 3, "DisplayName", "Path with 3rd highest persistence");
		else
			plot(ipath(:,1),ipath(:,2), "LineWidth", 3);
		end
	end
	M_non_paths = M(find( M(:,3) == 0 & M(:,5) > 50),:);
	non_paths = scatter( M_non_paths(:,1), M_non_paths(:,2),30, "blue" , "DisplayName", "unmatched bars with persistence > 50");
	legend([first_plt second_plt third_plt non_paths]);
	titulo= sprintf("Paths of topological features in time (Betti %d)", dimension);
	title(titulo,'FontSize', 24);
	xlim([0 1200]);
	ylim([0 1200]);
	xlabel("birth",'FontSize', 20);
	ylabel("death",'FontSize', 20);
	axis square;
end
if exist("Save","var")
	if Save == 1
		outfile = sprintf("Figures/Paths/WT/sim%d/%s/paths_WT_sim%d_%s_dim%d",sim_number,cell_type ,sim_number, cell_type, dimension);
		saveas(AllPathsPlot_bd, outfile);
	elseif Save ~= 0 & isstring(Save) == 1
		outfile = sprintf(Save+"paths_WT_sim%d_%s_dim%d",sim_number,cell_type ,sim_number, cell_type, dimension);
		saveas(AllPathsPlot_bd, outfile);
	end
end
end
   % (&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%&&&&&&%  
   % (&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%&&&&&&&&%%%%%%%%%%%%%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%&%  
   % (&%%%%%%%&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%&%%%%%%%%%%#&%  
   % (&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&%%%%%%%%%%%%%%%%%%%%%%%&&%%%%%%%%%%%%%%%%%&&%%%%%%%%&%  
   % (&&%%%%%%%%%%%%%%%%%%%%%%&&%%%%%%%%%%%%%%%%%%%&&&%%%%%%%%%%%%&&&%%%%%%%%%%%%%%%%%%%&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&&&&%%%%%%%&&&%%%%%%%%&&%%%%%%%%%%%%%%%&&&&&&%  
   % (&%%%%%%%%%%%%%%%%%%%%%%%%%&&%%%%%%%%%%%%%%%%%%%&&&%%%%%%%%%%%&&(#&&&&&&&&&&&&&%#(((((##&&&%%%%%%%%%%%%%%%%%%%%%&&&%%%%%%%%%%%%%%%%%%%%%&&&%%%%%%%%%%&&%%%%%%%%%%%%%%%%%%%%%#&%  
   % (&%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&%%%%%%%%%%%%%%%%%%%%%%%%%%&&&#((((((((((((((((((((((((((((((%&&&&&&&&&&&&&&&&#(////(&&&&&&&%%%%%&&&&&&(/&&%%%%%%%%%%%%%%%%%%%%&&&&&&%%%%%%%%%&%  
   % (&%%%%%%%%%%%%%%&&&&%%%%%%%%%&&&%&&&&%%%%%%%%%%%%%&&&&&&%(((((((((((((((((((((((((((((((((((((((((((((((((#&/////////////(((((((((((////&&#&&&%%%%%%%%%%%%%&&%((((((((#%&&&&%&%  
   % (&%%%%%%%%%%%%&&(((((#%%&&%#((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((%%///%%(((((((((##%%%%%%%&(//(&&%((((#%&&&&&&&&#(((((((((((((((((//&%  
   % (&#&&&&&&&&&%(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((&%(//%#((((((((((((((((((%#//%&&#((((((((((((((((((((((((((((((((//&%  
   % (&///((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((#####%%%((/((((((((((((((((((((((((((((((&%(//&(((((((&#((((&(((((%(((&&&#((((((((((((((((((((((((((((((((//&%  
   % (&/////(((((((((((((((((((((&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#(((((((((((((((((((((((((((((((&#///&((((((%%((((&((((((&(//&&&(((((((((((((((((((((((((((((((((//&%  
   % (&/////(((((((((((((((((((((&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#(((((((((((((((((((((((((((((((&#(//&(((((#&(((#&(((((((&/((&&&(((//((((((((((((((((((((((((((((//&%  
   % (&//////////////////////((((&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#/((///////((((((((((((((((//(((&#(//&(((((&(((#&((((((((&//(&&&((/////((((((((((((((((((((((((((//&%  
   % (&//////////(/&&(///////((((%&%%%&&&&&&&&%%%%&&%%&&&&&&&&&&&&&&&&&&&&&&&&&#/(//////////////////////////(((&(///&((((#(((#&((((((((#&/(%&&&((///////((((((((((((((((((((((((//&%  
   % (&//////////(%&(&((/////((((%&%%%%%%%%%%%%%%%&&&%%%%%%%%%%%%%%%%%%%%&&&&&&#/(///////////////////////////(/&(///&(((((((#&(((&(%#((%#/(%&&%(///////////////////(((////////////&%  
   % (&/////////(%&,*/&((////((((&&%%%%%%%%%%%%&#*(&%%%%%%%%%%%%%%&&&#((/%&&&&&#//////////////////////////////(&(///&((((((#&(((&(%#(((&(//&&&#(//////////////////((&%(///////////&%  
   % (&*////////&&,,,,#&///////((&&%%%%%%%%&%,,*,/&%%%%%%%%%%%%%%&%/////**///&&%(#&%//////////////////////////(&(///&(((((/&(((&((((&/(&//(&&&((////////////////////&%&%(/////////&%  
   % (&*///////&&,,,,,,&&//////((&&%%%%%%&*,,,,,,&%%%%%%%%%&%%%%%&////////*#&&&&&///%&/////////////////////////&(///&(((((((((%((((((((&/(#&&&(/////////////////////(&*#&#(///////&%  
   % (&*//////&#,,...,,*&(///////&&%%%#&/*,,,,,,%%%####%%&**&%%%%&///(&&&&&&&&&(/////&#////////////////////////&#///&%#####%%%%%%###(((&/(%&&&(//////////////////////&#*,&&//////*&%  
   % (&*////(&*,......,,&%///////&&###&.,,/&#,,,&#######%%,,,&%##%&&&&&&&&&&(//**///&&&&&#/////////////////////%%/////////////////////////&&&%(//////////////////////%&,,,&&/////*&%  
   % (&*///%&.,...,*#,,,&%///////&%#%(,,,&*%/,,*&#######&,*,,,&&&//////*////*****/%&#(////(&&(/////////////////#&(((((((((###%%%&&&&&&&&&&&&&#///////////////////////(&,,,,&%////*&%  
   % (&*//#&,....,/%&,,,&%///////&%&*,,.@,*%/,*,&######&,,%&,,#&///*********///(&&(///////////&%/////////////////#&&&&&&&&&&&&&&&&&&&&&&&&&%(////////////////////////(&,,,,%&////*&%  
   % (&*/*&*,...,*%/&,,,&#///////&&*,,,@,*,##,,,&#####&*,,&,&,,(&&(/////*/#&&&&////////////#&&(*(&&//////////////////////////////////////////////////////////////////&&,,,,/&///(&&%  
   % (&#/(&.,...*&,/&,,*&(/****/*&#,,,@,,,,,&.,,,&(((&*,,,%,,&,,,&((&&#//////////////////&(.. &&&&%,&//****//////////////////////////////////////////////////////////&*,.,,*&//&(.&%  
   % (&%&/&,,...&**/&,,#&//****//&(,,%,,,,,,#(,,,,&(&*,,.@,,*(*,,#(&%&/////#&&##&&&/////&. ./&&&&&&&.@//***//////////////////////////////////(&&///////////////*////&/,,**,*&&&...&%  
   % (&,/&&(,,.,@,*(%,,&&/*&#/*//&(,((,,,,,,*%/,,,,*,,,,&*,,,./,(&&&&&///#%........*&(/(&   &&&&&&&&(&(/***//////////////////////////**////&&&&////////////////*//(&*,,#%/,.&#,,..&%  
   % (&..,&&,,,,%,*##,*&///%&&*//&/,@,,,,,,,,,(#,,.,...&*,,,,,&&&&&&&&//##..&&&&&&*..&(/@  .#&&&&&&&%&//**//*/#&&&&&#*/******************&&,/&///****************&&.,,&/%(,.......&%  
   % (&...,*&*,,%,,&/,%&**/%#(%(/&/,&,,,,,,,,,,,@.,,,.&*,,*(&&&&&&&&&&//@ .#&&&&&&&,.#&//@ ../&&&&&&#(////&&&&&&&&&&&&&***************/&&,,,&&//**************//&/,,.&,*(%,......&&%  
   % (&,......,,%,,&*,&%***##,,&#&(*,,,,,,,,,,,,,/&.,@,**(&&&&&&&&&&&&//& .&&&&&&&&,.&%////%&&//////////#&&&&&&&&&&&&&&/***********//&&,,,,,%&//*************/#&.,,.&,,,,&....,(&.&%  
   % (&(%......,#,,@,*&(***%/,,,&&(*,,,,,,,,,,,,,,**&**,&&&&&&&&&&&&&&//%*..(&&&&&*,(&///(&/////////////%&&&&&&&&&&&&&#/***********&@,,#(&,,.&#/************/%&,,,.@,,,,,&,,,.&/,,&%  
   % (&.,&.,..,(,,*&,#&****&/..,*&(*,,,,,,,,,,,,,,,,,,*%&&&&&&&&&&&&&%///#&.....,./&#//////////&//////////#&&&&&&&&&#/***********&&,,/&,*/#,,,*&//**********(%,..,%*,,,,*(%,.&*,,.&%  
   % (&.,,#%..,/,,#%,&&****@,...,%@,,,,,,,,,,,,,,,,,,,,&&&&&&&&&&&&&&///////(%%%#////////////(&#/////////////////(&%***********#&/**&,,&&&&&&&&&&&&&&&&&%%*/&.,.,.@,,,,,,,&*&/,,,.&%  
   % (&.,,,,&.,,,,&*,&%***(%,...,,&(,,,,,,,,,,,,,,,,,**&&&&&&&&&&&&&//////////////////////((//#&&&(///////////&&&/***%&#*****/&&,,&/***&&&&&&&&&&&&&&&&&&((/,,,...@,,,,,,,,%#,,,,.&%  
   % (&.,,,,,,,,,,@,,&#***&,,.....&&,,,,,,,,,,,,,,,,#%///&&&&&&&&&//////////////////////////////#%**/#&&&&#(*******&&*&/****&&,,%%,,,,/&,,,&&&&&&&&&&/,*&,,,,@.,..,,,,,,,,,,,,,,,.&%  
   % (&.,,,,,,,,,/&,,&(**#%,,....,(&,,,,,,,,,,,,,,*#%/////%#((&%////////////////////////////////%&***************(&/*,&#***&%,,&*,,,,,*&*,,.,,,,,,,,,,,*&@&.,*&,,,,,,,,,,,,,,,,,,,&%  
   % (&.,,,,,,,,,(#,,&#**&,,,,#%,,/&*,,,,,,,,,,,,,,&//////&/(/&%////////////////////////////////(&**************/&/,,,%&**&#,,&*,,,,,,,&(,,.........,,,*&@&&,#&&&&&&&&&&&&&&&&&&&&&%  
   % (&.,,,,,,,,,(#,,%%*%#,.,,%&,,*&*,,,,,,,,,,,,*/&//////&//(&#///////&//(&/(///////////////////&&********,,***&&,.,,,&%&&,,/#*,,*,(&&&&,,...........,/&,,,&&******************,,&%  
   % (&.,,,,,,,,,*&,,,@,@...,.&%/,/&*,,,,,,,,,,,,*/&//////&((/&%///////&(//&%////////////////////&&&****,,,,,***&#,,*&#,&&/,,,,%&%*****(&.,...........,#&#%***********,,********,,&%  
   % (&.,,,,,,,,,,&,,,/&#,.,.%*##,/&*,,,,,,,,,,,,*/&//////&(//&&///////%#//&%////////////////////%&#%***,,,,,,**&#,,/%(&,,,,,&*******,**&*,,,,,,,,,..,,#%***,,,,***,,,,,,,,,,,,*,,&%  
   % (&.,,,,,,,,,,(%,,.,,..,/(*(#,/&,,,,,,,,,,,,,*/&#/////&((((&(//////(&///&(///////////////////&&/&,**,,,,,,**&#,,(%*,&(,,*&%,*****,***/((((((/**********,,,,,,,,,,,,,,,,,,,,,,,&%  
   % (&.,,,,,,,,,,,&/,....,(#**(%,#&,***,,,,*,****/&%(///(&((((%&///////&%///&#//////////////////&#/&/**,*,,,,*,&%,,*&,,,,#(/&**/%&&&%#/**,*************************************,,&%  
   % (&.,,,,,,,,,,,,&/,..,(#,,,##,%&,/***//*/*///*/&&(////%%((((%&#//////#&#/(/%&&#%%&&&&&%((/&%////&&&**,,,,,*,&&,,.@,,,,,,/&******,***********//(#%#%&&&&&&&&&&&&&&&&&&&&&&&&&&&&%  
   % (&.,,,,,,,,,,,,,&/,,(#,,,,%#,&&////*(/*/**///(/&%(/////%&((((&%#///////#&&&#////(////(/&&&(&&(////&(**,,,,*&&,,.@*,,,,,*&*****,,,,,,,,,,,,,,,,,,*****,,,********************/&%  
   % (&.,,,,,,,,,,,,,,#%,&,,,,,%/,&%/(/((((#%%%&%/(((%&&((/&(////////&&&#/////////#&&&&/////#&(((/&(///(%/%%#,**&&,,,&/,,,,,,,/&&&&%(**********,**,***(&(*,,,,,,,,,,,,,,,,,,,,*,&%&%  
   % (&.,,,,,,,%,,,,,,,*&%,,,,,&*,&&/((/((((((%(((((((((((%&///////////(//#&&&&&#///(%&//////&((((%#//&&(////&#*&&,,,%(,,,,,,,,,,,,,*,&&&&&&&&&&&&&&&&&&&**,******************(&(,&%  
   % (&.,,,,,,,*#,,,,,,,,,,,,,,&,,%&/(((//(((&&##(((//((((%&////////////////////////((&#////&&&(/(&%////////&@**&#,,,&*,,,,,,,,,,,,,,*&&&&%#(#&,,,**%&/(&**(&&&&&&&%%%%%&%#%#&&,..&%  
   % (&.,,,,,,,,@,,,,,,,,,,,,,,&*,*&((((((((&/(&((((/((((((&&&&%%#######%#%%&&&&&&&&&&&&///////(/%&&###%&&&(****&*,,.&,,,,,,,,,,,,,,,*&**,,**%&,*,%&(,,/&**%&&#&&%&,,,,,,,,@&,,,..&%  
   % (&.,,,,,,,,&,,,,,,,,,,,,,,%(,,&%/(/(/(&*,,&%(/((//(((((/((((&/(((&(%&&&&&&&%&&&&&&&%//(%&&&&&&&&&***,,,,,*%%,,,%/,,,,,,*(&/,,,,,*&**,,**%&,@&,,..,/&**%&/*#&&&,,,,,,#&*,.....&%  
   % (&.,,*.,,,,%/,,,,,,,,,,,,,/&,,*&((///&/,,,/&((//(/((((((((((&///(&(((((((((%&&&&&&(///@(((&&&&&&&/(/////((&,,.%/,,,,,,#&#&,,,,,,*&******&&(,,....,#&,*%(***&&&,,,*#%#,.....%#&%  
   % (&.,,*,,,,,(%,,,,,,,,,,,,,,&,,,#&(((%%,,#*,%&((///((((((((((&///(&#(//////(%&&&&&&////&(((&&&&&%&///////&(,,,,,,,,,&&((#&,,,,,,,*&/**&&(,,.......,&&,,%*,,,,,,,,#&(,.....*&*,&%  
   % (&.,,**,,,,/&,,,,,,,,,,,,,,,&.,,#&//&,,,&#,,&&/////(((//(/(/&#(//&#((((((((%&&&&&&////&(((&&&&&(&((///(/((((,#&&#(((((%&***,/&((/&&&(,,.........,.&#,,,,,,,,,,&&*,.....,%&,,.&%  
   % (&.,,,#,,,,,&,,,,,,,,,,,,,,,/&,,,,&(&,,*%&,,,&&///((%&&((/((%%((/&(((((((((#&&&&&&////&(((&&&&%#&/////(((((((((((((((&(&%(/(((/&&(,,.....#,.....,/&*,,,,,,,/&&,,.......&&,,,,&%  
   % (&.,,,&,,,,,%*,,,,*/,,,/,,,,,/&,,,,&&,,#/*&,,,&&/((%#&&((/((%&(((&(((/(((/(#&&&&&&////&(((&&&&(%&(////(((((((((((((((((//(((&&%,.....,,&*%*,....,(&,,,,,,&&(,.......,.&%,,,,.&%  
   % (&.,,,,#,,,,(%,,,,,&,,,(*,,,,,,&*,,....@,*,&,,,%&(#&,&%((/((#%(((&(((((((((#&&&&&&////&(((%&&&(%&/(///((((((%###((((((((((&@,,,,..,.%%,,,&,.....,(&,*,,&&,,..........&#,,,,,.&%  
   % (&.,,,,/,,,,,&,,,,,//,,,&,,,,,,,(&.,.,,%,,,,&/,,/&&(,&&////(%&/((&((((((((((&&&&&&////&(((#&&%(%&/(///(((((((((((((/((((&&.,....,.%(,,,,,%*,....,*&/*%&*,.........,/&*,,,,,,.&%  
   % (&.,,,,,(,,,,/%,,,,,&,,,,&,,,,,,,,*&/,%*,,,,,/%,,,,,,%&////(#&/((&((((((((((&&&&&&////&((((&&((%&((///((((((((((///(((&&.,.....,%(,,,,,,,*&.....,.&&&%,,........,.&%,,,,,,,,.&%  
   % (&.,,,,,,,,,,,%*,,,,,&,,,,%,,,,,,,,,,(#,,,,,,,,#&.,,,*&(///(#&&&&&((((((((((((((/&//((&((//&&&&&&/(///((((//((((((((&&,,.....,*&,,,,,,,,,,#(,.....,@*,........,/&**,,,,,,,,,.&%  
   % (&.,,,,,,,,,,,,&,,,,,,%,,,,,,,,,,,,,,,,,,,,,,,,,,(&*,,/%(//((((((%%%%&&&%%%((((((((((((((((((((((((((((((((((((//(/&%,,.....,%/,,,,,,,,,,,,(&,,............,,&#,,,,,,,,,,,,,.&%  
   % (&.,,,,,,,,,,,,,#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*,,,*/*/(((%%&&&&&&&&&&%%((((((((((((((((((%&&&&&&&&&%%%#(((///(&#,,......&,,,,,,,,,,,,,,,,&(,,.........&%,,,,,,,,,,,,,,,,.&%  
   % (&.,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,**,,,*(&&%((/(//(((((((((((//(((((((((((((((((/##%%%#&&&&%%((((((//(&#,,......%*,,,,,,,,,,,,,,,,,,&%.,...,(&,,,,,,,,,,,,,,,,,,,.&%  
   % (&.,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,(&&%#(((((/////((&&&&&%%%((((((((((((((((((((((((((((/((((((((((((((&&,,,....,/#,,,,,,,,,,,,,,,,,,,,,,&#.#&,,,,,,,,,,,,,,,,,,,,,.&%