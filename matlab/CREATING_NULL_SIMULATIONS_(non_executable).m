%							CREATING NULL SIMULATIONS
%
% Creating the NULL model (50 null simulations, with 46 time points each). We only need to keep the barcodes.
% 
% First figure out how many cells should go on each day. We do this by averaging over the number of cells on the 50 Volkening
% simulations. WE make sure that this is good enough by computing the variance and making sure it's small.

% Step 1: Load barcodes for simulations 1-50 (wild type, XanC), compute averages and variance of persistences for each day.
	cells_xan = cell(50,46);
	for j = 1 : 50;
	file_name = sprintf("WT_default/Out_WT_default_%d", j);
	load(file_name, "cellsXc", "boundaryX" , "boundaryY", "numXanc");
	% clean cell coordinates for XanC cells, cutting the domain, we will work on these cells in the remainder
	cutoff = 0.1*boundaryY ;
		for i = 2 : 1 : 46
			cells_xan{j,i} = cellsXc(find(cellsXc(1:numXanc(i), 2,i) > cutoff(i) &  cellsXc(1:numXanc(i), 2, i) < boundaryY(i) - cutoff(i)), :, i);
		end
	end
	% create a matrix M (46 x 2) where each row corresponds to a pair (mean # of cells, standard deviation # of cells)
	M = zeros(46,4);
	for i = 2: 46
		temp_vector = zeros(50,1);
		for j = 1 : 50
			temp_vector(j) =  length(cells_xan{j,i}) ;
		end
		M(i,1) = mean(temp_vector);
		M(i,2)=  sqrt(var(temp_vector)) ;
		M(i,3) = max(temp_vector);
		M(i,4) = min(temp_vector);
	end

% Step 2: Create 50 cell-coordinate matrices, one for each simulation. For a fixed simulation S, we decide the number 
% 			of cells per day by a normal distribution given by the data from the 50 Volkening simulations.
	num_cells = zeros(50, 46); % matrix containing the number of cells with columns (sim_number, day)
	for j =1 :50
		for i =2:46
			num_cells(j,i) = round(normrnd(M(i,1),M(i,2)));
		end
	end
	% now creating the cell-coordinate matrices
	NULL_Sims = cell(50,46);
	for j = 1 : 50
		for i = 1 : 46
			temp_x_coordinates =  (boundaryX(i)).*rand(num_cells(j,i),1);
			temp_y_coordinates =  cutoff(i) + (boundaryY(i) - 2*cutoff(i)).*rand(num_cells(j,i),1);
			NULL_Sims{j,i} = [temp_x_coordinates temp_y_coordinates];			
		end
	end
% Step 3: Save each simulation
for j = 1:50
	NULL_sim_temp = cell(1, 46);
	for i=1:46
		NULL_sim_temp{i} = NULL_Sims{j,i};
	end
	filename = sprintf("NULL_simulations/NULL_sim_%d.mat",j);
	save(filename,"NULL_sim_temp");
end

% Step 4: Compute and save barcodes
	barcodes_betti0 = cell(50,46);
	barcodes_betti1 = cell(50,46);
	Distance_matrices=cell(50,46);
	for j = 1:50 % fix a simulation
		for i = 2 : 46 
			[Distance_matrices{j,i}, ~, ~] = getPeriodicDistMats(NULL_Sims{j,i}, boundaryX(i));
			Distance_matrices_py = py.numpy.array(Distance_matrices{j,i});
			barcodes_temp = py.rip_on_ML.main(Distance_matrices_py);
			barcodes_temp_0 = barcodes_temp(1);
			barcodes_temp_1 = barcodes_temp(2);
			barcodes_betti0{j,i} = double(barcodes_temp_0{1});
			barcodes_betti1{j,i} = double(barcodes_temp_1{1});
			temp_var_0 = barcodes_betti0{j,i};
			temp_var_1 = barcodes_betti1{j,i};
			% save barcodes for sim j time i
			fname_0 = sprintf("NULL_simulations/barcodes/sim%d/BC_NULL_XanC_sim%d_time%d_dim%d", j,j,i,0);
			fname_1 = sprintf("NULL_simulations/barcodes/sim%d/BC_NULL_XanC_sim%d_time%d_dim%d", j,j,i,1);
			save(fname_0,"temp_var_0");
			save(fname_1,"temp_var_1");
			%{
fileID = fopen(fname_0,'w');
			fprintf(fileID,barcodes_betti0{j,i});
			fclose(fileID);
			fileID2 = fopen(fname_1,'w');
			fprintf(fileID2,barcodes_betti1{j,i});
			fclose(fileID2);
%}

		end
	end