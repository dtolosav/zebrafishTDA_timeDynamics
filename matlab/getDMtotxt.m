
% Author: Daniel Tolosa, Purdue University, 2022. Based on work by Melissa R. McGuirl, Brown University, 2019.

%----------------------------------- CHOOSE WHETHER USER INPUTS SIM NUMBER OR NOT  ---------------------------------

% MAKE SURE ALL THE NEEDED MODULES/FUNCTIONS ARE IN MATLAB PATH

function [] = getDMtotxt(sim_number)
if ~exist('sim_number','var')
    % if parameter does not exist, default it to something
    sim_number = 808;
end
file = sprintf('WT_default/Out_WT_default_%d', sim_number);

% --------------------------------------------------------------------------------------------------

load(file);

cd E:\Math\Zebrafish\Distance_Matrices;
dirname=sprintf("sim%d",sim_number);
mkdir(dirname);
cd(dirname);

% we loaded all the cell coordinates of all cell types at each time, from 1 to 46

% select only the Mel cells that live inside the boundary at each time from 1 to 46

cutoff = 0.1*boundaryY ;

% 	WORKING WITH MEL

	cells_mel = [];

	mkdir Mel;
	cd Mel;
	for i = 1 : 46
		cells_mel{i} = cellsM(find(cellsM(1:numMel(i), 2,i) > cutoff(i) &  cellsM(1:numMel(i), 2, i) < boundaryY(i) - cutoff(i)), :, i);
	end

	% Get periodic distance matrix D_mel . 

	D_mel=cell(46,1);
	for i = 1 : 3
		[D_mel{i}, ~, ~] = getPeriodicDistMats(cells_mel{i}, boundaryX(i));
		%D_mel_py = py.numpy.array(D_mel{i});
		ssavefile= sprintf("Dist_Matrix_Mel_sim%d_time%d.txt", sim_number, i);
		writematrix(D_mel{i}, ssavefile);
		%type ssavefile;
	end

	% Save periodic distance matrix to .m file

	sf= sprintf("Dist_Matrix_Mel_sim%d.mat", sim_number);
	save(sf,  "D_mel");

cd ..\

% WORKING WITH XAN C
	mkdir XanC;
	cd XanC;
	cells_XanC = [];
	for i = 1 : 46
		cells_XanC{i} = cellsXc(find(cellsXc(1:numXanc(i), 2,i) > cutoff(i) &  cellsXc(1:numXanc(i), 2, i) < boundaryY(i) - cutoff(i)), :, i);
	end

	% Get periodic distance matrix D_xanc . 

	D_xanc=cell(46,1);
	for i = 1 : 3
		[D_xanc{i}, ~, ~] = getPeriodicDistMats(cells_XanC{i}, boundaryX(i));
		sssavefile= sprintf("Dist_Matrix_XanC_sim%d_time%d.txt", sim_number, i);
		writematrix(D_xanc{i}, sssavefile);
	end

	% Save periodic distance matrix to .m file

	sf_xanc= sprintf("Dist_Matrix_XanC_sim%d.mat", sim_number);
	save(sf_xanc,  "D_xanc");

cd ..\
 
%WORKING WITH XAN SN
	mkdir XanSn;
	cd XanSn;
	cells_XanSn = [];
	for i = 1 : 46
		cells_XanSn{i} = cellsXsn(find(cellsXsn(1:numXansn(i), 2,i) > cutoff(i) &  cellsXsn(1:numXansn(i), 2, i) < boundaryY(i) - cutoff(i)), :, i);
	end

	% Get periodic distance matrix D_xansn . 

	D_xansn=cell(46,1);
	for i = 1 : 3
		[D_xansn{i}, ~, ~] = getPeriodicDistMats(cells_XanSn{i}, boundaryX(i));
		avefile= sprintf("Dist_Matrix_XanSn_sim%d_time%d.txt", sim_number, i);
		writematrix(D_xansn{i}, avefile);
	end

	% Save periodic distance matrix to .m file

	sf_xansn= sprintf("Dist_Matrix_XanSn_sim%d.mat", sim_number);
	save(sf_xansn,  "D_xansn");

cd ..\


% WORKING WITH IRI DENSE
	mkdir IriD;
	cd IriD;
	cells_IriD = [];
	for i = 1 : 46
		cells_IriD{i} = cellsId(find(cellsId(1:numIrid(i), 2,i) > cutoff(i) &  cellsId(1:numIrid(i), 2, i) < boundaryY(i) - cutoff(i)), :, i);
	end

	% Get periodic distance matrix D_xanc . 

	D_irid=cell(46,1);
	for i = 1 : 3
		[D_irid{i}, ~, ~] = getPeriodicDistMats(cells_IriD{i}, boundaryX(i));
		vefile= sprintf("Dist_Matrix_IriD_sim%d_time%d.txt", sim_number, i);
		writematrix(D_irid{i}, vefile);
	end

	% Save periodic distance matrix to .m file

	sf_irid= sprintf("Dist_Matrix_IriD_sim%d.mat", sim_number);
	save(sf_irid,  "D_irid");

cd ..\

% WORKING WITH IRI LOOSE
	mkdir IriL;
	cd IriL;
	cells_IriL = [];
	for i = 1 : 46
		cells_IriL{i} = cellsIl(find(cellsIl(1:numIril(i), 2,i) > cutoff(i) &  cellsIl(1:numIril(i), 2, i) < boundaryY(i) - cutoff(i)), :, i);
	end

	% Get periodic distance matrix D_xanc . 

	D_iril=cell(46,1);
	for i = 1 : 3
		[D_iril{i}, ~, ~] = getPeriodicDistMats(cells_IriL{i}, boundaryX(i));
		efile= sprintf("Dist_Matrix_IriL_sim%d_time%d.txt", sim_number, i);
		writematrix(D_iril{i}, efile);
	end

	% Save periodic distance matrix to .m file

	sf_iril= sprintf("Dist_Matrix_IriL_sim%d.mat", sim_number);
	save(sf_iril,  "D_iril");

cd ..\..\