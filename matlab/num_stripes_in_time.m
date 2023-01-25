% This is a MATLAB function that reads in .txt files containing the
% barcodes of an agent-based model for zebrafish in the form of a vector and counts stripes 
% as a function of time.
%
% Inputs:number of simulation;
% 
% Outputs: list of number of stripes for times 1 through 46.
% Author: Daniel Tolosa, Purdue University, 2022.
%----------------------------------- CHOOSE WHETHER USER INPUTS SIM NUMBER OR NOT  ---------------------------------
function [num_stripes, barcodes] = num_stripes(sim_number,cell_type)
	if ~exist('sim_number','var')
	    % if parameter does not exist, default it to something
	    sim_number = 808;
	end
	%tic
	barcodes=cell(46,1);

	for i = 1:46
		temp_file = sprintf('WT_barcodes/sim%d/%s/PD_%ssim%dtime%d_dim1', sim_number,cell_type , cell_type,sim_number, i);
		[barcodes{i}, ~, ~]= importdata(temp_file);
	end
	% --------------------------------------------------------------------------------------------------

	num_stripes = zeros(46,1);
	for j = 1 : 46
	    num_stripes(j) = num_stripes_diagonal_distance(barcodes{j}, j, 0 );
	end

	%num_stripes
	%toc
end