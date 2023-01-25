function [curvatures] = curvatures(sim_number)
    if ~exist('sim_number','var')
        % if parameter does not exist, default it to something
        sim_number = 808;
    end
% This is a Matlab script that computes stripe straightness of a single simulation, using Xc cells.
% Output: A vector that has the average curvature of stripes at each moment in time, and a plot of this information.
path= sprintf('WT_default/Out_WT_default_%d', sim_number);
load(path);

% load in bar codes of mel, xanc and xanl. The first coordinate corresponds to cell type (mel, xc, xl)
barcodes=cell(3,46);
    for i = 1:46
        temp_file = sprintf('barcodes/sim%d/Mel/PD_Melsim%dtime%d_dim1', sim_number, sim_number, i);
        [barcodes{1,i}, ~, ~]= importdata(temp_file);
    end
    for i = 2:46
        temp_file = sprintf('barcodes/sim%d/XanC/PD_XanCsim%dtime%d_dim1', sim_number, sim_number, i);
        [barcodes{2,i}, ~, ~]= importdata(temp_file);
    end    
    for i = 1:46
        temp_file = sprintf('barcodes/sim%d/XanSn/PD_XanSnsim%dtime%d_dim1', sim_number, sim_number, i);
        [barcodes{3,i}, ~, ~]= importdata(temp_file);
    end

%  get dimension 1 persistence; death-birth

mel_pers = cell(46,1);
xanC_pers = cell(46,1);
xanSn_pers = cell(46,1);

for i = 2:46
temp_mel = cell2mat(barcodes(1,i));
mel_pers{i}= temp_mel(:,2)-temp_mel(:,1);
temp_xanC = cell2mat(barcodes(2,i));
xanC_pers{i} = temp_xanC(:,2)-temp_xanC(:,1);
temp_xanSn = cell2mat(barcodes(3,i));
xanSn_pers{i} = temp_xanSn(:,2)- temp_xanSn(:,1);
end
% define persistence thresholds based on cell-to-cell measurements
pers_cutoff = 200;

% compute dimension 1 betti numbers

b1_mel = zeros(46);
b1_xanC = zeros(46);
b1_xanSn = zeros(46);

for i =2:46
b1_xanC(i) = length(find(xanC_pers{i} > pers_cutoff  & barcodes{2,i}(:,1) < 80));
b1_xanSn(i) = length(find(xanSn_pers{i} > pers_cutoff & barcodes{3,i}(:,1) < 100));
b1_mel(i) = length(find(mel_pers{i} > pers_cutoff & barcodes{1,i}(:,1) < 90));

end

num_stripes = b1_xanSn; 
num_Istripes = b1_xanC;


% This break part only works for fully formed patterns. 
%stripe_breaks = zeros(46,1);
%Istripe_breaks = zeros(46,1);
%for i =2:46
    %if b1_xanSn(i) < 2 && b1_mel(i) < 2
    %    stripe_breaks(i) = 1;
   % else 
     %   stripe_breaks(i) = 0;
    %end

    %if b1_xanC(i) < 3 
    %    Istripe_breaks(i) = 1;
   % else 
  %      Istripe_breaks(i) = 0;
 %   end
%end


% Use X^C cells to compute straightness measure of stripes 

% num_stripes is just betti 1 of xanC
num_stripes = zeros(46,1);
    for j = 2 : 46
        num_stripes(j) = length(find(barcodes{2,j}(:,1) < 150 & barcodes{2,j}(:,2) > 210));
    end

curvatures = zeros(46,1);
x_querys = cell(46,1);

% Now we clean the cell data of the weird -40000 cells. If we 
% don't do this our clustering goes crazy

cutoff4cells = 0.1*boundaryY ;

cells_xanC = [];
    for i = 1 : 46
        cells_xanC{i} = cellsXc(find(cellsXc(1:numXanc(i), 2,i) > cutoff4cells(i) &  cellsXc(1:numXanc(i), 2, i) < boundaryY(i) - cutoff4cells(i)), :, i);
    end

for i =2 :46
    x_querys{i} = 0:60:boundaryX(i);
    figname = sprintf("Quantify Stripes/Stripe Curvature/sim%d/time%d",sim_number,i);
    [top_cv, bottom_cv] = straightness_measure_Dan(cells_xanC{i}, x_querys{i},num_stripes(i) , i , b1_xanSn(i), b1_mel(i), figname, 1);
    curvatures(i) = mean([top_cv; bottom_cv]);
end

% plot curvatures
f= figure;
plot(curvatures,'LineWidth',2);
xlim([0 46]);
%axis equal;
titulo= sprintf("plot mean curvatures sim%d",sim_number);
title(titulo);
filename = sprintf("Quantify Stripes/Stripe Curvature/sim%d/curvature.mat", sim_number);
save(filename,"curvatures");

end