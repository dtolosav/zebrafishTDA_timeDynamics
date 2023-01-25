% This is a MATLAB function that takes a matrix of barcodes (the output of Ripser) and 
% plots the corresponding barcodes.
%
% Inputs: barcodes matrix, dimension of betti;
% 
% Outputs: Barcode figure.
% Author: Daniel Tolosa, Purdue University, 2022. Based on work by Melissa R. McGuirl, Brown University, 2019.

function [handle] = plot_barcodes_Dan(intervals, dimension)
num_intervals = length(intervals(:,1));

line_width = 1.5;
x_min=0;
x_max= max(intervals(:,2));
epsilon=1e-6;
point_width = 0.006 * (350 - 0);
endpoints = intervals(:,2);
x_max_int = max(endpoints(~isinf(endpoints)));
threshold = min(x_max_int *(10/9), 700);

% dimension = 0;

handle = figure;
hold on;
for i = 1:num_intervals
    start = intervals(i, 1);
    finish = intervals(i, 2);
            y = num_intervals - i + 1;
            
            if (finish >= threshold && start <= -threshold)
                line([x_min, x_max], [y, y], 'LineWidth', line_width);
                line([x_min, x_min], [y, y], 'Marker', '<', 'LineWidth', line_width);
                line([x_max, x_max], [y, y], 'Marker', '>', 'LineWidth', line_width);
            end
            
            if (finish >= threshold && start > -threshold)
                line([start, x_max_int *(7/6)], [y, y], 'LineWidth', line_width);
                line([x_max_int *(7/6), x_max_int *(7/6)], [y, y], 'Marker', '>', 'LineWidth', line_width);
            end
            
            if (finish < threshold && start <= -threshold)
                line([x_min, finish], [y, y], 'LineWidth', line_width);
                line([x_min, x_min], [y, y], 'Marker', '<', 'LineWidth', line_width);
            end
            
            if (finish < threshold && start > -threshold)
                if (abs(finish - start) < epsilon)
                    line([start - 0.5 * point_width, finish + 0.5 * point_width], [y, y], 'LineWidth', line_width);
                else
                    line([start, finish], [y, y], 'LineWidth', line_width);
                end
            end

        axis([x_min, threshold, 0, num_intervals + 1]);

if dimension == 0
    ylabel("connected components","FontSize", 25);
end

if dimension == 1
    ylabel("loops (elements of $H_1^\epsilon$)", "Interpreter", "Latex","FontSize", 25);
end
caption = 'Barcode of Betti ';
           title(sprintf('%s %d', caption, dimension));
xlabel("$\epsilon$ values in micrometers", "Interpreter", "Latex","FontSize", 25);
ax.FontSize=25;

 end
    
    hold off;

end
    