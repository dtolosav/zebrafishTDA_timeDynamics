function [handle] = plot_bettigraph_Dan(intervals,dimension,right_bound)

if ~exist('right_bound','var')
    % third parameter does not exist, so default it to something
    right_bound = 1200;
end

r_bound = [right_bound];

x_min=0;
endpoints= intervals(:,2);
num_increments = 500;
threshold = max(endpoints(~isinf(endpoints)))*(7/6);
increment = threshold/num_increments;
handle = figure;
hold on;

x_values = x_min : increment : threshold;
betti = zeros(1, num_increments+1);

intervals(:,2)

for i = 1 : length(x_values)
    betti(i) = nnz(intervals(:,1) <= x_min + increment*i & intervals(:,2) >= x_min + increment*i ) ;
end

nnz(betti)

%if dimension == 1
%    plot(cat(1,x_values,r_bound), cat(1,betti,[1]));
%
%else
%    plot(x_values , betti);
%end

plot(x_values,betti, "LineWidth", 2);
if dimension == 0
    title('$B(\epsilon)$ of Betti 0', "Interpreter", "Latex", "FontSize", 25);
end
if dimension == 1
    title('$B(\epsilon)$ of Betti 1', "Interpreter", "Latex", "Fontsize", 25);
end
ylabel('Betti', "FontSize", 25);
xlabel('epsilon', "FontSize", 25);
ax.FontSize=25;
%ylim([-1, min(10, max(betti)/10)]);
    
hold off;

end
    