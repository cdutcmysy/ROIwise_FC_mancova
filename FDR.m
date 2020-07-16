function [significant_index, final_threshold] = FDR(p, alpha)

% [significant_index] = FDR(p, alpha)
%
%

[p, I] = sort(p);
threshold = alpha/length(p);

final_threshold = 0;
significant_index = [];
for i = 1:length(p)
    if(p(i) <= i * threshold)
        final_threshold = i * threshold;
        significant_index = [significant_index; I(i)];
    else
       break; 
    end
end