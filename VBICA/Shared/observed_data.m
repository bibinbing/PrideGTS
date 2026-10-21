function [values, observed] = observed_data(data)
% NaN denotes an unobserved entry; zeros here are algebraic placeholders.
% Every likelihood term using VALUES must also use OBSERVED.
observed = ~isnan(data);
values = data;
values(~observed) = 0;
end
