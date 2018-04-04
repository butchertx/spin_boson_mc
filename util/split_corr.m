function [ new_corr ] = split_corr( corr, lx, ly )
%SPLIT_CORR split a vector into a matrix of dimension ly x lx

new_corr = zeros(ly, lx);
for x = 1:lx
    new_corr(:,x) = corr((x-1)*ly + 1 : x*ly);
end

end

