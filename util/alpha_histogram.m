function [ output_args ] = alpha_histogram( dumpdir, num_sims, num_bins )
%Plot a histogram of the actions for a given parallel simulation

figure()
hold on
for i=1:num_sims
    tempdat = csvread(strcat(dumpdir,'dump',num2str(i-1),'.csv'),0,1);
    if(nargin==2)
        histogram(tempdat(13,:)');
    else
        histogram(tempdat(13,:)', num_bins);
    end
end
print('histogram','-dpng')


end

