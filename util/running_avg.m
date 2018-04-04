function avgs = running_avg(data, bin_size)
%%compute running averages of samples held in the columns of 'data'
%%bin_size is the number of samples to average for each point in avgs

[rows, cols] = size(data);
avg_rows = floor(rows/bin_size);
if(avg_rows * bin_size ~= rows)
    disp('error: in "running_avg", bin size does not divide # of samples')
    disp(bin_size)
    disp(rows)
end
avgs = zeros(avg_rows,cols);

for i = 1:cols
    for j = 1:avg_rows
        avgs(j,i) = mean(data(((j-1)*bin_size + 1):(j*bin_size),i));
    end
end