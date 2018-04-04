%%Plot entanglement entropy given some dump files
N = 12;
fid = fopen('results.dat');
results = zeros(N,12);
tline = fgetl(fid);%column labels
for i = 1:N
    tline = fgetl(fid);
    C = textscan(tline, '%f');
    results(i, :) = C{1}';
end

sx = results(:,11);
sx_avg = mean(sx,2);
sz = results(:,2);
sz_avg = mean(sz,2);
p_plus = 0.5*(1 + sqrt(sx_avg.*sx_avg + sz_avg.*sz_avg));
p_minus = 0.5*(1 - sqrt(sx_avg.*sx_avg + sz_avg.*sz_avg));
E = squeeze(-p_plus.*log2(p_plus) - p_minus.*log2(p_minus));
plot(results(:,1), E)