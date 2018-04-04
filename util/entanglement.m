%%Plot entanglement entropy given some dump files

%*****************************************
% Read full_input.txt to obtain parameters
%*****************************************
N = 12;%number of different param sets
params = SB_Params('full_input.txt', N);

%*************************************************************
% Read all dump files and plot trajectory and autocorrelations
%*************************************************************
%read all dump files
dump0 = csvread('dump0.csv', 0, 1);
%fid = fopen('dump0.csv');
%data = textscan(fid,'%f','Delimiter',',','EndOfLine','\r\n');
%dump0 = data{1};
%fclose(fid);
full_meas = zeros(size(dump0, 1), size(dump0, 2), N);
size(full_meas)
for i = 1:N
    full_meas(:, :, i) = csvread(strcat('dump',num2str(i - 1),'.csv'), 0, 1);
    %fid = fopen(strcat('dump',num2str(i - 1),'.csv'));
    %data = textscan(fid,'%f','Delimiter',',','EndOfLine','\r\n');
    %full_meas(:, :, i) = data{1};
    %fclose(fid);
end

sx = full_meas(