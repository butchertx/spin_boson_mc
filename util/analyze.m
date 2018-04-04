%%make and save plots for all mc results in a given directory

%*****************************************
% Read full_input.txt to obtain parameters
%*****************************************
N = 12;%number of different param sets
params = SB_Params('full_input.txt', N)

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


% %%plot trajectories
% for i = 1:size(dump0,1)
%     dumpi = full_meas(i,:,:);
%     fig = figure();
%     hold on
%     for j = 1:N
%         plot(dumpi(1,:,j), 'DisplayName', num2str(j));
%     end
%     legend('show')
%     print(strcat('plots/traj',num2str(i)),'-dpng')
%     delete(fig)
% end

% %%plot autocorrelations
% for i = 1:N
%    dumpi = full_meas(:,:,i);
%    fig = figure();
%    hold on
%    for j = 1:size(dumpi, 1)
%        plot(autocorrelation(dumpi(j,:,1)), 'DisplayName', num2str(j));
%    end
%    legend('show')
%    print(strcat('plots/autocorr',num2str(i)),'-dpng')
%    delete(fig);
% end

%************************************
% Read correlation functions and plot
%************************************
full_corr_nosplit = zeros(params.lx*params.ly, N);
full_corr_split = zeros(params.ly,params.lx, N);
%make plot for each alpha with different line for each site
for i = 1:N
    full_corr_nosplit(:,i) = csvread(strcat('correlation',num2str(i - 1),'.csv'), 1,0);
    full_corr_split(:,:,i) = split_corr(full_corr_nosplit(:,i),params.lx, params.ly);
    fig = figure();
    hold on
    for j = 1:params.lx
        plot(full_corr_split(:,j,i),'DisplayName',num2str(j));
    end
    legend('show')
    print(strcat('plots/split_corr',num2str(i)),'-dpng')
    delete(fig)
end

%make plot for spatial correlation functions
r_corr = zeros(params.lx, N);
for i = 1:N
   norm = sum(full_corr_split(:, 1, i));
   for x = 1:params.lx
       r_corr(x, i) = 1/norm * sum(full_corr_split(:, x, i));
   end
end
fig = figure();
hold on
for i = 1:N
    plot(r_corr(:, i), 'DisplayName', strcat('alpha ', num2str(i)))
end
print('plots/r_corr', '-dpng')
delete(fig)

%********************************
% Read results data file and plot
%********************************
fid = fopen('results.dat');
results = zeros(N,12);
tline = fgetl(fid);%column labels
for i = 1:N
    tline = fgetl(fid);
    C = textscan(tline, '%f');
    results(i, :) = C{1}';
end
%plot loc2 and mag2
figure()
hold on
plot(results(:,1), results(:, 3))
legend('\langle l^2 \rangle')
print('plots/results','-dpng')

%********************************
% Resampling
%********************************
% 
% resample_alpha_range = [1.0 1.02 1.05 1.09 1.14 1.19 1.25 1.30 1.35 1.42 1.49 1.56];
% resample_loc2 = zeros(length(resample_alpha_range), N);
% 
% for i = 1:N
%     [S1, S2] = split_action(full_meas(7, 1:end, i), params.alpha(i), params.gamma, full_meas(8, 1:end, i), params.ly, params.lx);
%     resample_loc2(:, i) = swendsen_resample_point(params.alpha(i), S2, resample_alpha_range, full_meas(3,1:end,i));
% end
% plot(resample_alpha_range, resample_loc2, '-x')
















