%%make and save plots for all mc results in a given directory, using
%%multiple simulations with the same parameters
%%Directory structure:
%   1-----------------
%   ----full_input.txt
%   ----mc_mpi.slurm
%   ----results.dat
%   ----dump/---------
%       ----dump0.csv etc.
%   2-----------------
%   .
%   .
%   .

%*****************************************
% Read full_input.txt to obtain parameters
%*****************************************
N = 12;%number of different param sets
dir = 10; %number of directories
params = SB_Params('1/full_input.txt', N)

%*************************************************************
% Read all dump files and plot trajectory and autocorrelations
%*************************************************************
%read all dump files
dump0 = csvread('1/dump/dump0.csv', 0, 1);
full_meas = zeros(size(dump0, 1), size(dump0, 2)*dir, N);
size(full_meas)
for i = 1:N
    for j = 1:dir
        full_meas(:, (j - 1)*size(dump0, 2) + 1 : j*size(dump0, 2), i) = csvread(strcat(num2str(j), '/dump/dump',num2str(i - 1),'.csv'), 0, 1);
    end
end

%%plot trajectories
for i = 1:size(dump0,1)
    dumpi = full_meas(i,:,:);
    fig = figure();
    hold on
    for j = 1:N
        plot(dumpi(1,:,j), 'DisplayName', num2str(j));
    end
    legend('show')
    print(strcat('plots/traj',num2str(i)),'-dpng')
    delete(fig)
end
% 
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
    for j = 1:dir
        full_corr_nosplit(:,i) = full_corr_nosplit(:, i) + csvread(strcat(num2str(j), '/dump/correlation',num2str(i - 1),'.csv'), 1,0)';
    end
    full_corr_nosplit = full_corr_nosplit/dir;
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
%make plot for integrated distance dependence with all alpha on same plot
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
full_results = zeros(N,12);
for j = 1:dir
    results = zeros(N, 12);
    fid = fopen(strcat(num2str(j),'/results.dat'));
    tline = fgetl(fid);%column labels
    for i = 1:N
        tline = fgetl(fid);
        C = textscan(tline, '%f');
        results(i, :) = C{1}';
    end
    full_results = full_results + results/dir;
end
%plot loc2 and mag2
plot(full_results(:,1), full_results(:,3), full_results(:,1), full_results(:, 7))
legend('loc2','mag2')
print('plots/results','-dpng')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%Read full_input.txt to obtain parameters
% N = 12;%number of different param sets
% ndir = 1;%number of simulation directories
% params = SB_Params('onesite1/full_input.txt', N)
% 
% %%Read all output files to obtain and organize data
% %read all dump files
% dump0 = csvread('onesite1/dump/dump0.csv', 0, 1);
% full_meas = zeros(size(dump0, 1), size(dump0, 2), N, ndir);
% size(full_meas)
% for dir = 1:ndir
%     for i = 1:N
%         full_meas(:, :, i, dir) = csvread(strcat('onesite', num2str(dir),'/dump/dump',num2str(i - 1),'.csv'), 0, 1);
%     end
% end
% size(full_meas(3, :, 1, :))
% for i = 1:N
%     loc2(i) = mean(mean(full_meas(3, :, i, :)));
% end
% figure()
% hold on
% for i = 1:ndir
%     plot(full_meas(3,:,1,i))
% end
% legend('1','2','3','4','5','6','7','8','9','10')
% %%Plot results
% figure()
% hold on
% plot(params.alpha, loc2)
% xlabel('alpha')
% ylabel('loc2')
% 
% %%Plot results with resampling
% resample_alpha_range = 0.55:0.001:0.9;
% resample_loc2 = zeros(length(resample_alpha_range), N);
% %full_meas = cat_measures(full_meas);
% % size(full_meas)
% % figure()
% % plot(full_meas(3,:,3))
% % plot(full_meas(7,:,3))
% for i = 1:N
%     [S1, S2] = split_action(full_meas(7, 1:1000, i,1), params.alpha(i), params.gamma, full_meas(8, 1:1000, i,1), params.ly, params.lx);
%     resample_loc2(:, i) = swendsen_resample_point(params.alpha(i), S2, resample_alpha_range, full_meas(3,1:1000,i,1));
% end
% figure()
% plot(resample_alpha_range, resample_loc2, '-x')



%%Fit correlations
% full_corrs = zeros(params.lx*params.ly, N, ndir);
% for dir = 1:ndir
%     figure()
%     hold on
%     for i = 1:N
%         full_corrs(:,i,dir) = csvread(strcat('onesite',num2str(dir),'/dump/correlation',num2str(i-1),'.csv'), 1,0)';
%         plot(full_corrs(:,i,dir))
%     end
% end
% avg_corrs = zeros(params.lx*params.ly, N);
% for i = 1:N
%     avg_corrs(:,i) = mean(full_corrs(:,i,:),3);
% end
% mat_corrs = zeros(params.ly,params.lx,N);
% x_corrs = zeros(params.lx, N);
% for i = 1:N
%     mat_corrs(:,:,i) = split_corr(avg_corrs(:,i),params.lx, params.ly);
%     x_corrs(:,i) = integrate_corr(avg_corrs(:,i),params.lx,params.ly);
% end
% figure()
% plot(x_corrs)