function sim_data = plot_sim( input_type, input_data, varargin )
%Make plots of the results of a single simulation, possibly including
%multiple threads/MPI
%input_type is 'dir' or 'data'
%input_data is a simulation directory name (for 'dir') or a results data
%array for ('data')
%varargin is a list of plots to make

if(strcmp(input_type,'dir'))
    params = SB_Params(strcat(input_data,'../full_input.txt'));
    %create cell array for all dump data
    sim_data = cell(params.n_sims, 1);
    size(sim_data)
    sim_corrs = cell(params.n_sims, 1);
    sim_vort_corrs = cell(params.n_sims, 1);
    sim_ints = cell(params.n_sims,1);
    sim_results = csvread(strcat(input_data,'results.csv'),0,1);
    for n = 1:params.n_sims
        sim_data{n} = csvread(strcat(input_data,'dump',num2str(n-1),'.csv'),0,1);
        sim_corrs{n} = csvread(strcat(input_data,'correlation',num2str(n-1),'.csv'),1,0);
        sim_ints{n} = csvread(strcat(input_data,'interactions',num2str(n-1),'.csv'),1,0);
        sim_vort_corrs{n} = csvread(strcat(input_data,'vort_corr',num2str(n-1),'.csv'),1,0);
    end
    i = 1;
    while i < length(varargin)
        %first is a word giving the type of plot to make: traj, avg, calc,
        %corr, or vortcorr
        %for traj, avg, and calc, next arg is a number giving the number of
        %observables (nobv) to plot
        %the next nobv args are the names of observables. The next arg is
        %the number of samples to skip, if any
        argin = varargin{i};
        i = i + 1;
        if(strcmp(argin,'traj'))
            nobv = varargin{i};%must be >= 1
            if(nobv == 1)
                obvs = {varargin{i+1:i+nobv}};
            else
                obvs = varargin{i+1:i+nobv};
            end
            num_skip = varargin{i + nobv + 1};
            i = i + nobv + 2; %now points to next element after traj specifiers
            celldisp(obvs)
            for k = 1:length(obvs)
                if(strcmp(obvs{k},'loc2'))
                    plot_mat = zeros(size(sim_data{1}(2,:)'));
                    size(plot_mat)
                    for n = 1:params.n_sims
                        plot_mat(:,n) = sim_data{n}(2,:)';
                    end
                    figure()
                    plot(plot_mat);
                end
            end
            
        elseif(strcmp(argin,'avg'))
            
        elseif(strcmp(argin,'calc'))
            
        elseif(strcmp(argin,'int'))
            %next arg is a vector of the picks and the next is 'full' or 'integrated'
            picks = varargin{i};
            plot_type = varargin{i+1};
            i = i + 2;
            if(strcmp(plot_type,'integrated'))
                for p = 1:length(picks)
                    temp_mat = zeros(params.lx,1);
                    for x = 1:params.lx
                        temp_mat(x) = sum(sim_ints{picks(p)}(params.ly*(x-1)+1:params.ly*x));
                    end
                    figure()
                    plot(-temp_mat)
                end
            end 
            
        elseif(strcmp(argin,'corr') || strcmp(argin,'vortcorr'))
            %first arg after 'corr' gives 'real' or 'fourierkw'
            %next gives a vector of the simulation indices to plot
            %then a type: 'planar', 'linear', or 'integrated'
            %if 'integrated' is chosen, the next argument is the
            %integration dimension.
            space = varargin{i};
            picks = varargin{i+1};
            plot_type = varargin{i+2};
            if(strcmp(plot_type,'integrated'))
                int_dim = varargin{i+3};
                i = i + 4;
            else
                i = i + 3;
            end
            %make into matrices
            corr_mats = cell(length(picks),1);
            for p = 1:length(picks)
                temp_mat = zeros(params.ly,params.lx);
                for x = 1:params.lx
                    if(strcmp(argin,'corr'))
                        temp_mat(:,x) = sim_corrs{picks(p)}(params.ly*(x-1)+1:params.ly*x)';
                    else
                        temp_mat(:,x) = sim_vort_corrs{picks(p)}(params.ly*(x-1)+1:params.ly*x)';
                    end
                end
                corr_mats{p} = temp_mat;
            end            
            for p = 1:length(picks)
                if(strcmp(space,'fourierkw'))
                    corr_mats{p} = abs(fft2(corr_mats{p}));
                elseif(~strcmp(space,'real'))
                    disp('No valid space chosen for correlation plot')
                end
                if(mod(size(corr_mats{p},1),2) == 0)
                    corr_mats{p} = vertcat(corr_mats{p},corr_mats{p}(1,:));
                end
                if(mod(size(corr_mats{p},2),2) == 0)
                    corr_mats{p} = horzcat(corr_mats{p}, corr_mats{p}(:,1));
                end
                if(strcmp(plot_type,'linear'))
                    
                elseif(strcmp(plot_type,'planar'))
                    figure()
                    if(strcmp(space,'fourierkw'))
                        imagesc([0,2],[0,2],real(corr_mats{p}))
                    else
                        imagesc(corr_mats{p})
                    end
                    colorbar
                elseif(strcmp(plot_type,'integrated'))
                    if(strcmp(space,'fourierkw'))
                        if(int_dim == 1)
                            plot_mat = corr_mats{p}(:,1);
                        elseif(int_dim == 2)
                            plot_mat = ifft(corr_mats{p}(2,:)');
                        end
                    else
                        plot_mat = sum(corr_mats{p},3-int_dim)'/size(corr_mats{p},3-int_dim);
                    end
                    figure()
                    plot(1:length(plot_mat),(plot_mat))
                else
                    disp('No valid plot type chosen for correlation plot')
                end
            end
        else
            disp('No valid plot choice')
            disp(argin)
        end
    end
elseif(strcmp(input_type,'data'))
    
end

end

