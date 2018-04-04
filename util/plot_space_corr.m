function [full_corr, space_corr] = plot_space_corr(filename, varargin)
%Plot C(r,omega=0) for the given file (.csv) with dimensions listed in the first row
%Order of options
%'stagger'
%'logy','loglog', or 'logx'

dat = csvread(filename);%this needs to be altered if the size is larger than 100000 elements
lx = dat(1,1);
ly = dat(1,2);
full_corr = split_corr(dat(2,:),lx,ly);
space_corr = integrate_corr(dat(2,:),lx,ly);
plot_corr = space_corr;

if(nargin>=2)
    if(strcmp(varargin{1},'stagger'))
        indices = -(-1).^(1:lx)';
        plot_corr = indices.*space_corr;
        if(length(varargin)>2)
            varargin = varargin{2:end};
        elseif (length(varargin)==2)
            varargin{1} = varargin{2};
        end
    end

    if(strcmp(varargin{1},'logy'))
        plot_corr = log(plot_corr);
        if(length(varargin)>2)
            varargin = varargin{2:end};
        elseif (length(varargin)==2)
            varargin{1} = varargin{2};
        end
    end
end

figure()
plot(plot_corr);

end

