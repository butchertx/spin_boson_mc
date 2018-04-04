function [ combined ] = cat_measures( full_meas )
%CAT_MEASURES combine full measurements for multiple simulations into one
%   assume full_meas is a matrix made from the dump files listed in
%   multiple different simulation directories.  full_meas has dimensions
%   [measure_type,measurements,num_alphas,num_sim_dirs], so this just takes
%   the fourth dimension and concatenates each to return a 3-d matrix

combined = zeros(size(full_meas,1),size(full_meas,2)*size(full_meas,4),size(full_meas,3));
for n = 1:size(full_meas,4)
    combined(:,(n-1)*size(full_meas,2) + 1 : n*size(full_meas,2),:) = full_meas(:,:,:,n);
end

end

