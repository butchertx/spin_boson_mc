function [ resampled_measurements ] = swendsen_resample_point( alpha, S1, a_resample, measurements )
%SWENDSEN_RESAMPLE_POINT resample a set of measurements
%   action is of the form S = S_0 + alpha*S1
%   a_resample is a list of the desired alphas to resample at
%   measurements is a list of measurements of a given order parameter
%   corresponding to the values listed in S1

%   Assume length(S1) = length(measurements)

%%create a histogram h(S1)
num_bins = 100;
histogram = zeros(num_bins, 1);
bin_length = (max(S1) - min(S1))/(num_bins(1) - 1);
for i = 1:length(S1)
    histogram(floor((S1(i) - min(S1))/bin_length) + 1) = 1 + histogram(floor((S1(i) - min(S1))/bin_length) + 1);
end
histogram = histogram/length(measurements);

resampled_measurements = zeros(length(a_resample), 1);
for i = 1:length(a_resample)
    temp_meas = 0;
    part_func = 0;
    %approximate the partition function for the new value of alpha
%     for j = 1:length(measurements)
%         temp_prob = histogram(floor((S(j) - min(S))/bin_length) + 1)/length(measurements);
%         part_func = part_func + temp_prob*exp(-(a_resample(i) - alpha)*S1(j));
%     end
    for j = 1:num_bins
        part_func = part_func + histogram(j)*exp(-(a_resample(i) - alpha)*(bin_length*(j - 0.5) + min(S1)));
    end
    
    %calculate the new histogram
    new_hist = zeros(num_bins,1);
    for j = 1:num_bins
        new_hist(j) = histogram(j)*exp(-(a_resample(i) - alpha)*(bin_length*(j - 0.5) + min(S1)))/part_func;
    end
    norm = trapz(new_hist);
    new_hist = new_hist/norm;
    
    %average the measurements according to this partition function
    for j = 1:length(measurements)
        temp_prob = new_hist(floor((S1(j) - min(S1))/bin_length) + 1)/histogram(floor((S1(j) - min(S1))/bin_length) + 1)/length(measurements);
        temp_meas = temp_meas + measurements(j)*temp_prob;
    end
        
    resampled_measurements(i) =  temp_meas;
end

end