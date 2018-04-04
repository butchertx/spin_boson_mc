function [ autocorr ] = autocorrelation( vals )
%AUTOCORRELLATION calculate the autocorrelation function of a list of
%values

autocorr = fft(vals);
autocorr = abs(autocorr).^2;
autocorr = ifft(autocorr)/length(vals);

end

