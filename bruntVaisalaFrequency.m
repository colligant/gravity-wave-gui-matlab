function [N2] = bruntVaisalaFrequency(potentialTemperature, heightSamplingFrequency)
% Calculates BV freq. for a whole sounding.
%

g = ones(size(potentialTemperature), 'like', potentialTemperature) * 9.8;
N2 = (g / potentialTemperature) * gradient(potentialTemperature, heightSamplingFrequency); 
end

