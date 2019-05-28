function [N2] = bruntVaisalaFrequency(potentialTemperature, heightSamplingFrequency)
% Calculates BV frequency for a whole sounding, using
% N^2 = g/theta * dtheta / dz, where theta is potential temperature.
g = ones(size(potentialTemperature), 'like', potentialTemperature) * 9.8;
N2 = (g / potentialTemperature) * gradient(potentialTemperature, heightSamplingFrequency); 
end

