  function [D, P, Q, theta, axialRatio, degreeOfPolarization] = estimateParametersFromWavePacket(uWavePacket, vWavePacket, tempWavePacket)
%UNTITLED2 Summary of this function goes here
%   Estimate gravity wave parameters from a wave packet using the Stokes
%   parameters method.
u = real(uWavePacket);
v = real(vWavePacket);
vWavePacketHilbertTransformed = imag(vWavePacket);
% Stokes parameters from Murphy et al, 2014.
I = mean(u.^2) + mean(v.^2);
D = mean(u.^2) - mean(v.^2);
P = mean(2*u.*v);
Q = real(mean(2*u.*vWavePacketHilbertTransformed));
degreeOfPolarization = sqrt((P^2 + Q^2 + D^2)) / I;
if abs(Q) < 0.05 || abs(P) < 0.05 || degreeOfPolarization < 0.5 || degreeOfPolarization > 1
%if Q < 0.05 || P < 0.05 || degreeOfPolarization < 0.5 || degreeOfPolarization > 1
   theta = 0;
   axialRatio = 0;
   degreeOfPolarization = 0;
   Q = 0;
   return;
else
    theta = 0.5 * atan(P / D);
    % Rotate real parts of u, v data and calculate intrinstic frequency
    rotationMatrix = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    uv = [u; v];
    rotated = rotationMatrix * uv;
    ur = rotated(1, :);
    vr = rotated(2, :);
    axialRatio = abs(mean(ur) / mean(vr));
    uWavevWave = [uWavePacket; vWavePacket]; % zink 2000 eqn 3.17.
    % rotate complex reconstructed wave packets and use to calculate
    % the phase of the coherence function.
    rotatedSpectral = rotationMatrix * uWavevWave;
    ur = rotatedSpectral(1, :);
    numerator = mean(ur.*conj(tempWavePacket));
    denominator = sqrt(mean(abs(ur).^2)*mean(abs(tempWavePacket).^2));
    gamma = numerator / denominator;
    if gamma < 0
        theta = theta + pi;
    end
    if axialRatio < 1
        theta = 0;
        axialRatio = 0;
        degreeOfPolarization = 0;
        Q = 0;
        return;
    end
end
end

