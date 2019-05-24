  function [D, P, theta, axialRatio, degreeOfPolarization, Q, AR] = estimateParametersFromWavePacket(uWavePacket, vWavePacket, tempWavePacket)
%UNTITLED2 Summary of this function goes here
%   Estimate gravity wave parameters from a wave packet using the Stokes
%   parameters method.
u = real(uWavePacket);
v = real(vWavePacket);
vWavePacketHilbertTransformed = imag(vWavePacket);
% Stokes parameters from Murphy et al, 2014.
AR = mean(u) / mean(v);
I = mean(u.^2) + mean(v.^2);
D = mean(u.^2) - mean(v.^2);
P = mean(2*u.*v);
Q = real(mean(2*u.*vWavePacketHilbertTransformed));
degreeOfPolarization = sqrt((P^2 + Q^2 + D^2)) / I;
% Phase difference between u, T from Zink, 2000, eqn 3.17
numerator = mean(uWavePacket.*conj(tempWavePacket));
denominator = sqrt(mean(abs(uWavePacket).^2)*mean(abs(tempWavePacket).^2));
gamma = numerator / denominator;
if (Q < 0.05 && P < 0.05) || degreeOfPolarization < 0.5 || degreeOfPolarization > 1
%if degreeOfPolarization < 0.5 || degreeOfPolarization > 1
   theta = 0;
   axialRatio = 0;
   degreeOfPolarization = 0;
   Q = 0;
   return;
else
    theta = 0.5 * atan(P / D);
    if gamma < 0
        theta = theta + pi;
    end
    rotationMatrix = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    uv = [u; v];
    rotated = rotationMatrix * uv;
    ur = rotated(1, :);
    vr = rotated(2, :);
    arRotated = mean(ur) / mean(vr);
    if abs(arRotated) < 1
        %fprintf("Axial ratio < 1 meaning omega < f. Does not satisfy criteria\n");
        theta = 0;
        axialRatio = 0;
        degreeOfPolarization = 0;
        Q = 0;
        return;
    end
%     plot(ur, vr)
%     hold on;
%     plot(u, v, 'r');
%     uiwait();
    axialRatio = cot(0.5*asin(Q/(degreeOfPolarization*I)));
    fprintf("%f, %f, %f, %f, %f\n", abs(arRotated), abs(1/axialRatio), degreeOfPolarization, Q, P);
end
end

