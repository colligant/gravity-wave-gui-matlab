function doAnalysis(f)
% does g-wave analysis for a radiosonde sounding.
data = readRadioSondeData(f);
warning('off','all')
latitudeFortLaramie = 42.2127;
[maxAlt, mai] = max(data.Alt);
[~, lai] = min(abs(data.Alt - 12000));

ws = data.Ws(lai:mai);
wd = data.Wd(lai:mai);
pressure = data.P(lai:mai);
alt = data.Alt(lai:mai);
temp = data.T(lai:mai) + 273.15;
time = data.Time(lai:mai);
potentialTemperature = (1000.0^0.286)*temp./(pressure.^0.286); %kelvin
u = ws.*cosd(wd);
v = ws.*sind(wd);

u = fitAndRemovePolynomial(time, u);
v = fitAndRemovePolynomial(time, v);
temp = fitAndRemovePolynomial(time, temp);

heightSamplingFrequency = 50;
u = averageToAltitudeResolution(u, alt, heightSamplingFrequency);
v = averageToAltitudeResolution(v, alt, heightSamplingFrequency);
temp = averageToAltitudeResolution(temp, alt, heightSamplingFrequency);
potentialTemperature = averageToAltitudeResolution(potentialTemperature, alt, heightSamplingFrequency);
alt = averageToAltitudeResolution(alt, alt, heightSamplingFrequency);

bvFreq = bruntVaisalaFrequency(potentialTemperature, heightSamplingFrequency); % this is N^2, not N.
coriolisFreq = coriolisFrequency(latitudeFortLaramie);
wt = WaveletTransform(u, v, temp, heightSamplingFrequency);


contourf(alt, wt.fourierPeriod, wt.powerSurface);
set(gca,'YScale', 'log')
title(f);
hold on;
plot(alt, wt.coi, 'k')
[rows, cols] = find(imregionalmax(wt.powerSurface));
ylim([wt.s0 Inf])


for i=1:size(rows)
     [s1, s2, a1, a2] = wt.clipWindowedTransformToValue(rows(i), cols(i));
     if s1 == 0 || s2 == 0 || a1 == 0 || a2 == 0
         continue
     end
     wwt = WindowedWaveletTransform(s1, s2, a1, a2);
     [ui, vi, tempi, lambda_z] = wt.invertWindowedTransform(wwt);
     [D, P, Q, theta, axialRatio, d] = estimateParametersFromWavePacket(ui, vi, tempi);
     if theta == 0
        continue;
     end
     scatter(alt(cols(i)), wt.fourierPeriod(rows(i)), 'ro');
     intrinsicFreq = coriolisFreq*axialRatio;
     bvMean = mean(bvFreq(a1:a2));
     m = 2*pi / lambda_z;
     k_h = sqrt(((intrinsicFreq^2 - coriolisFreq^2)*m^2) / (bvMean^2));
     intrinsicVerticalGroupVel = -(1 / (intrinsicFreq*m))*(intrinsicFreq^2 - coriolisFreq^2);
     zonalWaveNumber = k_h*sin(theta);
     meridionalWaveNumber = k_h*cos(theta);
     intrinsicVerticalPhaseSpeed = intrinsicFreq / m;
     intrinsicHorizPhaseSpeed = intrinsicFreq / k_h;
     intrinsicZonalGroupVel = zonalWaveNumber * bvMean / (intrinsicFreq * m^2);
     intrinsicMeridionalGroupVel = meridionalWaveNumber * bvMean / (intrinsicFreq * m^2);
     intrinsicHorizGroupVel = sqrt(intrinsicZonalGroupVel^2 + intrinsicMeridionalGroupVel^2);
     lambda_h = 2*pi / k_h;
     %fprintf("%f, %f, %f, %f, %f, %f, %f, %f\n", lambda_z, lambda_h, rad2deg(theta), axialRatio, intrinsicVerticalGroupVel, intrinsicHorizGroupVel, intrinsicVerticalPhaseSpeed, intrinsicHorizPhaseSpeed);
     fprintf("%f, %f\n", lambda_h, lambda_z);
end
uiwait()
end

