f = 'EDBackup/08212017S5_SouthWY_BLaunch_Murdock181828.txt';
data = readRadioSondeData(f);
warning('off','all')

[maxAlt, mai] = max(data.Alt);
[~, lai] = min(abs(data.Alt - 12000));

ws = data.Ws(lai:mai);
wd = data.Wd(lai:mai);
alt = data.Alt(lai:mai);
temp = data.T(lai:mai) + 273.15;
time = data.Time(lai:mai);
u = ws.*cosd(wd);
v = ws.*sind(wd);
u = fitAndRemovePolynomial(time, u);
v = fitAndRemovePolynomial(time, v);
temp = fitAndRemovePolynomial(time, temp);
heightSamplingFrequency = 50;
u = averageToAltitudeResolution(u, alt, heightSamplingFrequency);
v = averageToAltitudeResolution(v, alt, heightSamplingFrequency);
temp = averageToAltitudeResolution(temp, alt, heightSamplingFrequency);
alt = averageToAltitudeResolution(alt, alt, heightSamplingFrequency);
wt = WaveletTransform(u, v, temp, heightSamplingFrequency);
[rows, cols] = find(imregionalmax(wt.powerSurface));

% contourf(alt, wt.fourierPeriod, wt.powerSurface);
% set(gca, 'Yscale', 'log');
% hold on;
% scatter(alt(cols), wt.fourierPeriod(rows));


for i=1:size(rows)
     [s1, s2, a1, a2] = wt.clipWindowedTransformToValue(rows(i), cols(i));
     if s1 == 0 || s2 == 0 || a1 == 0 || a2 == 0
         continue
     end
     wwt = WindowedWaveletTransform(s1, s2, a1, a2);
     [ui, vi, tempi, m] = wt.invertWindowedTransform(wwt);
     lambda_z = m;
     fprintf("%f\n", lambda_z);
     [D, P, t, ar, dp, Q, AR] = estimateParametersFromWavePacket(ui, vi, tempi);
     
end