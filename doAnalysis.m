function doAnalysis(f, save, saveDir, show, lowerCutOffAltitude)
% does g-wave analysis for a radiosonde sounding.]
if nargin < 5
    lowerCutOffAltitude = 12000;
end
if nargin < 4
    show = false;
end
if nargin < 3
    saveDirExists = false;
else
    saveDirExists = true;
end
if nargin < 2
    save = false;
end
if save
   [~, saveFileName, ~] = fileparts(f);
   saveFileName = strcat(saveFileName, '_gravity_wave_parameters.csv');
   if saveDirExists
      saveFileName = fullfile(saveDir, saveFileName);
   end
end
data = readRadioSondeData(f);
latitudeFortLaramie = 42.2127;
[maxAlt, mai] = max(data.Alt);
[~, lai] = min(abs(data.Alt(1:mai) - lowerCutOffAltitude)); 
if maxAlt < lowerCutOffAltitude
    fprintf("Flight %s did not reach %d m\n", f, lowerCutOffAltitude);
    return;
end

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

% enforce uniform spatial sampling.
heightSamplingFrequency = 50;
u = averageToAltitudeResolution(u, alt, heightSamplingFrequency);
v = averageToAltitudeResolution(v, alt, heightSamplingFrequency);
temp = averageToAltitudeResolution(temp, alt, heightSamplingFrequency);
potentialTemperature = averageToAltitudeResolution(potentialTemperature, alt, heightSamplingFrequency);
alt = averageToAltitudeResolution(alt, alt, heightSamplingFrequency);

bvFreq = bruntVaisalaFrequency(potentialTemperature, heightSamplingFrequency); % this is N^2, not N.
coriolisFreq = coriolisFrequency(latitudeFortLaramie);
wt = WaveletTransform(u, v, temp, heightSamplingFrequency);
[rows, cols] = find(imregionalmax(wt.powerSurface));
if show
    contourf(alt, wt.fourierPeriod, wt.powerSurface);
    set(gca,'YScale', 'log')
    title(f);
    hold on;
    plot(alt, wt.coi, 'k')
    ylim([wt.s0 Inf])
end


for i=1:size(rows)
     [s1, s2, a1, a2] = wt.clipWindowedTransformToValue(rows(i), cols(i));
     if s1 == 0 || s2 == 0 || a1 == 0 || a2 == 0
         continue
     end
     wwt = WindowedWaveletTransform(s1, s2, a1, a2);
     [ui, vi, tempi, lambda_z] = wt.invertWindowedTransform(wwt);
     [~, ~, Q, theta, axialRatio, degreeOfPolarization] = estimateParametersFromWavePacket(ui, vi, tempi);
     if theta == 0
        continue;
     end
     if show
        scatter(alt(cols(i)), wt.fourierPeriod(rows(i)), 'ro');
     end
     intrinsicFreq = coriolisFreq*axialRatio;
     bvMean = mean(bvFreq(a1:a2));
     m = 2*pi / lambda_z;
     k_h = sqrt(((intrinsicFreq^2 - coriolisFreq^2)*m^2) / (bvMean));
     intrinsicVerticalGroupVel = -(1 / (intrinsicFreq*m))*(intrinsicFreq^2 - coriolisFreq^2);
     zonalWaveNumber = k_h*sin(theta);
     meridionalWaveNumber = k_h*cos(theta);
     intrinsicVerticalPhaseSpeed = intrinsicFreq / m;
     intrinsicHorizPhaseSpeed = intrinsicFreq / k_h;
     intrinsicZonalGroupVel = zonalWaveNumber * bvMean / (intrinsicFreq * m^2);
     intrinsicMeridionalGroupVel = meridionalWaveNumber * bvMean / (intrinsicFreq * m^2);
     intrinsicHorizGroupVel = sqrt(intrinsicZonalGroupVel^2 + intrinsicMeridionalGroupVel^2);
     lambda_h = 2*pi / k_h;
     altitudeOfDetection = mean(alt(a1:a2));
     if ~(sqrt(bvMean) > intrinsicFreq) && (intrinsicFreq > coriolisFreq)
         continue
     end
     if save
         header = {'altOfDetection (km)', 'vert. wavelength (km)', 'horiz. wavelength (km)', 'propagation dir (deg ccw from E)', 'axial ratio (unitless)', 'int. vert. group vel (m/s)', 'int. horiz group vel (m/s)', 'int. vert. phase spd (m/s)', 'int. horiz phase spd (m/s)', 'degreeofpolarization', 'stokes param Q'};
         data = [altitudeOfDetection/1000 lambda_z/1000 lambda_h/1000 rad2deg(theta), axialRatio, intrinsicVerticalGroupVel, intrinsicHorizGroupVel, intrinsicVerticalPhaseSpeed, intrinsicHorizPhaseSpeed, degreeOfPolarization, Q];
         if ~isfile(saveFileName)
             writecell([header; num2cell(data)], saveFileName);
         else
             dlmwrite(saveFileName, data, 'delimiter', ',', '-append');
         end
     end
     fprintf("Alt: %f, L_z: %f, L_h: %f, Theta: %f, w/f: %f, Vert. group vel: %f, Horiz. group vel: %f, Vert. phase spd: %f, Horiz. phase spd: %f\n", altitudeOfDetection/1000, lambda_z/1000, lambda_h/1000, rad2deg(theta), axialRatio, intrinsicVerticalGroupVel, intrinsicHorizGroupVel, intrinsicVerticalPhaseSpeed, intrinsicHorizPhaseSpeed);
end
if show
    uiwait()
end
end

