function doAnalysis(f, save, saveDir, show, lowerCutOffAltitude, upperCutOffAltitude, latitude)
% does g-wave analysis for a radiosonde sounding.
if nargin < 7
    latitude = 42.2127;
end
if nargin < 6
    upperCutOffAltitude = 40000; % flights never reach 40km.
end
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
   % default save in same directory.
   saveFileName = strcat(saveFileName, '_gravity_wave_parameters.csv');
   % if a data directory is provided, use it.
   if saveDirExists
      saveFileName = fullfile(saveDir, saveFileName);
   end
end
data = readRadioSondeData(f);
[maxAlt, mai] = max(data.Alt);
[~, lai] = min(abs(data.Alt(1:mai) - lowerCutOffAltitude)); 
[~, mai] = min(abs(data.Alt(lai:mai) - upperCutOffAltitude));
mai = mai + lai;
if maxAlt < lowerCutOffAltitude && upperCutOffAltitude == 40000
    fprintf("Flight %s did not reach %d m\n", f, lowerCutOffAltitude);
    return;
elseif upperCutOffAltitude ~= 40000 && data.Alt(mai) < (upperCutOffAltitude - 500)
    fprintf("Flight %s did not reach %d m\n", f, upperCutOffAltitude);
    return
elseif mai == lai + 1
    return
end

% Prepare data 
ws = data.Ws(lai:mai);
ws = ws(~isnan(ws));
wd = data.Wd(lai:mai);
wd = wd(~isnan(wd));
pressure = data.P(lai:mai);
pressure = pressure(~isnan(pressure));
alt = data.Alt(lai:mai);
alt = alt(~isnan(alt));
temp = data.T(lai:mai) + 273.15;
temp = temp(~isnan(temp));
time = data.Time(lai:mai);
time = time(~isnan(time));
potentialTemperature = (1000.0^0.286)*temp./(pressure.^0.286); % from Jaxen
u = ws.*cosd(wd);
v = ws.*sind(wd);


% remove background winds with a cubic polynomial.
u = fitAndRemovePolynomial(time, u);
v = fitAndRemovePolynomial(time, v);
temp = fitAndRemovePolynomial(time, temp);

% enforce uniform spatial sampling.
heightSamplingFrequency = 50; % 50m.
u = averageToAltitudeResolution(u, alt, heightSamplingFrequency);
v = averageToAltitudeResolution(v, alt, heightSamplingFrequency);
temp = averageToAltitudeResolution(temp, alt, heightSamplingFrequency);
potentialTemperature = averageToAltitudeResolution(potentialTemperature, alt, heightSamplingFrequency);
alt = averageToAltitudeResolution(alt, alt, heightSamplingFrequency);

% calculate constant values to use in analysis
bvFreqSquared = bruntVaisalaFrequency(potentialTemperature, heightSamplingFrequency); % returns the squared BV frequency.
if size(bvFreqSquared(bvFreqSquared < 0))
    altIndices = bvFreqSquared < 0;
    alts = alt(altIndices);
    disp(alts);
    plot(bvFreqSquared(altIndices))
    uiwait()
end

coriolisFreq = coriolisFrequency(latitude);

% finally, do the wavelet transform.
wt = WaveletTransform(u, v, temp, heightSamplingFrequency);
% get local maxima that (could) correspond to gravity wave packets
% "Peaks were identified as a function of scale and altitude" - MAL2014.
[rows, cols] = find(imregionalmax(wt.powerSurface, 8));
%[rows, cols] = filterLocalMaxOutsideOfCOI(wt.coi, wt.fourierWavelength, rows, cols);
if show
    figure()
    contourf(alt, wt.fourierWavelength, wt.powerSurface);
    [~, titleName, ~] = fileparts(f);
    set(gca,'YScale', 'log')
    titleString = sprintf("%s", titleName);
    title(titleString, 'Interpreter', 'none');
    hold on;
    plot(alt, wt.coi, 'k')
    ylim([wt.s0 Inf])
end
gWaveDetected = false;
for i=1:size(rows)
     % clip the wavelet transform to a box (s1, s2, a1, a2) that 
     % corresponds to where the power surface equals 1/4Smax
     [s1, s2, a1, a2] = wt.clipWindowedTransformToValue(rows(i), cols(i));
     if s1 == 0 || s2 == 0 || a1 == 0 || a2 == 0
         % Too close to an edge.
         continue
     end
     wwt = WindowedWaveletTransform(s1, s2, a1, a2); % helper object to ease passing parameters in to invertWaveletTransform
     % invert the wavelet transform in the windowed transform to get a wave
     % packet.
     % "Wavelet coefficients in the vicinity of the peak were used to
     % reconstruct temperature and wind perturbations associated with each
     % wave around the wave altitude ``z'' and to estimate the vertical
     % wavelength, lambda_z" - MAL2014.
     [ui, vi, tempi, lambda_z] = wt.invertWindowedTransform(wwt);
     % estimateParametersFromWavePacket thresholds wave candidates based on
     % criteria laid out in Murphy et al, 2014.
     [~, ~, Q, theta, axialRatio, degreeOfPolarization] = estimateParametersFromWavePacket(ui, vi, tempi);
     if theta == 0
         % theta = 0 when wave packet does not pass filtering criteria.
        continue;
     end
     if show
        scatter(alt(cols(i)), wt.fourierWavelength(rows(i)), 'ro');
     end
     % all equations below from Murphy et al, 2014:
     % "Radiosonde observations of gravity waves in the lower stratosphere
     %   over Davis, Antartica", table 2.
     intrinsicFreq = coriolisFreq*axialRatio; 
     bvMean = mean(bvFreqSquared(a1:a2)); % Get the mean squared Brunt-Vaisala frequency over the height range of the gravity wave packet
     if ~((sqrt(bvMean) > intrinsicFreq) && (intrinsicFreq > coriolisFreq))
         % if the intrinsic frequency is greater than the buoyancy
         % frequency or if it's less than the coriolis frequency, the
         % gravity wave is not physical, so skip it.
         continue
     end
     gWaveDetected = true;
     m = 2*pi / lambda_z; % vertical wavenumber (1 / meters)
     k_h = sqrt(((coriolisFreq^2*m^2)/(bvMean))*(intrinsicFreq^2/coriolisFreq^2 - 1)); % horizontal wavenumber (1 / meters)
     intrinsicVerticalGroupVel = -(1 / (intrinsicFreq*m))*(intrinsicFreq^2 - coriolisFreq^2); % m/s
     zonalWaveNumber = k_h*sin(theta);% 1/m
     meridionalWaveNumber = k_h*cos(theta); % 1 / m
     intrinsicVerticalPhaseSpeed = intrinsicFreq / m; % m/s
     intrinsicHorizPhaseSpeed = intrinsicFreq / k_h; % m/s
     intrinsicZonalGroupVel = zonalWaveNumber * bvMean / (intrinsicFreq * m^2); % m/s
     intrinsicMeridionalGroupVel = meridionalWaveNumber * bvMean / (intrinsicFreq * m^2); % m/s
     intrinsicHorizGroupVel = sqrt(intrinsicZonalGroupVel^2 + intrinsicMeridionalGroupVel^2); % m/s
     lambda_h = 2*pi / k_h; % horizontal wavelength (m)
     altitudeOfDetection = mean(alt(a1:a2)); % mean of the altitude range as the central altitude of the gravity wave.
     if save
         % finally, save the data.
         header = {'altOfDetection_km', 'vert_wavelength_km', 'horiz_wavelength_km', 'propagation_dir', 'axial_ratio', 'int_vert_group vel_ms)', 'int_horiz_group_vel_ms', 'int_vert_phase_spd_ms', 'int_horiz_phase_spd_ms)', 'degreeofpolarization', 'stokes_param_Q'};
         data = [altitudeOfDetection/1000 lambda_z/1000 lambda_h/1000 rad2deg(theta), axialRatio, intrinsicVerticalGroupVel, intrinsicHorizGroupVel, intrinsicVerticalPhaseSpeed, intrinsicHorizPhaseSpeed, degreeOfPolarization, Q];
         if ~isfile(saveFileName)
             % check if the file exists - if it doesn't, write a header to
             % the file to ease further analysis.
             writecell([header; num2cell(data)], saveFileName);
         else
             % if the file does exist, append to it.
             % There is no overwrite functionality, so it's possible to run
             % the analysis multiple times and write duplicates to a file.
             % This must be taken care of in later analysis.
             dlmwrite(saveFileName, data, 'delimiter', ',', '-append');
         end
     end
     fprintf("Alt: %f, L_z: %f, L_h: %f, Theta: %f, w/f: %f, period (hours): %f, Vert. group vel: %f, Horiz. group vel: %f, Vert. phase spd: %f, Horiz. phase spd: %f\n", altitudeOfDetection/1000, lambda_z/1000, lambda_h/1000, rad2deg(theta), axialRatio, (2*pi/intrinsicFreq)/3600, intrinsicVerticalGroupVel, intrinsicHorizGroupVel, intrinsicVerticalPhaseSpeed, intrinsicHorizPhaseSpeed);
end
if ~gWaveDetected
    fprintf("No gravity waves detected.\n");
end
end

