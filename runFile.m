warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
addpath('wave_matlab/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           User defined variables                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataDirectory = 'eclipseData/';
dataDirectory = '/Users/thomascolligan/Practice_Flight_Data/';
saveDirectory = 'gravityWaveData/';
showPowerSurfaces = true; % Do you want to show the wavelet transform power surfaces?
save = false; % Do you want to save the data? It will save in saveDirectory.
lowerCutOffAltitude = 0; % Altitude where you want to start analysis
upperCutOffAltitude = 40000; % Altitude where you want to end analysis - 
% a value of 40000 will go to the highest point in the profile.
latitude = 46; % Latitude of launch location.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            End of user editing                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: Make quiver arrows the same color as the launch color.
% GIF of wavelet transforms over time with imwrite
% TODO get tropopause values from radiosonde flights
textFiles = fullfile(dataDirectory, "*.txt");
files = dir(textFiles);
f1 = figure;
f2 = figure;
set(0, 'CurrentFigure', f1);
hold on;
set(0, 'CurrentFigure', f2);
hold on;
counter = 0;
first = true;
minLat = Inf;
minLon = Inf;
maxLat = -Inf;
maxLon = -Inf;
% Iterate over files in dataDirectory
for i=1:size(files)
    current = files(i).name;
    fprintf("Current file: %s\n", current);
    current = fullfile(dataDirectory, current);
    if contains(current, 'W4_L1') && ~contains(current, 'bad_data_removed')
        % removed the bad data manually
        continue;
    end
    if ~contains(current, 'W6')
        continue;
    end
    try
        % All analysis logic is in doAnalysis
        [latitudeArray, longitudeArray, altitude, data, ~, ~, ~] = doAnalysis(current, save, saveDirectory, showPowerSurfaces, lowerCutOffAltitude, upperCutOffAltitude);
        % the rest of the code here is plotting and error checking.
        if isempty(data)
            continue;
        end
        % Get latitude of bounding box around all radiosondes
        mLat = min(latitudeArray);
        mLon = min(longitudeArray);
        mxLat = max(latitudeArray);
        mxLon = max(longitudeArray);
        if mxLon > maxLon
            maxLon = mxLon;
        end
        if mxLat > maxLat
            maxLat = mxLat;
        end
        if mLat < minLat
            minLat = mLat;
        end
        if mLon < minLon
            minLon = mLon;
        end
        % plot the latitude, longitude, and altitude in 3d
        set(0, 'CurrentFigure', f1);
        plot3(longitudeArray, latitudeArray, altitude);
        % 3D plot magnitudes of quiver have to be 0 or else the plot is 
        % weirdly scaled. The 3D plot is just for altitude of detection and
        % propagation direction.
        magnitudes = 0*data.axial_ratio + 0.1;
        angle = data.propagation_dir;
        quiver3(data.lon_of_detection, data.lat_of_detection, ...
            data.alt_of_detection_km*1000, magnitudes.*cosd(angle), ...
            magnitudes.*sind(angle), zeros(size(angle), 'like', angle), 0);
        set(0, 'CurrentFigure', f2);
        offset = 10;
        % Plot a vertical line with placeholder, each one offset by 10 from
        % each other on the x-axis.
        placeholder = zeros(size(altitude), 'like', altitude) + (counter*offset);
        plot(placeholder, altitude/1000, 'k'); % plot altitude in km.
        hold on;
        magnitudes = data.axial_ratio;
        % The magnitudes of the red lines are the axial ratios of the
        % gravity waves - axial ratio = intrinsic frequency / coriolis
        % frequency.
        % The for loop below just plots the red lines in the direction of
        % propagation of the gravity wave
        for q=1:size(data.alt_of_detection_km)
            x1 = counter*offset;
            x2 = counter*offset + magnitudes(q)*cosd(angle(q));
            y1 = data.alt_of_detection_km(q);
            y2 = data.alt_of_detection_km(q) + magnitudes(q)*sind(angle(q));
            plot([x1 x2], [y1 y2], 'r');
        end
        scaleFactor = 15;
        if first
            indicesForFilenames = i;
            offsets = counter*offset;
            first = false;
        else
            indicesForFilenames = cat(2, indicesForFilenames, i);
            offsets = cat(2, offsets, counter*offset);
        end
        counter = counter + 1;
    catch e
        if (strcmp(e.identifier, 'MATLAB:table:UnrecognizedVarName'))
            fprintf("Data file %s does not contain a variable needed for analysis. Rerun sounding.\n", files(i).name);
            fprintf("----------------------------\n");
            continue;
        end
        rethrow(e);
    end
        fprintf("----------------------------\n");
end
set(0, 'CurrentFigure', f1);
% xlim([minLon maxLon])
% ylim([minLat maxLat])
xlabel('Longitude (deg)');
ylabel("Latitude (deg)");
zlabel("Altitude (m)");
title("Gravity wave detection altitudes and directions, all summer launches");
tiffPath = 'private/montana_dem.tif';
[mt, R] = geotiffread(tiffPath);
info = geotiffinfo(tiffPath);
mt = double(mt);
[x, y] = pixcenters(info);
h = surf(x, y, mt);
set(h,'LineStyle', 'none')
set(0, 'CurrentFigure', f2);
allFiles = dir(textFiles);
filenames = allFiles(indicesForFilenames)';
xticks(offsets);
xticklabels({filenames.name}); % these curly brackets took 2 hours of my life.
set(gca,'XTickLabelRotation', 45)
ylabel("Altitude (km)")
xlabel("Launch")
title("Propagation direction and altitude of detection of gravity wave");

