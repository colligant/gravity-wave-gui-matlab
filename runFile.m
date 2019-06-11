warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
addpath('wave_matlab/');
%data_dir = 'eclipseData/';
data_dir = '/Users/thomascolligan/box/Eclipse 2019/Practice_Flight_Data/Profile';
% GIF of wavelet transforms over time with imwrite
t = fullfile(data_dir, "*.txt");
files = dir(t);
saveDirectory = 'gravityWaveData/';
show = false;
save = true;
lowerCutOffAltitude = 12000;
upperCutOffAltitude = 40000;
latitude = 42.89;
f1 = figure;
f2 = figure;
set(0, 'CurrentFigure', f1);
hold on
minLat = Inf;
minLon = Inf;
maxLat = -Inf;
maxLon = -Inf;
set(0, 'CurrentFigure', f2);
hold on;
counter = 0;
first = true;
for i=1:size(files)
    current = files(i).name;
    fprintf("Current file: %s\n", current);
    current = fullfile(data_dir, current);
    if contains(current, 'W4_L1') && ~contains(current, 'bad_data_removed')
        % removed the bad data manually
        continue;
    end
    try
        [latitudeArray, longitudeArray, altitude, data] = doAnalysis(current, save, saveDirectory, show, lowerCutOffAltitude, upperCutOffAltitude);
        if isempty(data)
            continue;
        end
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
        set(0, 'CurrentFigure', f1);
        mapshow(longitudeArray, latitudeArray);
        magnitudes = data.axial_ratio;
        angle = data.propagation_dir;
        quiver(data.lon_of_detection, data.lat_of_detection, magnitudes.*cosd(angle), magnitudes.*sind(angle))
        set(0, 'CurrentFigure', f2);
        offset = 10;
        placeholder = zeros(size(altitude), 'like', altitude) + (counter*offset);
        x = zeros(size(data.alt_of_detection_km), 'like', data.alt_of_detection_km) + (counter*offset);
        plot(placeholder, altitude/1000, 'k');
        hold on;
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
xlim([minLon maxLon])
ylim([minLat maxLat])
set(0, 'CurrentFigure', f2);
whatIsGoingOn = dir(t);
filenames = whatIsGoingOn(indicesForFilenames)';
xticks(offsets);
xticklabels({filenames.name}); % these curly brackets took 2 hours of my life.
set(gca,'XTickLabelRotation', 45)
ylabel("Altitude (km)")
xlabel("Launch")
title("Propagation direction and altitude of detection of gravity wave");

