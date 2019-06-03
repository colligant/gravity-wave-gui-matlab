warning('off','all')

data_dir ='/Users/thomascolligan/box/Eclipse 2019/Practice_Flight_Data';
%data_dir = '/Users/thomascolligan/Box/matlab-analysis/EclipseData';
t = fullfile(data_dir, "*.txt");
files = dir(t);
saveDirectory = 'gravityWaveData';
show = false;
save = true;
lowerCutOffAltitude = 12000;
for i=1:size(files)
    current = files(i).name;
    fprintf("Current file: %s\n", current);
    current = fullfile(data_dir, current);
    if contains(current, 'W4_L1') && ~contains(current, 'bad_data_removed')
        continue;
    end
    try
        doAnalysis(current, save, saveDirectory, show, lowerCutOffAltitude);
    catch e
        if (strcmp(e.identifier, 'MATLAB:table:UnrecognizedVarName'))
            fprintf("----------------------------\n");
            fprintf("Exception occured when reading file %s with message:\n%s\n", current, e.message);
            fprintf("----------------------------\n");
            continue;
        end
        rethrow(e);
    end
        fprintf("----------------------------\n");
end