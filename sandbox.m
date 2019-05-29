data_dir = '/Users/thomascolligan/Box/matlab-analysis/EclipseData';
t = fullfile(data_dir, "*.txt");
files = dir(t);
saveDirectory = 'gravityWaveData';
show = false;
save = true;
lowerCutOffAltitude = 11000;
for i=1:size(files)
    current = files(i).name;
    fprintf("Current file: %s\n", current);
    current = fullfile(data_dir, current);
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