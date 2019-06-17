data_dir = '/Users/thomascolligan/Box/BOREALIS-2019/gravity-wave-analysis/matlab-code/gravityWaveData/';
t = fullfile(data_dir, "*.csv");
files = dir(t);
first = true;
for i=1:size(files)
    current = files(i).name;
    fprintf("Current file: %s\n", current);
    current = fullfile(data_dir, current);
    if contains(current, 'W6')
        continue;
    end
    if first
        data = readGravityWaveData(current); 
        first = false;
    else
        data = [data; readGravityWaveData(current)];
    end
end
ax = polarhistogram(data.propagation_dir, 10);
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');