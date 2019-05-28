%f = 'EDBackup/'
f = '/Users/thomascolligan/box/Eclipse 2019/Practice_Flight_Data/';
t = fullfile(f, "*.txt");
files = dir(t);

for i=1:size(files)
    current = files(i).name;
    fprintf("Current file: %s\n", current);
    current = fullfile(f, current);
    try
        doAnalysis(current);
    catch
        continue;
    end
        fprintf("----------------------------\n");
end















