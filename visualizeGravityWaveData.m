data_dir ='/Users/thomascolligan/Box/BOREALIS-2019/gravity-wave-analysis/matlab-code/gravityWaveData/';
t = fullfile(data_dir, "*.csv");
files = dir(t);
first = true;
for i=1:size(files)
    current = files(i).name;
    fprintf("Current file: %s\n", current);
    current = fullfile(data_dir, current);
    if first
        data = readGravityWaveData(current); 
        first = false;
    else
        data = [data; readGravityWaveData(current)];
    end
end
figure
plot3(data.lon_of_detection, data.lat_of_detection, data.alt_of_detection_km, 'ro');
hold on;
ar = data.axial_ratio;
quiver3(data.lon_of_detection, data.lat_of_detection, data.alt_of_detection_km, ar.*cosd(data.propagation_dir), ar.*sind(data.propagation_dir), 0*data.alt_of_detection_km)
grid on;
figure
plot(data.lon_of_detection, data.lat_of_detection, 'ro');
hold on;
ar = data.axial_ratio;
quiver(data.lon_of_detection, data.lat_of_detection, ar.*cosd(data.propagation_dir), ar.*sind(data.propagation_dir))
grid on;







