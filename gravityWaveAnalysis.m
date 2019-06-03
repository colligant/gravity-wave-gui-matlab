data_dir ='/Users/thomascolligan/Box/BOREALIS-2019/gravity-wave-analysis/matlab-code/gravityWaveData';
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
%data.Properties.VariableNames
%histogram(data.axial_ratio, 5)
%size(data)
%mean(data.horiz_wavelength_km)
figure
subplot(2, 2, 2)
polarhistogram(data.propagation_dir, 5);
thetaticks([])
title("Propagation directions, n=23")
subplot(2, 2, 1)
histogram(data.int_horiz_group_vel_ms, 5)
title("Horiz. group velocity, n=23")
xlabel("group velocity, (m/s)")
subplot(2, 2, [3, 4])
histogram(data.int_horiz_phase_spd_ms_, 5)
title("Horiz. phase speed, n=23")
xlabel("phase speed, (m/s)")
saveas(gcf, "/Users/thomascolligan/Desktop/analysis.png")
