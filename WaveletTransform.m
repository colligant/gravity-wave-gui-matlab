classdef WaveletTransform
    %WaveletTransform 
    %   This class acts as a logical grouping of variables for a wavelet
    %   transform. The transform and inverse transform rely on the
    %   same parameters (s0, dj, dt) so keeping the parameters grouped with
    %   the transform will reduce the mental overhead incurred by trying to
    %   pass the same set of parameters to the forward and backward
    %   transform.
    
    properties
        s0;
        dj;
        dt;
        u_wind_component;
        v_wind_component;
        powerSurface;
        uWavelet;
        vWavelet;
        waveletScales;
        fourierPeriod;
        tempWavelet;
        coi;
    end
    
    methods
        function obj = WaveletTransform(u_wind_component, v_wind_component, temperature, dt)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            pad = 1;
            dj = 0.05;
            s0 = 2*dt;
            obj.u_wind_component = u_wind_component;
            obj.v_wind_component = v_wind_component;
            obj.s0 = s0;
            obj.dj = dj;
            obj.dt = dt;
            [obj.uWavelet, ~, obj.waveletScales, ~] = wavelet(u_wind_component, dt, pad, dj, s0); % wave, period, scale, COI
            [obj.vWavelet, ~, ~, obj.coi] = wavelet(v_wind_component, dt, pad, dj, s0);
            [obj.tempWavelet, ~, ~, ~] = wavelet(temperature, dt, pad, dj, s0);
            obj.powerSurface = abs(obj.uWavelet).^2 + abs(obj.vWavelet).^2;
            obj.fourierPeriod = 1.03 * obj.waveletScales;
        end
        
        function [u_wind_reconstructed, v_wind_reconstructed, tempReconstructed, dominantVerticalWavelength] = invertWindowedTransform(obj, windowedWaveletTransform)
            scale_index_1 = windowedWaveletTransform.scale_index_1;
            scale_index_2 = windowedWaveletTransform.scale_index_2;
            alt_index_1 = windowedWaveletTransform.alt_index_1;
            alt_index_2 = windowedWaveletTransform.alt_index_2;
            constant_coef = obj.dj * sqrt(obj.dt) / (0.776*pi^(1/4)); % Magic from T&C.
            windowed_scales = obj.waveletScales(scale_index_1:scale_index_2);
            windowed_u_wavelet = obj.uWavelet(scale_index_1:scale_index_2, alt_index_1:alt_index_2);
            u_wind_reconstructed = constant_coef*sum(windowed_u_wavelet ./ sqrt(windowed_scales)', 1);
            windowed_v_wavelet = obj.vWavelet(scale_index_1:scale_index_2, alt_index_1:alt_index_2);
            windowedTempWavelet = obj.tempWavelet(scale_index_1:scale_index_2, alt_index_1:alt_index_2);
            tempReconstructed = constant_coef*sum(windowedTempWavelet ./ sqrt(windowed_scales)', 1);
            v_wind_reconstructed = constant_coef*sum(windowed_v_wavelet ./ sqrt(windowed_scales)', 1);
            fruig = obj.fourierPeriod(scale_index_2) + obj.fourierPeriod(scale_index_1);
            dominantVerticalWavelength = fruig / 2;
        end
        
        function [uWindInverted, vWindInverted] = invertWaveletTransform(obj)
           constantCoef = obj.dj * sqrt(obj.dt) / (0.776*pi^(1/4)); % Magic from T&C.
           uWindInverted = constantCoef*sum(real(obj.uWavelet) ./ sqrt(obj.waveletScales)', 1); % sum over scales
           vWindInverted = constantCoef*sum(real(obj.vWavelet) ./ sqrt(obj.waveletScales)', 1); 
           
        end
        
        
        function [a, b, c, d] = clipWindowedTransformToValue(obj, localMaxRow, localMaxCol)
            % clipWindowedTransformToValue either clips the windowed
            % transform to ``value'', or clips it to the next inflection
            % point.
            % localMaxRow, localMaxCol are relative to the first indices
            % of the windowedTransform. When we want to find the nearest
            % local minimum, we have to scan the entire power surface. This
            % means we need the "absolute" indices of the local max - i.e.
            % add the initial offset (scale_index_1, row_index_1) to the
            % indices returned by finding the maximum value in the windowed
            % power surface. 
            maxValue = obj.powerSurface(localMaxRow, localMaxCol);    
            thresholdValue = 0.25*maxValue;
            column = obj.powerSurface(:, localMaxCol); % extract whole column
            row = obj.powerSurface(localMaxRow, :); % extract whole row
            
            % need the locations of the first minima on each side of the
            % local max in both row and column. These locations will define
            % the windowed power surface before we clip it to
            % valueToClipTo. s
            
            % Error checking on zero indices returned.
            
            [a, b] = findMinimaClosestToIndex(abs(column-thresholdValue), localMaxRow);
            [c, d] = findMinimaClosestToIndex(abs(row-thresholdValue), localMaxCol);

        end

    end
end

