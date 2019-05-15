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
        power_surface;
        u_wavelet;
        v_wavelet;
        waveletScales;
        coi;
    end
    
    methods
        function obj = WaveletTransform(u_wind_component, v_wind_component, dt)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            pad = 1;
            dj = 0.01;
            s0 = 2*dt;
            obj.u_wind_component = u_wind_component;
            obj.v_wind_component = v_wind_component;
            obj.s0 = s0;
            obj.dj = dj;
            obj.dt = dt;
            
            [obj.u_wavelet, ~, obj.waveletScales, ~] = wavelet(u_wind_component, dt, pad, dj, s0); % wave, period, scale, COI
            [obj.v_wavelet, ~, ~, obj.coi] = wavelet(v_wind_component, dt, pad, dj, s0);
            obj.power_surface = abs(obj.u_wavelet).^2 + abs(obj.v_wavelet).^2;
        end
        
        function [u_wind_inverted, v_wind_inverted, v_wind_hilbert_transformed] = invertWindowedTransform(obj, windowedWaveletTransform)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            scale_index_1 = windowedWaveletTransform.scale_index_1;
            scale_index_2 = windowedWaveletTransform.scale_index_2;
            alt_index_1 = windowedWaveletTransform.alt_index_1;
            alt_index_2 = windowedWaveletTransform.alt_index_2;
            constant_coef = obj.dj * sqrt(obj.dt) / (0.776*pi^(-1/4)); % Magic from T&C.
            windowed_scales = obj.waveletScales(scale_index_1:scale_index_2);
            windowed_u_wavelet = obj.u_wavelet(scale_index_1:scale_index_2, alt_index_1:alt_index_2);
            u_wind_reconstructed = constant_coef*sum(windowed_u_wavelet ./ sqrt(windowed_scales)');
            windowed_v_wavelet = obj.v_wavelet(scale_index_1:scale_index_2, alt_index_1:alt_index_2);
            v_wind_reconstructed = constant_coef*sum(windowed_v_wavelet ./ sqrt(windowed_scales)');
            v_wind_inverted = real(v_wind_reconstructed);
            u_wind_inverted = real(u_wind_reconstructed);
            v_wind_hilbert_transformed = complex(v_wind_reconstructed);
            
        end
        
        function [a, b, c, d] = clipWindowedTransformToValue(obj, windowedWaveletTransform, localMaxRow, localMaxCol)
            % clipWindowedTransformToValue either clips the windowed
            % transform to ``value'', or clips it to the next inflection
            % point.
            scale_index_1 = windowedWaveletTransform.scale_index_1;
            scale_index_2 = windowedWaveletTransform.scale_index_2;
            alt_index_1 = windowedWaveletTransform.alt_index_1;
            alt_index_2 = windowedWaveletTransform.alt_index_2;
            windowedTransform = obj.power_surface(scale_index_1:scale_index_2, alt_index_1:alt_index_2);
            if nargin < 3
               % get local maxima in the candidate surface passed in by
               % user rectangle
               maxOfCols = max(windowedTransform);
               [maxValue, ~] = max(maxOfCols);
               [localMaxRow, localMaxCol] = find(windowedTransform == maxValue);
               localMaxRow = localMaxRow + scale_index_1;
               localMaxCol = localMaxCol + alt_index_1;
            end
            % localMaxRow, localMaxCol are relative to the first indices
            % of the windowedTransform. When we want to find the nearest
            % local minimum, we have to scan the entire power surface. This
            % means we need the "absolute" indices of the local max - i.e.
            % add the initial offset (scale_index_1, row_index_1) to the
            % indices returned by finding the maximum value in the windowed
            % power surface. 
            maxValue = obj.power_surface(localMaxRow, localMaxCol);    
            thresholdValue = 0.25*maxValue;
            column = obj.power_surface(:, localMaxCol); % extract whole column
            row = obj.power_surface(localMaxRow, :); % extract whole row
            
            % need the locations of the first minima on each side of the
            % local max in both row and column. These locations will define
            % the windowed power surface before we clip it to
            % valueToClipTo. 
            
            % Error checking on zero indices returned.
            
            [a, b] = findMinimaClosestToIndex(abs(column-thresholdValue), localMaxRow);
            [c, d] = findMinimaClosestToIndex(abs(row-thresholdValue), localMaxCol);

        end
        
        
        
        
    end
end

