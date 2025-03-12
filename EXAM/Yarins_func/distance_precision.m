function precision_mm = distance_precision(distance_meter, mm_precision, ppm_precision)
    
    % Convert distance to kilometers
    distance_km = distance_meter / 1000;
    
    % Calculate precision in km
    precision_mm = mm_precision + (ppm_precision* distance_km);
end


