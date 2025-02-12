function C = x_y_dist(dist_nums, coordinates)
    X = zeros(size(dist_nums));
    Y = zeros(size(dist_nums));

    for i = 1:size(dist_nums, 1)
        for j = 1:size(dist_nums, 2)
            idx = find(coordinates(:, 1) == dist_nums(i, j), 1);
            if ~isempty(idx)
                Y(i, j) = coordinates(idx, 2);
                X(i, j) = coordinates(idx, 3);
            end
        end
    end

    C = [X Y];    
end
