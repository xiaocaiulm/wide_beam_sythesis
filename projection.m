function x_projected = projection(x, s)
    % Ensure x and s are column vectors
    s = s(:);

    % Get the size of x
    [rows, cols] = size(x);

    % Initialize the output matrix
    x_projected = zeros(size(x));

    % Calculate the closest point in s for each element in x
    for i = 1:rows
        for j = 1:cols
            % Calculate the Euclidean distance from the current element of x to each element of s
            distances = abs(x(i, j) - s.');

            % Find the index of the element in s that is closest to the current element of x
            [~, closest_point_idx] = min(distances);

            % Retrieve the closest point
            x_projected(i, j) = s(closest_point_idx);
        end
    end
end

