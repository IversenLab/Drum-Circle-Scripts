function delay_matrix = calculate_delay_matrix(diameter, num_nodes)
    % Define delay factor (1ms/ft) and speed of sound
    delay_factor = diameter / 1.12533;  % Example delay scale factor in ms/ft
    pi = acos(-1.0);

    % Calculate the angular separation between adjacent nodes
    angle_between_nodes = 2.0 * pi / num_nodes;

    % Create a delay matrix to store delays between nodes
    delay_matrix = zeros(num_nodes);

    % Store x and y coordinates for each node
    x = zeros(num_nodes, 1);
    y = zeros(num_nodes, 1);

    % Calculate the positions of all nodes in a circular arrangement
    for i = 1:num_nodes
        angle_pos = (i - 1) * angle_between_nodes;  % Angular position for node i
        x(i) = cos(angle_pos);  % X-coordinate of node i
        y(i) = sin(angle_pos);  % Y-coordinate of node i
    end

    % Calculate the Euclidean distance and delay between each pair of nodes
    for i = 1:num_nodes
        for j = 1:num_nodes
            if i ~= j
                % Calculate the Euclidean distance between nodes i and j
                dx = x(i) - x(j);  % Difference in x-coordinates
                dy = y(i) - y(j);  % Difference in y-coordinates
                euclid_distance = sqrt(dx^2 + dy^2);  % Euclidean distance between nodes i and j

                % Calculate delay based on Euclidean distance and delay factor
                delay_matrix(i, j) = delay_factor % following part for circles* euclid_distance;
            end
        end
    end
end
