%with pairs
% File naming
fpath = '/Users/jonathankirsh/Documents/DrumCircle/PILOT2/pairs/';  % Directory to save the file
fextension = '.netdef';  % File extension

% Format
clear S  % Clears any existing structure S

S.description = 'Test experiment matlab';
S.experimentDuration = 20;
S.notes = '';

S.networkSize.nodes = 6;  % Number of nodes (e.g., drummers)
S.networkSize.stimulusChannels = 0;
stimulusListeners{1} = [];
stimulusListeners{2} = [];

% Define experiment parts
ep{1}.type = 'network';  % Recording only, network with no connections
ep{1}.description = 'No Connections (record only)';
ep{1}.connections = {[]};
ep{1}.startTime = 0;

%% Define a network with 3 pairs of nodes and self-feedback
netType = 'pairs_with_self_feedback';  % Define new network type for pairs and self-feedback
clear net

switch netType
    case 'pairs_with_self_feedback'  % Manually define 3 pairs of nodes with self-feedback
        % Define pairs of nodes (each row defines a pair of nodes)
        pairs = [1 2; 3 4; 5 6];  % Example pairs: (1,2), (3,4), (5,6)
        
        % Initialize connectivity matrix
        net = zeros(S.networkSize.nodes);
        
        % Set pairs to be connected
        for i = 1:size(pairs, 1)
            net(pairs(i, 1), pairs(i, 2)) = 1;  % Connect first node in pair to second
            net(pairs(i, 2), pairs(i, 1)) = 1;  % Connect second node in pair to first
        end
        
        % Add self-feedback (self-connections)
        for i = 1:S.networkSize.nodes
            net(i, i) = 1;  % Add self-feedback connection for each node
        end
        
        selfFeedback = 1:S.networkSize.nodes;  % Self-feedback for all nodes

    otherwise
        error('Unknown network type')
end

% Loop through different diameters for each experiment part
diameters = [0, 45.0132, 90.0264, 135.0396, 180.0528, 225.066, 270.0792];  % Example diameters for different experiment parts
for k = 1:length(diameters)
    diameter = diameters(k);  % Get the current diameter
    
    % Calculate delay matrix based on the current diameter
    delay_matrix = calculate_delay_matrix(diameter, S.networkSize.nodes);
    
    % Adjust delay matrix by dividing by 2 to correct for scaling issues
    %delay_matrix = delay_matrix / 2;
    
    % Run the connection string
    con = connectionString(net, selfFeedback, stimulusListeners, delay_matrix);
    
    % Define the experiment part for the current diameter
    ep{2}.type = 'network';
    ep{2}.description = [num2str(diameter) 'ft delay pairs'];  % Description with the current diameter
    ep{2}.diameter = diameter;
    ep{2}.connections = con;
    ep{2}.startTime = 10;  % Increment start time for each experiment part
    ep{2}.isDelay = true;

    % Save each experiment part to a separate XML file
    fname_diameter = [fpath 'diameterpair_' num2str(diameter) fextension];  % File name for the current diameter
    S.experimentParts = ep;  % Update the main structure with current experiment parts
    xml = structToXMLPlist(S);
    saveXMLPlist(fname_diameter, S);
end

% Set the experiment parts in the main structure
S.experimentParts = ep;

disp(con)
disp(delay_matrix)
