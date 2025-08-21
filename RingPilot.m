% Create network definition files for RhythmNetwork experiment

% File naming
fpath = '/Users/jonathankirsh/Documents/DrumCircle/PILOT2/5persontest/';  % Directory to save the file
fextension = '.netdef';  % File extension

% Format
clear S  % Clears any existing structure S

S.description = 'Test experiment matlab';
S.experimentDuration = 20;
S.notes = '';

S.networkSize.nodes = 5;  % Number of nodes (e.g., drummers)
S.networkSize.stimulusChannels = 0;
stimulusListeners{1} = [];
stimulusListeners{2} = [];

% Define experiment parts
ep{1}.type = 'network';  % Recording only, network with no connections
ep{1}.description = 'No Connections (record only)';
ep{1}.connections = {[]};
ep{1}.startTime = 0;

%% Define a network
% Define which nodes hear which stimulus
selfFeedback = [1 2 3 4 5];  % Nodes with self-feedback

netType = 'ring';
clear net

% Connectivity matrix
switch netType
    case 'random', % generates a random connectivity matrix with a threshold of 0.5 (i.e., each connection between nodes has a 50% chance of existing)
        thresh = 0.5;
        net = (rand(S.networkSize.nodes) > thresh);

    case 'none',
        network = zeros(S.networkSize.nodes);
        selfFeedback = []; %no connections
        
    case 'self'
        net = zeros(S.networkSize.nodes);
        selfFeedback = 1:S.networkSize.nodes; %only self connections
        
    case 'ring' %ring network
        net = zeros(S.networkSize.nodes);
        for i = 1:S.networkSize.nodes,
            from = i;
            to = mod(from,S.networkSize.nodes)+1;
            net(from,to) = 1;
        end 

    case 'all_to_all'  % Fully connected network
        net = ones(S.networkSize.nodes) - eye(S.networkSize.nodes);  % Set all off-diagonal elements to 1
        selfFeedback = 1:S.networkSize.nodes;  % Include self feedback
    otherwise
        error('Unknown network type')
end

% Loop through different diameters for each experiment part
diameters = [0, 45.0132, 90.0264, 135.0396, 180.0528, 225.066, 270.0792];  % Example diameters for different experiment parts, 45.0132, 90.0264, 135.0396, 180.0528, 225.066, 270.0792
for k = 1:length(diameters)
    diameter = diameters(k);  % Get the current diameter
    
    % Calculate delay matrix based on the current diameter
    delay_matrix = calculate_delay_matrix(diameter, S.networkSize.nodes);
    
    delay_matrix = delay_matrix/2;

for i = 1:S.networkSize.nodes,
   if any(selfFeedback == i),
       net(i,i) = 1;
   else
       net(i,i) = 0;
   end
end
    % Run the connection string
    con = connectionString(net, selfFeedback, stimulusListeners, delay_matrix);
    
    % Define the experiment part for the current diameter
    ep{k+1}.type = 'network';
    ep{k+1}.description = [num2str(diameter) 'ft delay ring'];  % Description with the current diameter
    ep{k+1}.diameter = diameter;
    ep{k+1}.connections = con;
    ep{k+1}.startTime = 10 + (k-1)*60;  % Increment start time for each experiment part
    ep{k+1}.isDelay = true;
    fname_diameter = [fpath 'diameter_' num2str(diameter) fextension];  % File name for the current diameter
    S.experimentParts = ep;  % Update the main structure with current experiment parts
    xml = structToXMLPlist(S);
    saveXMLPlist(fname_diameter, S);
end

% Set the experiment parts in the main structure
S.experimentParts = ep;

% Save the experiment definition to XML

disp(con)
disp(delay_matrix)
