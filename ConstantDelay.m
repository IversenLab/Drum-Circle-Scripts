%Attempt to make circle with delays
% -	Simulate a physical drum circle of different diameters – different delays based on proximity of tappers
% Need a program that you could give it a list of participants (e.g. 1-6 to define size), 
% and then space them equally around a circle. For a diameter of 1, calculate the delay factors. 
% Then have a diameter factor (1ms/foot), then write out a list of those connection strings, end up saving some kind of a text file

% Real circle is about 6ft in diameter*

% create network definition files for RhythmNetwork experiment
%
%   saves a .netdef file that defines an experiment
%
%   an example file is appended to end
%
%   constraints: must have some network defined before stimuli are defined
%
%   JRI 3/16/05
% TO-DO - BUILD OUT FULL EXPERIMENT, try 3 different delay settings*
%file naming
fpath = '/Users/jonathankirsh/Documents/DrumCircle/'; %specifies the directory where the file will be saved
fextension = '.netdef'; % file extension
fname = [fpath 'solnaug11exp100ft' fextension]; %cobines the path and file name to define the full path for the file

%format
clear S % clears any existing structure S and sets up basic info for the experiment

S.description = 'Test experiment matlab';
S.experimentDuration = 20;
S.notes = '';


S.networkSize.nodes = 6;
S.networkSize.stimulusChannels = 0;
stimulusListeners{1} = [];
stimulusListeners{2} = [];

%define experiment parts

ep{1}.type = 'network'; %recording only, network with no connections
ep{1}.description = 'No Connections (record only)';
ep{1}.connections = {[]};
ep{1}.startTime = 0;




%% define a network
%define which nodes hear which stimulus

%set self connections
%nodes with selfFeedback
selfFeedback = [1 2 3 4 5 6]; 

%get a connectivity matrix
netType = 'all_to_all';
clear net



%%
% connectivity matrix continued

switch netType,
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
        selfFeedback = 1:S.networkSize.nodes;  % include self feedback

    otherwise
        error('unknown network type')
end

%% make a circle
% Define parameters
diameter = 100;  % diameter in feet
delay_factor = diameter / 1.12533;  % Delay scale factor in ms/ft

% Speed of sound is 1125.33 feet per second (1/1125.33 gives seconds per foot)
pi = acos(-1.0);
% angle = pi / 2.0;  % Start angle is 12 o'clock

% Calculate the angular separation between adjacent nodes
angle_between_nodes = 2.0 * pi / S.networkSize.nodes;  % 360 degrees over the number of nodes

% Create a delay matrix to store delays between nodes
delay_matrix = zeros(S.networkSize.nodes);

% Store x and y coordinates for each node
x = zeros(S.networkSize.nodes, 1);
y = zeros(S.networkSize.nodes, 1);

% Calculate the positions of all nodes in a circular arrangement
for i = 1:S.networkSize.nodes
    angle_pos = (i - 1) * angle_between_nodes;  % Angular position for node i
    x(i) = cos(angle_pos);  % X-coordinate of node i
    y(i) = sin(angle_pos);  % Y-coordinate of node i
end

% Calculate the Euclidean distance and delay between each pair of nodes
for i = 1:S.networkSize.nodes
    for j = 1:S.networkSize.nodes
        if i ~= j
            % Calculate the Euclidean distance between nodes i and j
            dx = x(i) - x(j);  % Difference in x-coordinates
            dy = y(i) - y(j);  % Difference in y-coordinates
            euclid_distance = sqrt(dx^2 + dy^2);  % Euclidean distance between nodes i and j

            % Calculate delay based on Euclidean distance and delay factor
            delay_matrix(i, j) = delay_factor * euclid_distance;
        end
    end
end

%% ASK JOHN ABOUT THIS PART - DEFINING SELF-CONNECTIONS
for i = 1:S.networkSize.nodes,
   if any(selfFeedback == i),
       net(i,i) = 1;
   else
       net(i,i) = 0;
   end
end

%% Run Connection String
con = connectionString(net, selfFeedback, stimulusListeners, delay_matrix);

%% Use the 4th event that we just created, write about it
ep{2}.type = 'network';
ep{2}.description = '100ftdelay';
ep{2}.diameter = 10;
ep{2}.connections = con;
ep{2}.startTime = 10;
ep{2}.isDelay = true;

netType = 'ring';
clear net
switch netType,
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
        selfFeedback = 1:S.networkSize.nodes;  % include self feedback

    otherwise
        error('unknown network type')
end
for i = 1:S.networkSize.nodes,
   if any(selfFeedback == i),
       net(i,i) = 1;
   else
       net(i,i) = 0;
   end
end
% Define the connection strings
con1 = connectionString(net, selfFeedback, stimulusListeners, delay_matrix);
ep{3}.type = 'network'; % define first stimulus
ep{3}.description = 'ring';
ep{3}.connections = con1
ep{3}.startTime = 0; 



S.experimentParts = ep;

disp(con)
disp(delay_matrix)

xml = structToXMLPlist(S);
saveXMLPlist(fname, S)




