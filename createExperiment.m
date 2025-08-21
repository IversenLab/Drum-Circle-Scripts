% create network definition files for RhythmNetwork experiment
%
%   saves a .netdef file that defines an experiment
%
%   an example file is appended to end
%
%   constraints: must have some network defined before stimuli are defined
%
%   JRI 3/16/05

%file naming
fpath = '/Users/jri/';
fextension = '.netdef';
fname = [fpath 'matlabtest' fextension];

%format
clear S

S.description = 'Test experiment matlab';
S.experimentDuration = 20;
S.notes = '';


S.networkSize.nodes = 6;
S.networkSize.stimulusChannels = 2;

%define experiment parts

ep{1}.type = 'network';
ep{1}.description = 'No Connections (record only)';
ep{1}.connections = {[]};
ep{1}.startTime = 0;

ep{2}.type = 'stimulus';
ep{2}.description = 'test stim 1';
ep{2}.startTime = 0;
ep{2}.stimulus = '1: 16(64), IOI=1000.0, events=20, phase=0.0, jitter=500.0';

ep{3}.type = 'stimulus';
ep{3}.description = 'test stim 2';
ep{3}.startTime = 0;
ep{3}.stimulus = '2: 15(64), IOI=700.0, events=28, phase=400.0';

%%define a network
%define which nodes hear which stimulus
stimulusListeners{1} = [1 1.0; 2 0.5];
stimulusListeners{2} = [3 1.0];

%set self connections
%nodes with selfFeedback
selfFeedback = [1 3 5]; 

%get a connectivity matrix
netType = 'ring';
clear net

switch netType,
    case 'random',
        thresh = 0.5;
        net = (rand(S.networkSize.nodes) > thresh);

    case 'none',
        network = zeros(S.networkSize.nodes);
        selfFeedback = [];
        
    case 'self'
        net = zeros(S.networkSize.nodes);
        selfFeedback = 1:S.networkSize.nodes;
        
    case 'ring'
        net = zeros(S.networkSize.nodes);
        for i = 1:S.networkSize.nodes,
            from = i;
            to = mod(from,S.networkSize.nodes)+1;
            net(from,to) = 1;
        end

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

con = connectionString(net, selfFeedback, stimulusListeners);

ep{4}.type = 'network';
ep{4}.description = 'random net';
ep{4}.connections = {'{0.01, 1}', '{0.02,2}', '{2, 3}', '{3, 4}', '{4, 4}', '{4, 5}', '{4, 6}'};
ep{4}.startTime = 10;

S.experimentParts = ep;

xml = structToXMLPlist(S);
saveXMLPlist(fname, S)




