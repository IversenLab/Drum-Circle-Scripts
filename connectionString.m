function connStr = connectionString(net, selfFeedback, stimulusListeners, delay_matrix)
% connectionString  construct connection list for Rhythm Network
%
%   connStr = connectionString(matrix, stimulusListeners);
%
%   INPUT
%       matrix          n x n connection matrix, possibly weighted [0 1]
%                       need not be symmetric. diagonal is ignored)
%       stimulusListeners   cell array of m x 2 matrices indicating which nodes
%                       can hear which stimulus, with weight. e.g. [node weight;...]
%                       [1  1.0; 2 0.5] -> node 1, weight 1, node 2, weight 0.5
%
%   OUTPUT
%       connStr         cell array of connection strings in format
%                       {from, to}weight, e.g. {3, 5}1.0
%                       And to indicate stimulus channel, node 0, hundreds
%                       {0.01, 2}1.0 --> channel 1 to node 2, weight 1
%
%   JRI 3/22/05


connStr = {};

%make stimulus listeners
if nargin >= 2,
    for iStim = 1:length(stimulusListeners),
        conn = stimulusListeners{iStim};
        for iNode = 1:size(conn,1),
            connStr{end+1} = sprintf('{%.2f, %d}%f',iStim/100, conn(iNode,1), conn(iNode,2));
        end
    end
end

%regular connections
for i = 1:size(net,1),
    for j = 1:size(net,2),
        weight = net(i,j);
        if weight > 0,
            delay = delay_matrix(i,j);
            connStr{end+1} = sprintf('{%d, %d}%f,%f', i, j, weight, delay);
        end
    end
end



