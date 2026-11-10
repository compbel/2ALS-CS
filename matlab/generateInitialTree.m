function [tree, source] = generateInitialTree(n, userTree, treeModel)
%GENERATEINITIALTREE Validate a supplied tree or generate an independent one.
%   TREE = GENERATEINITIALTREE(N, USERTREE, MODEL) returns USERTREE after
%   validating that it is a full rooted binary digraph on taxa 1..N. If
%   USERTREE is an empty digraph, MODEL selects a new random tree:
%       'YULE'    - Yule-Harding topology
%       'UNIFORM' - uniform rooted binary labeled topology
%   This function is independent of instance and alignment simulation.

if ~isscalar(n) || ~isnumeric(n) || ~isreal(n) || ...
        ~isfinite(n) || n < 2 || n ~= floor(n)
    error('padic:taxonCount', 'n must be an integer at least 2.');
end
if nargin < 2 || isempty(userTree)
    userTree = digraph();
end
if nargin < 3 || isempty(treeModel)
    treeModel = 'YULE';
end
if ~isa(userTree, 'digraph')
    error('padic:initialTreeFormat', 'userTree must be a digraph.');
end

if numnodes(userTree) > 0
    information = validateRootedDigraph(userTree);
    if information.n ~= n
        error('padic:initialTreeTaxonCount', ...
            'The supplied tree has %d taxa but n is %d.', information.n, n);
    end
    tree = userTree;
    source = 'user';
    return;
end

model = upper(strrep(char(string(treeModel)), '-', '_'));
switch model
    case 'YULE'
        tree = randomYuleRootedTopology(n);
    case 'UNIFORM'
        tree = randomRootedTopology(n);
    otherwise
        error('padic:treeModel', 'treeModel must be YULE or UNIFORM.');
end
source = 'generated';
end
