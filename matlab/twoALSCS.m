function result = twoALSCS(alignment, varargin)
%TWOALSCS Rooted Camin-Sokal inference by exact 2-adic local search.
%   RESULT = TWOALSCS(M) minimizes rooted Camin-Sokal parsimony for the
%   binary taxon-by-character matrix M. The default start is a random
%   Yule-Harding tree and the default neighborhood contains both single-
%   column and allelic-state column-pair transformations at every iteration.
%
%   Name-value options:
%     InitialTree          rooted binary digraph; empty generates a tree
%     InitialTreeModel     'YULE' (default) or 'UNIFORM'
%     PairMode             'ALWAYS' (default), 'LOCAL', or 'NEVER'
%     MaxIterations        positive integer or Inf (default Inf)
%     ParallelFamilies     logical scalar (default true)
%     ParallelWorkers      exact worker count (default 4)
%     ParallelPoolType     'processes' (default) or 'threads'
%     ValidateNeighborhood enable expensive internal checks (default false)
%
%   The search accepts strict improvements only and scans every active
%   transformation family completely. It uses no restarts, plateau moves,
%   ratchets, annealing, or hybrid moves.

alignment = validateAlignment(alignment);
n = size(alignment, 1);

parser = inputParser;
parser.FunctionName = 'twoALSCS';
addParameter(parser, 'InitialTree', digraph(), ...
    @(value) isempty(value) || (isa(value, 'digraph') && isscalar(value)));
addParameter(parser, 'InitialTreeModel', 'YULE', @textScalar);
addParameter(parser, 'PairMode', 'ALWAYS', @textScalar);
addParameter(parser, 'MaxIterations', inf, @validMaximum);
addParameter(parser, 'ParallelFamilies', true, @logicalScalar);
addParameter(parser, 'ParallelWorkers', 4, @positiveInteger);
addParameter(parser, 'ParallelPoolType', 'processes', @validPoolType);
addParameter(parser, 'ValidateNeighborhood', false, @logicalScalar);
parse(parser, varargin{:});
options = parser.Results;

pairMode = normalizeColumnPairFamilyMode(options.PairMode);
[initialTree, initialTreeSource] = generateInitialTree( ...
    n, options.InitialTree, options.InitialTreeModel);
instance = struct('alignment', alignment);
[primaryFamilies, escapeFamilies, singleFamilies, pairFamilies] = ...
    configureCSFunctionFamilies(instance, ...
        @allelicStateTreeCSFunctionFamilies, pairMode);
objective = @(tree) caminSokalScore(tree, alignment);

result = polynomialLocalSearchFast(initialTree, objective, primaryFamilies, ...
    'EscapeFamilies', escapeFamilies, ...
    'MaxIterations', options.MaxIterations, ...
    'ImprovementTolerance', 0, ...
    'ValidateNeighborhood', options.ValidateNeighborhood, ...
    'NeighborhoodArithmetic', 'integer', ...
    'ParallelFamilies', options.ParallelFamilies, ...
    'ParallelWorkers', options.ParallelWorkers, ...
    'ParallelPoolType', options.ParallelPoolType, ...
    'CaminSokalAlignment', alignment);

result.algorithm = '2ALS-CS';
result.version = twoALSVersion();
result.initialTreeSource = initialTreeSource;
result.pairMode = pairMode;
result.taxonCount = n;
result.characterCount = size(alignment, 2);
result.singleColumnFamilyCount = numel(singleFamilies);
result.pairFamilyCount = numel(pairFamilies);
result.initialTree = result.initialTopology;
result.finalTree = result.finalTopology;
end

function alignment = validateAlignment(input)
if ~(isnumeric(input) || islogical(input)) || ~ismatrix(input) || ...
        isempty(input) || size(input, 1) < 2 || size(input, 2) < 1 || ...
        any(input(:) ~= 0 & input(:) ~= 1)
    error('twoALS:binaryAlignment', ...
        ['alignment must be a nonempty binary taxon-by-character matrix ' ...
         'with at least two taxa.']);
end
alignment = logical(input);
end

function result = textScalar(value)
result = (ischar(value) && isrow(value)) || ...
    (isstring(value) && isscalar(value));
end

function result = logicalScalar(value)
result = islogical(value) && isscalar(value);
end

function result = positiveInteger(value)
result = isnumeric(value) && isscalar(value) && isreal(value) && ...
    isfinite(value) && value >= 1 && value == floor(value);
end

function result = validMaximum(value)
result = isnumeric(value) && isscalar(value) && isreal(value) && ...
    ((isinf(value) && value > 0) || ...
    (isfinite(value) && value >= 1 && value == floor(value)));
end

function result = validPoolType(value)
result = textScalar(value) && ...
    any(strcmpi(string(value), ["processes", "threads"]));
end
