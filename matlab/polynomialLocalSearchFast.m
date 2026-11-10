function result = polynomialLocalSearchFast(initialTree, objectiveFunction, ...
        familySpecification, varargin)
%POLYNOMIALLOCALSEARCHFAST Embedding-scored Camin-Sokal local search.
%   Candidate scores are computed directly from pairwise 2-adic valuations.
%   No candidate tree is constructed. The selected improving embedding is
%   converted to a tree once, after the complete family scan, and only that
%   accepted tree is canonically re-embedded for the next iteration.
%
%   FAMILIES may be a fixed structure array with unary function handles .g
%   and .h (and optional .id), or an adaptive provider handle @(state).
%
%   Supported name-value options:
%     'MaxIterations'          positive integer or Inf (default Inf)
%     'ImprovementTolerance'   nonnegative scalar (default 0)
%     'ValidateNeighborhood'   logical scalar (default false)
%     'NeighborhoodArithmetic' 'auto' or 'integer' (default auto)
%     'ParallelFamilies'       logical scalar (default true)
%     'ParallelWorkers'        exact requested worker count (default 4)
%     'ParallelPoolType'       'processes' or 'threads' (default processes)
%     'EscapeFamilies'         secondary families used at a primary optimum
%     'CaminSokalAlignment'    required binary alignment
%
%   OBJECTIVE must be a function handle but candidate scores are computed
%   directly from CAMINSOKALALIGNMENT.
%   INITIALTREE may be a rooted binary digraph or compact topology structure.

parser = inputParser;
parser.FunctionName = 'polynomialLocalSearchFast';
addParameter(parser, 'MaxIterations', inf, @validMaximum);
addParameter(parser, 'ImprovementTolerance', 0, @validTolerance);
addParameter(parser, 'ValidateNeighborhood', false, ...
    @(value) islogical(value) && isscalar(value));
addParameter(parser, 'NeighborhoodArithmetic', 'auto', ...
    @(value) (ischar(value) || (isstring(value) && isscalar(value))) && ...
    any(strcmpi(string(value), ["auto", "integer", "symbolic"])));
addParameter(parser, 'ParallelFamilies', true, ...
    @(value) islogical(value) && isscalar(value));
addParameter(parser, 'ParallelWorkers', 4, ...
    @(value) isnumeric(value) && isscalar(value) && isreal(value) && ...
    isfinite(value) && value >= 1 && value == floor(value));
addParameter(parser, 'ParallelPoolType', 'processes', ...
    @(value) (ischar(value) || (isstring(value) && isscalar(value))) && ...
    any(strcmpi(string(value), ["processes", "threads"])));
addParameter(parser, 'EscapeFamilies', struct.empty(0, 1), ...
    @(value) isstruct(value) || isempty(value));
addParameter(parser, 'CaminSokalAlignment', [], @validAlignment);
parse(parser, varargin{:});
options = parser.Results;
options.NeighborhoodArithmetic = char(lower(string(options.NeighborhoodArithmetic)));
options.ParallelPoolType = char(lower(string(options.ParallelPoolType)));
if isempty(options.CaminSokalAlignment)
    error('padic:fastSearchAlignment', ...
        'polynomialLocalSearchFast requires CaminSokalAlignment.');
end
if strcmp(options.NeighborhoodArithmetic, 'symbolic')
    error('padic:fastSearchIntegerRequired', ...
        'polynomialLocalSearchFast supports exact machine-integer families only.');
end
if ~isa(objectiveFunction, 'function_handle')
    error('padic:objectiveHandle', 'objectiveFunction must be a function handle.');
end
[familyMode, fixedFamilies, familyProvider] = ...
    normalizeFamilySpecification(familySpecification);
escapeFamilies = normalizeOptionalFamilies(options.EscapeFamilies);
escapeEnabled = ~isempty(escapeFamilies);
options.EscapeFamilies = struct.empty(0, 1);

current = normalizeInitialTree(initialTree);
currentKey = current.key;
embedding = canonicalFullEmbeddingCompact(current);
scorePlan = caminSokalEmbeddingPlan(options.CaminSokalAlignment);
currentScore = caminSokalScoreEmbedding(embedding, scorePlan);
if currentScore ~= caminSokalScoreKey(currentKey, options.CaminSokalAlignment)
    error('padic:embeddingScoreMismatch', ...
        'Initial embedding and topology scores disagree.');
end
initialTopology = current;
initialKey = currentKey;
initialScore = currentScore;
trajectory = struct.empty(0, 1);
acceptedMoves = struct.empty(0, 1);
iteration = 0;
converged = false;
parallelEverUsed = false;
maximumParallelWorkers = 0;
parallelPoolTypeUsed = '';
searchPhase = 'primary';
sharedCollisionLookup = struct.empty(0, 1);
sharedCollisionEmbedding = zeros(0, 1);

while true
    if iteration >= options.MaxIterations
        terminationReason = 'maximum_iterations';
        break;
    end
    iteration = iteration + 1;
    state = struct( ...
        'iteration', iteration, ...
        'currentTopology', current, ...
        'currentTopologyData', current, ...
        'currentKey', currentKey, ...
        'currentScore', currentScore, ...
        'embedding', embedding, ...
        'initialTopology', initialTopology, ...
        'initialKey', initialKey, ...
        'initialScore', initialScore, ...
        'trajectory', trajectory, ...
        'acceptedMoves', acceptedMoves, ...
        'familyPhase', searchPhase, ...
        'options', options);
    if strcmp(searchPhase, 'escape')
        families = escapeFamilies;
    elseif strcmp(familyMode, 'fixed')
        families = fixedFamilies;
    else
        families = normalizeFamilies(familyProvider(state));
        assertNonemptyFamilies(families, iteration);
    end
    if familiesNeedSharedCollisionLookup(families)
        if isempty(sharedCollisionLookup) || ...
                ~isequal(sharedCollisionEmbedding, embedding)
            sharedCollisionLookup = buildSharedCollisionLookup(embedding);
            sharedCollisionEmbedding = embedding;
        end
    end
    familyBestCandidates = struct.empty(0, 1);
    [useParallel, parallelWorkerCount, iterationPoolType] = ...
        prepareParallelExecution(options.ParallelFamilies, ...
        options.ParallelWorkers, options.ParallelPoolType, numel(families));
    parallelEverUsed = parallelEverUsed || useParallel;
    maximumParallelWorkers = max(maximumParallelWorkers, parallelWorkerCount);
    if useParallel
        parallelPoolTypeUsed = iterationPoolType;
    end
    familyScans = cell(numel(families), 1);
    if useParallel
        parfor (familyIndex = 1:numel(families), parallelWorkerCount)
            familyScans{familyIndex} = scanFamilyFast( ...
                families(familyIndex), familyIndex, embedding, currentScore, ...
                scorePlan, options, sharedCollisionLookup);
        end
    else
        for familyIndex = 1:numel(families)
            familyScans{familyIndex} = scanFamilyFast( ...
                families(familyIndex), familyIndex, embedding, currentScore, ...
                scorePlan, options, sharedCollisionLookup);
        end
    end
    % Merge in input-family order so worker completion order cannot affect
    % ties or the retained replay witness.
    for familyIndex = 1:numel(families)
        scan = familyScans{familyIndex};
        if ~isempty(scan.bestCandidate)
            if isempty(familyBestCandidates)
                familyBestCandidates = scan.bestCandidate;
            else
                familyBestCandidates(end + 1, 1) = scan.bestCandidate; %#ok<AGROW>
            end
        end
    end

    [bestIndex, familyBestCandidates] = selectBestFast( ...
        familyBestCandidates, scorePlan);
    accepted = false;
    bestKey = currentKey;
    bestTopology = current;
    bestScore = currentScore;
    acceptedFamilyIndex = [];
    acceptedFamilyId = '';
    acceptedTheta = [];
    if ~isempty(bestIndex)
        best = familyBestCandidates(bestIndex);
        bestKey = best.key;
        bestScore = best.score;
        if bestScore < currentScore - options.ImprovementTolerance
            accepted = true;
            bestTopology = topologyDataFromCoordinates( ...
                best.transformedCoordinates);
            if ~isequal(bestTopology.key, bestKey) || ...
                    caminSokalScoreKey(bestTopology.key, ...
                    options.CaminSokalAlignment) ~= bestScore
                error('padic:embeddingScoreMismatch', ...
                    'Accepted embedding and reconstructed tree disagree.');
            end
            acceptedFamilyIndex = best.familyIndex;
            acceptedFamilyId = best.familyId;
            acceptedTheta = best.theta;
        end
    end
    continueToEscape = ~accepted && strcmp(searchPhase, 'primary') && ...
        escapeEnabled && iteration < options.MaxIterations;

    row = struct( ...
        'iteration', iteration, ...
        'familyPhase', searchPhase, ...
        'sourceKey', currentKey, ...
        'sourceTopology', current, ...
        'sourceScore', currentScore, ...
        'familiesScanned', numel(families), ...
        'parallelFamiliesUsed', useParallel, ...
        'parallelWorkerCount', parallelWorkerCount, ...
        'parallelPoolType', iterationPoolType, ...
        'bestKey', bestKey, ...
        'bestTopology', bestTopology, ...
        'bestScore', bestScore, ...
        'accepted', accepted, ...
        'acceptedFamilyIndex', acceptedFamilyIndex, ...
        'acceptedFamilyId', acceptedFamilyId, ...
        'theta', acceptedTheta, ...
        'continuedToEscape', continueToEscape, ...
        'completeScan', true);
    if isempty(trajectory)
        trajectory = row;
    else
        trajectory(end + 1, 1) = row; %#ok<AGROW>
    end

    if ~accepted
        if strcmp(searchPhase, 'primary') && escapeEnabled
            if iteration >= options.MaxIterations
                terminationReason = 'maximum_iterations';
                break;
            end
            searchPhase = 'escape';
            continue;
        end
        converged = true;
        terminationReason = 'local_optimum';
        break;
    end

    move = struct( ...
        'iteration', iteration, ...
        'sourceKey', currentKey, ...
        'sourceTopology', current, ...
        'sourceScore', currentScore, ...
        'sourceEmbedding', best.sourceEmbedding, ...
        'familyPhase', searchPhase, ...
        'familyIndex', best.familyIndex, ...
        'familyId', best.familyId, ...
        'gValues', best.gValues, ...
        'hValues', best.hValues, ...
        'theta', best.theta, ...
        'transformedCoordinates', best.transformedCoordinates, ...
        'targetKey', best.key, ...
        'targetTopology', bestTopology, ...
        'targetScore', best.score);
    if isempty(acceptedMoves)
        acceptedMoves = move;
    else
        acceptedMoves(end + 1, 1) = move; %#ok<AGROW>
    end
    current = bestTopology;
    currentKey = best.key;
    currentScore = best.score;
    embedding = canonicalFullEmbeddingCompact(current);
    sharedCollisionLookup = struct.empty(0, 1);
    sharedCollisionEmbedding = zeros(0, 1);
    searchPhase = 'primary';
end

result = struct( ...
    'initialTopology', initialTopology, ...
    'initialKey', initialKey, ...
    'initialScore', initialScore, ...
    'finalTopology', current, ...
    'finalKey', currentKey, ...
    'finalScore', currentScore, ...
    'iterations', iteration, ...
    'acceptedMoveCount', numel(acceptedMoves), ...
    'converged', converged, ...
    'terminationReason', terminationReason, ...
    'trajectory', trajectory, ...
    'acceptedMoves', acceptedMoves, ...
    'parallelFamiliesUsed', parallelEverUsed, ...
    'parallelWorkerCount', maximumParallelWorkers, ...
    'parallelPoolType', parallelPoolTypeUsed, ...
    'candidateScoring', 'embedding_pairwise_valuations', ...
    'options', options);
result = materializeReturnedTopologies(result);
end

function families = normalizeOptionalFamilies(input)
if isempty(input)
    families = struct.empty(0, 1);
else
    families = normalizeFamilies(input);
end
end

function scan = scanFamilyFast(family, familyIndex, embedding, ...
        currentScore, scorePlan, options, sharedLookup)
[useSharedLookup, collisionDepth, collisionZ] = ...
    familySharedCollisionMetadata(family, numel(embedding));
if ~useSharedLookup
    familyLookup = struct.empty(0, 1);
else
    familyLookup = sharedLookup;
end
try
    core = fastIntegerFamilyBest(embedding, family.g, family.h, ...
        scorePlan, currentScore, options.ImprovementTolerance, ...
        options.ValidateNeighborhood, familyLookup, ...
        collisionDepth, collisionZ);
catch cause
    if ~strcmp(cause.identifier, 'padic:integerFastPathUnavailable')
        rethrow(cause);
    end
    error('padic:fastSearchIntegerRequired', ...
        'Family %d is unavailable to the exact integer fast search: %s', ...
        familyIndex, cause.message);
end

bestCandidate = struct.empty(0, 1);
if ~isempty(core)
    bestCandidate = struct( ...
        'key', core.key, ...
        'keyKnown', core.keyKnown, ...
        'score', core.score, ...
        'familyIndex', familyIndex, ...
        'familyId', family.id, ...
        'gValues', core.gValues, ...
        'hValues', core.hValues, ...
        'theta', core.theta, ...
        'sourceEmbedding', embedding, ...
        'transformedCoordinates', core.transformedCoordinates);
end
scan = struct('bestCandidate', bestCandidate);
end

function result = familiesNeedSharedCollisionLookup(families)
result = false;
if ~isfield(families, 'collisionDepth') || ...
        ~isfield(families, 'collisionZ')
    return;
end
for index = 1:numel(families)
    if ~isempty(families(index).collisionDepth) || ...
            ~isempty(families(index).collisionZ)
        familySharedCollisionMetadata(families(index), []);
        result = true;
        return;
    end
end
end

function [usable, d, z] = familySharedCollisionMetadata(family, n)
usable = false;
d = [];
z = [];
if ~isfield(family, 'collisionDepth') || ...
        ~isfield(family, 'collisionZ')
    return;
end
if isempty(family.collisionDepth) && isempty(family.collisionZ)
    return;
end
d = double(family.collisionDepth);
z = double(family.collisionZ(:));
if ~isscalar(d) || ~isfinite(d) || d ~= floor(d) || d < 1 || d > 3
    error('padic:sharedCollisionMetadata', ...
        'collisionDepth must be one of 1, 2, or 3.');
end
if ~isempty(n) && numel(z) ~= n
    error('padic:sharedCollisionMetadata', ...
        'collisionZ must contain one value per taxon.');
end
if any(~isfinite(z)) || any(z ~= fix(z)) || ...
        any(z < 0) || any(z >= 2^d)
    error('padic:sharedCollisionMetadata', ...
        'collisionZ values must be integers in [0,2^collisionDepth).');
end
usable = true;
end

function [useParallel, workerCount, poolType] = prepareParallelExecution( ...
        requested, requestedWorkers, requestedPoolType, familyCount)
useParallel = false;
workerCount = 0;
poolType = '';
if ~requested || familyCount < 2
    return;
end
if isempty(ver('parallel')) || ~license('test', 'Distrib_Computing_Toolbox')
    error('padic:parallelUnavailable', ...
        ['Parallel Computing Toolbox is unavailable, so the requested ' ...
        '%d-worker pool cannot be created.'], requestedWorkers);
end
targetWorkers = requestedWorkers;
pool = gcp('nocreate');
if ~isempty(pool)
    existingPoolType = identifyPoolType(pool);
    if pool.NumWorkers ~= targetWorkers || ...
            ~strcmp(existingPoolType, requestedPoolType)
        delete(pool);
        pool = [];
    end
end
if isempty(pool)
    try
        profile = [upper(requestedPoolType(1)), requestedPoolType(2:end)];
        pool = parpool(profile, targetWorkers);
    catch cause
        error('padic:parallelPool', ...
            'Could not start the requested %d-worker %s pool: %s', ...
            targetWorkers, requestedPoolType, cause.message);
    end
end
poolType = identifyPoolType(pool);
workerCount = pool.NumWorkers;
if workerCount ~= targetWorkers || ~strcmp(poolType, requestedPoolType)
    error('padic:parallelPoolMismatch', ...
        ['MATLAB created a %d-worker %s pool instead of the requested ' ...
        '%d-worker %s pool.'], workerCount, poolType, targetWorkers, ...
        requestedPoolType);
end
useParallel = true;
end

function poolType = identifyPoolType(pool)
if isa(pool, 'parallel.ThreadPool')
    poolType = 'threads';
else
    poolType = 'processes';
end
end

function [mode, fixedFamilies, provider] = normalizeFamilySpecification(input)
if isstruct(input)
    mode = 'fixed';
    fixedFamilies = normalizeFamilies(input);
    assertNonemptyFamilies(fixedFamilies, []);
    provider = [];
elseif isa(input, 'function_handle')
    mode = 'adaptive';
    fixedFamilies = struct.empty(0, 1);
    provider = input;
else
    error('padic:familySpecification', ...
        'families must be a fixed struct array or an iteration provider handle.');
end
end

function assertNonemptyFamilies(families, iteration)
if ~isempty(families)
    return;
end
if isempty(iteration)
    error('padic:noFamilies', 'At least one function pair is required.');
else
    error('padic:noFamilies', ...
        'The family provider returned no function pairs at iteration %d.', iteration);
end
end

function families = normalizeFamilies(input)
if ~isstruct(input) || ~all(isfield(input, {'g', 'h'}))
    error('padic:familyFormat', ...
        'families must be a struct array containing g and h function handles.');
end
families = input(:);
for index = 1:numel(input)
    if isfield(input, 'id') && ~isempty(input(index).id)
        families(index).id = char(string(input(index).id));
    else
        families(index).id = sprintf('family_%d', index);
    end
    if ~isa(input(index).g, 'function_handle') || ...
            ~isa(input(index).h, 'function_handle')
        error('padic:functionFamilies', ...
            'Family %d must contain function handles g and h.', index);
    end
    families(index).g = input(index).g;
    families(index).h = input(index).h;
end
end

function topology = normalizeInitialTree(input)
if isa(input, 'digraph')
    validateRootedDigraph(input);
    topology = topologyDataFromCoordinates( ...
        canonicalFullEmbedding(input, false));
elseif isstruct(input) && isscalar(input) && ...
        all(isfield(input, {'sources', 'targets', 'taxon', 'key'}))
    topology = topologyDataFromCoordinates( ...
        canonicalFullEmbeddingCompact(input));
else
    error('padic:initialTreeFormat', ...
        'initialTree must be a rooted binary digraph or compact topology data.');
end
end

function [bestIndex, candidates] = selectBestFast(candidates, scorePlan)
bestIndex = [];
for index = 1:numel(candidates)
    if isempty(bestIndex) || candidates(index).score < candidates(bestIndex).score
        bestIndex = index;
    elseif candidates(index).score == candidates(bestIndex).score
        candidates(index) = ensureCandidateKey(candidates(index), scorePlan);
        candidates(bestIndex) = ensureCandidateKey( ...
            candidates(bestIndex), scorePlan);
        if lexicographicallyLess( ...
                candidates(index).key, candidates(bestIndex).key)
            bestIndex = index;
        end
    end
end
if ~isempty(bestIndex)
    candidates(bestIndex) = ensureCandidateKey(candidates(bestIndex), scorePlan);
end
end

function candidate = ensureCandidateKey(candidate, scorePlan)
if candidate.keyKnown
    return;
end
[score, valuations] = caminSokalScoreEmbeddingPrepared( ...
    candidate.transformedCoordinates, scorePlan);
if score ~= candidate.score
    error('padic:embeddingScoreMismatch', ...
        'Repeated direct scoring changed a candidate score.');
end
candidate.key = topologyKeyFromValuations(valuations);
candidate.keyKnown = true;
end

function result = lexicographicallyLess(left, right)
leftCodes = double(topologyKeyText(left));
rightCodes = double(topologyKeyText(right));
limit = min(numel(leftCodes), numel(rightCodes));
different = find(leftCodes(1:limit) ~= rightCodes(1:limit), 1);
if isempty(different)
    result = numel(leftCodes) < numel(rightCodes);
else
    result = leftCodes(different) < rightCodes(different);
end
end

function result = validMaximum(value)
result = isscalar(value) && isreal(value) && ...
    ((isinf(value) && value > 0) || (value >= 1 && value == floor(value)));
end

function result = validTolerance(value)
result = isnumeric(value) && isscalar(value) && isreal(value) && ...
    isfinite(value) && value >= 0;
end

function result = validAlignment(value)
result = isempty(value) || ((isnumeric(value) || islogical(value)) && ...
    ismatrix(value) && ~isempty(value) && all(value(:) == 0 | value(:) == 1));
end

function result = materializeReturnedTopologies(result)
result.initialTopologyData = result.initialTopology;
result.finalTopologyData = result.finalTopology;
result.initialTopology = topologyDataToDigraph(result.initialTopologyData);
result.finalTopology = topologyDataToDigraph(result.finalTopologyData);
for index = 1:numel(result.trajectory)
    sourceData = result.trajectory(index).sourceTopology;
    bestData = result.trajectory(index).bestTopology;
    result.trajectory(index).sourceTopologyData = sourceData;
    result.trajectory(index).bestTopologyData = bestData;
    result.trajectory(index).sourceTopology = topologyDataToDigraph(sourceData);
    result.trajectory(index).bestTopology = topologyDataToDigraph(bestData);
end
for index = 1:numel(result.acceptedMoves)
    sourceData = result.acceptedMoves(index).sourceTopology;
    targetData = result.acceptedMoves(index).targetTopology;
    result.acceptedMoves(index).sourceTopologyData = sourceData;
    result.acceptedMoves(index).targetTopologyData = targetData;
    result.acceptedMoves(index).sourceTopology = topologyDataToDigraph(sourceData);
    result.acceptedMoves(index).targetTopology = topologyDataToDigraph(targetData);
end
end
