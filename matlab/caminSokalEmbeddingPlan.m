function plan = caminSokalEmbeddingPlan(alignment)
%CAMINSOKALEMBEDDINGPLAN Precompute alignment subsets for embedding scoring.

if ~(isnumeric(alignment) || islogical(alignment)) || ...
        ~ismatrix(alignment) || isempty(alignment) || ...
        any(alignment(:) ~= 0 & alignment(:) ~= 1)
    error('padic:binaryAlignment', ...
        'alignment must be a nonempty binary matrix.');
end
alignment = logical(alignment);
[patterns, ~, patternIndex] = unique(alignment.', 'rows');
allZero = ~any(patterns, 2);
allOne = all(patterns, 2);
singletonOne = sum(patterns, 2) == 1;
variablePatterns = find(~allZero & ~allOne & ~singletonOne);
allPatternWeights = accumarray(patternIndex, 1, [size(patterns, 1), 1]);
patternWeights = allPatternWeights(variablePatterns);
oneIndices = cell(numel(variablePatterns), 1);
zeroIndices = cell(numel(variablePatterns), 1);
for index = 1:numel(variablePatterns)
    pattern = variablePatterns(index);
    oneIndices{index} = find(patterns(pattern, :)).';
    zeroIndices{index} = find(~patterns(pattern, :)).';
end
plan = struct('nTaxa', size(alignment, 1), ...
    'siteCount', size(alignment, 2), ...
    'patternCount', size(patterns, 1), ...
    'scoredPatternCount', numel(variablePatterns), ...
    'constantScore', sum(allPatternWeights(allOne | singletonOne)), ...
    'patternWeights', patternWeights, ...
    'oneIndices', {oneIndices}, 'zeroIndices', {zeroIndices});
end
