function [score, valuations] = caminSokalScoreEmbedding(coordinates, alignmentOrPlan)
%CAMINSOKALSCOREEMBEDDING Camin-Sokal score directly from a 2-adic embedding.
%   For each character c, the contribution is the number of unique rows in
%   V(I_c,J_c), where V(i,j)=v2(s_i-s_j). Repeated alignment columns are
%   grouped once by CAMINSOKALEMBEDDINGPLAN.

if isstruct(alignmentOrPlan) && isscalar(alignmentOrPlan) && ...
        all(isfield(alignmentOrPlan, {'nTaxa', 'patternWeights', ...
        'oneIndices', 'zeroIndices'}))
    plan = alignmentOrPlan;
else
    plan = caminSokalEmbeddingPlan(alignmentOrPlan);
end
coordinates = exactIntegerVector(coordinates);
n = numel(coordinates);
if n ~= plan.nTaxa
    error('padic:alignmentTaxa', ...
        'The embedding and alignment must have the same number of taxa.');
end
if numel(unique(coordinates)) ~= n
    error('padic:coordinateCollision', ...
        'Embedding coordinates must be distinct.');
end
[score, valuations] = caminSokalScoreEmbeddingPrepared(coordinates, plan);
end

function coordinates = exactIntegerVector(input)
if isnumeric(input) || islogical(input)
    coordinates = double(input(:));
    usable = isreal(coordinates) && all(isfinite(coordinates)) && ...
        all(coordinates == fix(coordinates)) && all(abs(coordinates) < 2^50);
else
    [coordinates, usable] = exactIntegerDoubles(sym(input(:)));
    usable = usable && all(abs(coordinates) < 2^50);
end
if ~usable
    error('padic:embeddingScoreInteger', ...
        'Embedding scoring requires exact machine integers below 2^50.');
end
end
