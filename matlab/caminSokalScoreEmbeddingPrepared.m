function [score, valuations] = caminSokalScoreEmbeddingPrepared(coordinates, plan)
%CAMINSOKALSCOREEMBEDDINGPREPARED Unchecked inner embedding-score kernel.
%   The fast local search guarantees distinct bounded integer coordinates and
%   supplies a validated CAMINSOKALEMBEDDINGPLAN.

n = numel(coordinates);
differences = abs(coordinates - transpose(coordinates));
unsignedDifferences = uint64(differences);
unsignedDifferences(1:n+1:end) = uint64(1);
lowbit = bitxor(unsignedDifferences, ...
    bitand(unsignedDifferences, unsignedDifferences - uint64(1)));
valuations = log2(double(lowbit));
valuations(1:n+1:end) = inf;

score = double(plan.constantScore);
patternWeights = plan.patternWeights;
oneIndexSets = plan.oneIndices;
zeroIndexSets = plan.zeroIndices;
for pattern = 1:numel(patternWeights)
    oneTaxa = oneIndexSets{pattern};
    zeroTaxa = zeroIndexSets{pattern};
    oneCount = numel(oneTaxa);
    zeroCount = numel(zeroTaxa);

    % Avoid materializing valuations(oneTaxa, zeroTaxa). This kernel is
    % called for every parameter region, and for these small taxon subsets a
    % JIT-compiled scalar reduction is substantially cheaper than allocating
    % and reducing a new indexed matrix for every character pattern.
    zeroThreshold = zeros(oneCount, 1);
    for position = 1:oneCount
        threshold = 0;
        taxon = oneTaxa(position);
        for zeroPosition = 1:zeroCount
            threshold = max(threshold, ...
                valuations(taxon, zeroTaxa(zeroPosition)));
        end
        zeroThreshold(position) = threshold;
    end

    representativePositions = zeros(1, oneCount);
    contribution = 0;
    for position = 1:oneCount
        taxon = oneTaxa(position);
        duplicate = false;
        for existing = 1:contribution
            representativePosition = representativePositions(existing);
            if valuations(taxon, ...
                    oneTaxa(representativePosition)) > ...
                    max(zeroThreshold(position), ...
                    zeroThreshold(representativePosition))
                duplicate = true;
                break;
            end
        end
        if ~duplicate
            contribution = contribution + 1;
            representativePositions(contribution) = position;
        end
    end
    score = score + patternWeights(pattern) * contribution;
end
end
