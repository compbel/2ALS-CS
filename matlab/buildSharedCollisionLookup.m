function lookup = buildSharedCollisionLookup(embedding)
%BUILDSHAREDCOLLISIONLOOKUP Precompute exact CS-family collision parameters.
%   LOOKUP = BUILDSHAREDCOLLISIONLOOKUP(S) builds the immutable lookup for
%   every taxon pair, d=1,2,3, and delta=-(2^d-1),...,2^d-1. A zero qID
%   denotes no collision in Z_2. Distinct reduced rational collision values
%   are stored once in LOOKUP.numerators and LOOKUP.denominators.

embedding = exactIntegerVector(embedding);
n = numel(embedding);
if n < 2 || numel(unique(embedding)) ~= n
    error('padic:sharedCollisionEmbedding', ...
        'The shared collision lookup requires distinct coordinates.');
end

maximumDifference = max(embedding) - min(embedding);
if maximumDifference > floor((flintmax - 7) / 7)
    error('padic:sharedCollisionRange', ...
        'The embedding is too large for exact shared collision arithmetic.');
end

% Delta is addressed by delta+8, so the same compact dense array supports
% all three depths. Unused entries remain the NO_COLLISION value zero.
qIds = zeros(n, n, 3, 15, 'uint32');
numerators = zeros(0, 1);
denominators = zeros(0, 1);
keyToId = containers.Map('KeyType', 'char', 'ValueType', 'uint32');

for i = 1:n - 1
    for j = i + 1:n
        deltaS = embedding(i) - embedding(j);
        for d = 1:3
            maximumDelta = 2^d - 1;
            factor = maximumDelta;
            for delta = -maximumDelta:maximumDelta
                denominator = factor * deltaS + delta;
                if denominator == 0
                    continue;
                end
                numerator = -deltaS;
                common = gcd(abs(numerator), abs(denominator));
                numerator = numerator / common;
                denominator = denominator / common;
                if denominator < 0
                    numerator = -numerator;
                    denominator = -denominator;
                end
                if integerValuation(numerator) < ...
                        integerValuation(denominator)
                    continue;
                end

                key = sprintf('%.0f/%.0f', numerator, denominator);
                if isKey(keyToId, key)
                    qId = keyToId(key);
                else
                    numerators(end + 1, 1) = numerator; %#ok<AGROW>
                    denominators(end + 1, 1) = denominator; %#ok<AGROW>
                    qId = uint32(numel(numerators));
                    keyToId(key) = qId;
                end
                qIds(i, j, d, delta + 8) = qId;
            end
        end
    end
end

if ~isempty(numerators) && max(abs(numerators)) * max(denominators) > ...
        flintmax / 2
    error('padic:sharedCollisionRange', ...
        'Reduced collision rationals exceed the exact cross-product bound.');
end
if isempty(numerators)
    globalRadixTree = struct.empty(0, 1);
else
    allCollisions = struct('numerators', numerators, ...
        'denominators', denominators);
    globalRadixTree = integerCollisionRadixTree(allCollisions, false);
end
lookup = struct('embedding', embedding, 'nTaxa', n, ...
    'qIds', qIds, 'numerators', numerators, ...
    'denominators', denominators, 'noCollision', uint32(0), ...
    'globalRadixTree', globalRadixTree);
end

function values = exactIntegerVector(input)
if isnumeric(input) || islogical(input)
    values = double(input(:));
    usable = isreal(values) && all(isfinite(values)) && ...
        all(values == fix(values)) && all(abs(values) < 2^49);
else
    [values, usable] = exactIntegerDoubles(sym(input(:)));
    usable = usable && all(abs(values) < 2^49);
end
if ~usable
    error('padic:sharedCollisionEmbedding', ...
        'Shared collision coordinates must be exact bounded integers.');
end
end

function valuation = integerValuation(value)
value = abs(value);
if value == 0
    error('padic:v2Zero', 'v2(0) is undefined.');
end
valuation = 0;
while mod(value, 2) == 0
    value = value / 2;
    valuation = valuation + 1;
end
end
