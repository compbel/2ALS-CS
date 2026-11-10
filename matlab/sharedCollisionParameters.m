function collision = sharedCollisionParameters(lookup, d, z)
%SHAREDCOLLISIONPARAMETERS Resolve one family's collision set by qID.
%   COLLISION = SHAREDCOLLISIONPARAMETERS(LOOKUP,D,Z) performs no rational
%   collision arithmetic. It collects and deduplicates qIDs from the lookup,
%   then returns the exact reduced numerators and denominators in the same
%   rational order used by direct collision construction.

validateLookup(lookup);
if ~(isnumeric(d) || islogical(d)) || ~isscalar(d) || ...
        ~isfinite(d) || d ~= floor(d) || d < 1 || d > 3
    error('padic:sharedCollisionDepth', 'd must be one of 1, 2, or 3.');
end
d = double(d);
z = exactZVector(z, lookup.nTaxa, d);

pairCount = lookup.nTaxa * (lookup.nTaxa - 1) / 2;
ids = zeros(pairCount, 1, 'uint32');
count = 0;
for i = 1:lookup.nTaxa - 1
    for j = i + 1:lookup.nTaxa
        delta = z(i) - z(j);
        qId = lookup.qIds(i, j, d, delta + 8);
        if qId ~= lookup.noCollision
            count = count + 1;
            ids(count) = qId;
        end
    end
end

if count == 0
    qIds = zeros(0, 1, 'uint32');
    numerators = zeros(0, 1);
    denominators = zeros(0, 1);
else
    ids = unique(ids(1:count));
    numerators = lookup.numerators(double(ids));
    denominators = lookup.denominators(double(ids));
    order = rationalOrder(numerators, denominators);
    qIds = ids(order);
    numerators = numerators(order);
    denominators = denominators(order);
end
collision = struct('valid', true, 'values', numerators ./ denominators, ...
    'numerators', numerators, 'denominators', denominators, ...
    'qIds', qIds, ...
    'degeneratePairs', zeros(0, 2));
end

function validateLookup(lookup)
required = {'embedding', 'nTaxa', 'qIds', 'numerators', ...
    'denominators', 'noCollision', 'globalRadixTree'};
if ~isstruct(lookup) || ~isscalar(lookup) || ...
        ~all(isfield(lookup, required))
    error('padic:sharedCollisionLookup', ...
        'lookup is not a shared collision lookup.');
end
if ~isequal(size(lookup.qIds), ...
        [lookup.nTaxa, lookup.nTaxa, 3, 15]) || ...
        numel(lookup.numerators) ~= numel(lookup.denominators)
    error('padic:sharedCollisionLookup', ...
        'The shared collision lookup has inconsistent dimensions.');
end
end

function values = exactZVector(input, n, d)
if ~(isnumeric(input) || islogical(input)) || numel(input) ~= n
    error('padic:sharedCollisionZ', ...
        'z must contain one machine integer per taxon.');
end
values = double(input(:));
if ~isreal(values) || any(~isfinite(values)) || ...
        any(values ~= fix(values)) || any(values < 0) || any(values >= 2^d)
    error('padic:sharedCollisionZ', ...
        'Every z value must be an integer in [0,2^d).');
end
end

function order = rationalOrder(numerators, denominators)
order = 1:numel(numerators);
for left = 1:numel(order)
    minimum = left;
    for right = left + 1:numel(order)
        if numerators(order(right)) * denominators(order(minimum)) < ...
                numerators(order(minimum)) * denominators(order(right))
            minimum = right;
        end
    end
    order([left, minimum]) = order([minimum, left]);
end
end
