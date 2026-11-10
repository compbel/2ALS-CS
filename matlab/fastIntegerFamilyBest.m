function [best, generatedEmbeddings, diagnostic] = fastIntegerFamilyBest( ...
        embedding, g, h, scorePlan, currentScore, improvementTolerance, ...
        validate, sharedLookup, collisionDepth, collisionZ)
%FASTINTEGERFAMILYBEST Score one family directly from candidate embeddings.
%   No candidate tree, compact topology, or canonical embedding is built.
%   Nonimproving candidates are discarded immediately. A topology key is
%   derived from a valuation matrix only when required to reproduce the
%   deterministic lexicographic tie rule.
%
%   When G is the identity on the source embedding, the unique parameter
%   region containing theta=0 is omitted during representative generation:
%   it can only reproduce the source topology. [~, EMBEDDINGS] =
%   FASTINTEGERFAMILYBEST(..., [], ...) is a diagnostic mode. It returns the
%   remaining representative embeddings as an n-by-r matrix and does not
%   score candidates. Production local search never requests this output and
%   therefore pays no candidate-retention cost.

if nargin < 7
    validate = false;
end
if nargin < 8
    sharedLookup = struct.empty(0, 1);
end
if nargin < 9
    collisionDepth = [];
end
if nargin < 10
    collisionZ = [];
end
diagnostic = struct.empty(0, 1);
embedding = integerVector(embedding, 'embedding', 2^49);
n = numel(embedding);
if n < 2
    error('padic:smallEmbedding', 'An embedding needs at least two leaves.');
end
if numel(unique(embedding)) ~= n
    error('padic:coordinateCollision', 'Embedding coordinates must be distinct.');
end

useSharedLookup = ~isempty(sharedLookup);
if useSharedLookup
    collisions = sharedCollisionParameters( ...
        sharedLookup, collisionDepth, collisionZ);
    if ~isequal(sharedLookup.embedding, embedding)
        error('padic:sharedCollisionEmbedding', ...
            'The shared collision lookup belongs to another embedding.');
    end
    collisionDepth = double(collisionDepth);
    collisionZ = double(collisionZ(:));
    a = embedding;
    b = -(2^collisionDepth - 1).*embedding - collisionZ;
    if validate
        if ~isequal(evaluateIntegerFunction(g, embedding, 'g'), a) || ...
                ~isequal(evaluateIntegerFunction(h, embedding, 'h'), b)
            error('padic:sharedCollisionMetadata', ...
                'Family d,z metadata disagrees with its g,h handles.');
        end
    end
else
    a = evaluateIntegerFunction(g, embedding, 'g');
    b = evaluateIntegerFunction(h, embedding, 'h');
end
excludeZeroRegion = isequal(a, embedding);
if any(abs(a) >= 2^24) || any(abs(b) >= 2^24)
    unavailable(['Integer g/h values must have magnitude below 2^24 for ' ...
        'provably exact collision arithmetic.']);
end

if ~useSharedLookup
    collisions = integerCollisionParameters(a, b);
end
best = struct.empty(0, 1);
generatedEmbeddings = zeros(n, 0);
if ~collisions.valid
    if nargout >= 3
        diagnostic = diagnosticResult(collisions, zeros(0, 1), ...
            useSharedLookup);
    end
    return;
end

tree = [];
if ~isempty(collisions.values)
    if useSharedLookup
        tree = sharedCollisionRadixTree( ...
            sharedLookup, collisions.qIds, validate);
    else
        tree = integerCollisionRadixTree(collisions, validate);
    end
end
if isempty(tree) && excludeZeroRegion
    representatives = zeros(0, 1);
elseif isempty(tree)
    representatives = 0;
else
    representatives = integerRepresentativeParameters( ...
        tree, validate, excludeZeroRegion);
end
if nargout >= 3
    diagnostic = diagnosticResult( ...
        collisions, representatives, useSharedLookup);
end

enumerateOnly = isempty(scorePlan);
if enumerateOnly
    generatedEmbeddings = zeros(n, numel(representatives));
end

function result = diagnosticResult(collisions, representatives, usedSharedLookup)
result = struct('collisions', collisions, ...
    'representatives', representatives, ...
    'usedSharedLookup', usedSharedLookup);
end

for index = 1:numel(representatives)
    theta = representatives(index);
    if any(abs(theta .* b) >= 2^49) || any(abs(a) >= 2^49)
        unavailable('A transformed coordinate would exceed the exact integer bound.');
    end
    coordinates = a - theta .* b;
    if any(abs(coordinates) >= 2^50)
        unavailable('A transformed coordinate would exceed the exact integer bound.');
    end
    if numel(unique(coordinates)) ~= n
        error('padic:generatedCollision', ...
            'A representative is an actual collision.');
    end

    if enumerateOnly
        generatedEmbeddings(:, index) = coordinates;
        continue;
    end

    [score, valuations] = caminSokalScoreEmbeddingPrepared(coordinates, scorePlan);
    if score >= currentScore - improvementTolerance
        continue;
    end
    if isempty(best) || score < best.score
        best = candidate(score, false(0, n), false, a, b, theta, coordinates);
    elseif score == best.score
        key = topologyKeyFromValuations(valuations);
        if ~best.keyKnown
            [~, bestValuations] = caminSokalScoreEmbeddingPrepared( ...
                best.transformedCoordinates, scorePlan);
            best.key = topologyKeyFromValuations(bestValuations);
            best.keyKnown = true;
        end
        if lexicographicallyLess(key, best.key)
            best = candidate(score, key, true, a, b, theta, coordinates);
        end
    end
end
end

function result = candidate(score, key, keyKnown, a, b, theta, coordinates)
result = struct('key', key, 'keyKnown', keyKnown, 'score', score, ...
    'gValues', a, 'hValues', b, 'theta', theta, ...
    'transformedCoordinates', coordinates);
end

function values = integerVector(input, name, limit)
if isnumeric(input) || islogical(input)
    values = double(input(:));
    usable = ~isempty(values) && all(isfinite(values)) && ...
        all(values == fix(values)) && all(abs(values) < limit);
else
    symbolic = sym(input(:));
    [values, usable] = exactIntegerDoubles(symbolic);
    usable = usable && ~isempty(values) && all(abs(values) < limit);
end
if ~usable
    unavailable(sprintf('%s is not a safely representable integer vector.', name));
end
end

function values = evaluateIntegerFunction(functionHandle, argument, name)
n = numel(argument);
try
    raw = functionHandle(argument);
    [values, vectorWorked] = numericIntegerOutput(raw, n);
catch
    vectorWorked = false;
end
if vectorWorked
    return;
end
values = zeros(n, 1);
for index = 1:n
    try
        raw = functionHandle(argument(index));
    catch cause
        error('padic:integerFastPathUnavailable', ...
            'The %s handle cannot be evaluated by the integer engine: %s', ...
            name, cause.message);
    end
    [scalar, usable] = numericIntegerOutput(raw, 1);
    if ~usable
        unavailable(sprintf('%s did not return machine integers.', name));
    end
    values(index) = scalar;
end
end

function [values, usable] = numericIntegerOutput(raw, expectedCount)
usable = (isnumeric(raw) || islogical(raw)) && numel(raw) == expectedCount;
if ~usable
    values = zeros(0, 1);
    return;
end
values = double(raw(:));
usable = isreal(values) && all(isfinite(values)) && ...
    all(values == fix(values)) && all(abs(values) < 2^49);
end

function collision = integerCollisionParameters(a, b)
n = numel(a);
degenerate = zeros(0, 2);
for i = 1:n
    for j = i + 1:n
        if a(i) == a(j) && b(i) == b(j)
            degenerate(end + 1, :) = [i, j]; %#ok<AGROW>
        end
    end
end
collision = struct('valid', isempty(degenerate), ...
    'values', zeros(0, 1), 'numerators', zeros(0, 1), ...
    'denominators', zeros(0, 1), 'degeneratePairs', degenerate);
if ~collision.valid
    return;
end

keyToIndex = containers.Map('KeyType', 'char', 'ValueType', 'double');
numerators = zeros(0, 1);
denominators = zeros(0, 1);
for i = 1:n
    for j = i + 1:n
        denominator = b(i) - b(j);
        if denominator == 0
            continue;
        end
        numerator = a(i) - a(j);
        common = gcd(abs(numerator), abs(denominator));
        numerator = numerator / common;
        denominator = denominator / common;
        if denominator < 0
            numerator = -numerator;
            denominator = -denominator;
        end
        if numerator ~= 0 && ...
                integerValuation(numerator) < integerValuation(denominator)
            continue;
        end
        key = sprintf('%.0f/%.0f', numerator, denominator);
        if ~isKey(keyToIndex, key)
            numerators(end + 1, 1) = numerator; %#ok<AGROW>
            denominators(end + 1, 1) = denominator; %#ok<AGROW>
            keyToIndex(key) = numel(numerators);
        end
    end
end

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
collision.numerators = numerators(order);
collision.denominators = denominators(order);
collision.values = collision.numerators ./ collision.denominators;
end

function theta = integerRepresentativeParameters(tree, validate, excludeZeroRegion)
if isscalar(tree.values)
    % With one collision, the only valid parameter region is Z_2 minus the
    % collision point. It contains zero unless the collision itself is zero.
    if excludeZeroRegion && tree.numerators(1) ~= 0
        theta = zeros(0, 1);
    else
        theta = rationalPhi(tree.numerators(1), tree.denominators(1), 0);
    end
    return;
end
theta = zeros(0, 1);
seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
visit(tree.root, 0, excludeZeroRegion);
if validate
    if numel(theta) > 2 * numel(tree.values) - 1
        error('padic:representativeBound', ...
            'The 2m-1 representative bound was exceeded.');
    end
    for i = 1:numel(theta)
        for j = 1:numel(tree.values)
            if tree.denominators(j) == 1 && theta(i) == tree.numerators(j)
                error('padic:collisionRepresentative', ...
                    'Generated an actual collision value.');
            end
        end
    end
end

    function add(value, containsZero)
        if containsZero
            return;
        end
        key = sprintf('%.0f', value);
        if ~isKey(seen, key)
            seen(key) = true;
            theta(end + 1, 1) = value;
        end
    end

    function visit(nodeIndex, ambientDepth, ambientContainsZero)
        node = tree.nodes(nodeIndex);
        if node.isLeaf
            % A leaf represents its ambient ball with the collision point
            % removed. If that ball contains zero and the collision is not
            % zero, this is the unique zero-containing parameter region.
            residualContainsZero = ambientContainsZero && node.numerator ~= 0;
            add(rationalPhi(node.numerator, node.denominator, ambientDepth), ...
                residualContainsZero);
            return;
        end

        if node.epsilon > ambientDepth
            nodeContainsZero = ambientContainsZero && containsZeroBall( ...
                node.numerator, node.denominator, node.epsilon);
            residualContainsZero = ambientContainsZero && ~nodeContainsZero;
            add(rationalPhi(node.numerator, node.denominator, node.epsilon - 1), ...
                residualContainsZero);
        else
            nodeContainsZero = ambientContainsZero;
        end

        childDepth = node.epsilon + 1;
        childContainsZero = false(1, 2);
        for childPosition = 1:2
            child = tree.nodes(node.children(childPosition));
            childContainsZero(childPosition) = nodeContainsZero && ...
                containsZeroBall(child.numerator, child.denominator, childDepth);
        end
        if validate && nnz(childContainsZero) > 1
            error('padic:multipleZeroChildren', ...
                'More than one radix child was marked as containing zero.');
        end
        visit(node.children(1), childDepth, childContainsZero(1));
        visit(node.children(2), childDepth, childContainsZero(2));
    end
end

function result = containsZeroBall(numerator, denominator, depth)
% B = numerator/denominator + 2^depth Z_2 contains zero iff v_2(q)>=depth.
if mod(denominator, 2) ~= 1
    error('padic:notIntegral', ...
        'Zero containment requires an odd reduced denominator.');
end
result = depth == 0 || numerator == 0 || integerValuation(numerator) >= depth;
end

function value = rationalPhi(numerator, denominator, depth)
if depth < 0 || depth + 1 > 50
    unavailable('A collision depth exceeds the exact integer fast-path bound.');
end
modulus = 2^(depth + 1);
residue = rationalResidue(numerator, denominator, depth + 1);
value = mod(residue + 2^depth, modulus);
end

function residue = rationalResidue(numerator, denominator, digits)
if digits == 0
    residue = 0;
    return;
end
if digits > 50
    unavailable('A collision depth exceeds the exact integer fast-path bound.');
end
modulus = 2^digits;
inverse = inverseOddPowerOfTwo(denominator, digits);
residue = multiplyModulo(mod(numerator, modulus), inverse, modulus);
end

function inverse = inverseOddPowerOfTwo(value, digits)
inverse = 1;
for bit = 2:digits
    modulus = 2^bit;
    if multiplyModulo(mod(value, modulus), inverse, modulus) ~= 1
        inverse = inverse + 2^(bit - 1);
    end
end
end

function result = multiplyModulo(left, right, modulus)
left = mod(left, modulus);
right = mod(right, modulus);
result = 0;
while right > 0
    if mod(right, 2) == 1
        result = mod(result + left, modulus);
    end
    left = mod(2 * left, modulus);
    right = floor(right / 2);
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

function unavailable(message)
error('padic:integerFastPathUnavailable', '%s', message);
end
