function data = topologyDataFromCoordinates(coordinates)
%TOPOLOGYDATAFROMCOORDINATES Recursive exact 2-adic topology construction.
%   Coordinates are split at their first differing least-significant bit.

coordinates = coordinates(:);
if isempty(coordinates)
    error('padic:emptyEmbedding', 'An embedding must contain a coordinate.');
end
n = numel(coordinates);
if n < 2
    error('padic:smallEmbedding', 'An embedding requires at least two coordinates.');
end

[numericCoordinates, useNumeric] = exactNumericCoordinates(coordinates);
if useNumeric
    if numel(unique(numericCoordinates)) ~= n
        error('padic:coordinateCollision', ...
            'Embedding coordinates must be distinct.');
    end
    symbolicCoordinates = [];
else
    symbolicCoordinates = sym(coordinates);
    assertDistinctSymbolic(symbolicCoordinates);
end

sources = zeros(2*n - 2, 1);
targets = zeros(2*n - 2, 1);
internalClades = false(n - 1, n);
internalCount = 0;
edgeCount = 0;
[root, ~] = build(1:n);
if internalCount ~= n - 1 || edgeCount ~= 2*n - 2 || root ~= 2*n - 1
    error('padic:topologySize', ...
        'Recursive 2-adic reconstruction produced an invalid tree size.');
end

key = internalClades(1:max(n - 2, 0), :);
if ~isempty(key)
    key = logical(sortrows(double(key)));
end
data = struct('sources', sources, 'targets', targets, ...
    'taxon', [transpose(1:n); zeros(n - 1, 1)], 'key', key);

    function [node, clade] = build(labels)
        if isscalar(labels)
            node = labels(1);
            clade = false(1, n);
            clade(node) = true;
            return;
        end

        epsilon = inf;
        if useNumeric
            pivot = numericCoordinates(labels(1));
            for index = 2:numel(labels)
                epsilon = min(epsilon, valuationInteger( ...
                    numericCoordinates(labels(index)) - pivot));
            end
            bits = mod(floor(numericCoordinates(labels) / 2^epsilon), 2);
        else
            pivot = symbolicCoordinates(labels(1));
            for index = 2:numel(labels)
                epsilon = min(epsilon, v2Exact( ...
                    symbolicCoordinates(labels(index)) - pivot));
            end
            bits = zeros(size(labels));
            for index = 1:numel(labels)
                bits(index) = bit2adic( ...
                    symbolicCoordinates(labels(index)), epsilon);
            end
        end

        zeroLabels = labels(bits == 0);
        oneLabels = labels(bits == 1);
        if isempty(zeroLabels) || isempty(oneLabels)
            error('padic:radixSplit', ...
                'First differing digit failed to split a node.');
        end
        [zeroNode, zeroClade] = build(zeroLabels);
        [oneNode, oneClade] = build(oneLabels);

        internalCount = internalCount + 1;
        node = n + internalCount;
        clade = zeroClade | oneClade;
        internalClades(internalCount, :) = clade;
        edgeCount = edgeCount + 1;
        sources(edgeCount) = node;
        targets(edgeCount) = zeroNode;
        edgeCount = edgeCount + 1;
        sources(edgeCount) = node;
        targets(edgeCount) = oneNode;
    end
end

function [numeric, usable] = exactNumericCoordinates(coordinates)
if isnumeric(coordinates) || islogical(coordinates)
    numeric = double(coordinates);
    usable = isreal(numeric) && all(isfinite(numeric)) && ...
        all(numeric == fix(numeric)) && all(abs(numeric) < 2^50);
else
    [numeric, usable] = exactIntegerDoubles(sym(coordinates));
    usable = usable && all(abs(numeric) < 2^50);
end
if ~usable
    numeric = zeros(0, 1);
end
end

function assertDistinctSymbolic(coordinates)
for left = 1:numel(coordinates) - 1
    for right = left + 1:numel(coordinates)
        if isAlways(coordinates(left) == coordinates(right))
            error('padic:coordinateCollision', ...
                'Embedding coordinates must be distinct.');
        end
    end
end
end

function valuation = valuationInteger(value)
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
