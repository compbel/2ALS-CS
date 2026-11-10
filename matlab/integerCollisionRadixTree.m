function tree = integerCollisionRadixTree(collision, validate)
%INTEGERCOLLISIONRADIXTREE Build the exact compressed 2-adic collision trie.
%   Used both for individual integer families and as the one-time global
%   builder for shared Camin-Sokal collision values.

if nargin < 2, validate = false; end
numerators = collision.numerators(:);
denominators = collision.denominators(:);
if isempty(numerators)
    error('padic:emptyTree', 'A collision radix tree requires at least one value.');
end
if any(mod(denominators, 2) ~= 1)
    error('padic:notIntegral', 'Every radix-tree value must lie in Z_2.');
end
nodeTemplate = struct('isLeaf', false, 'start', 0, 'stop', 0, ...
    'representative', 0, 'epsilon', NaN, 'children', zeros(1, 0), ...
    'value', 0, 'numerator', NaN, 'denominator', NaN);
nodes = repmat(nodeTemplate, 0, 1);
[root, order, internalCount] = buildNode(1:numel(numerators), 1);
orderedNumerators = numerators(order);
orderedDenominators = denominators(order);
tree = struct('values', orderedNumerators ./ orderedDenominators, ...
    'numerators', orderedNumerators, 'denominators', orderedDenominators, ...
    'qIds', reshape(order, [], 1), ...
    'root', root, 'nodes', nodes, 'internalCount', internalCount, ...
    'vertexCount', numel(order) + internalCount, 'arithmeticMode', 'integer');
if validate
    validateNode(root);
    if internalCount > numel(numerators) - 1
        error('padic:treeSize', 'Compressed tree has too many internal nodes.');
    end
end

    function [nodeIndex, ordered, count] = buildNode(group, offset)
        if isscalar(group)
            valueIndex = group(1);
            node = struct('isLeaf', true, 'start', offset, ...
                'stop', offset + 1, ...
                'representative', numerators(valueIndex) / denominators(valueIndex), ...
                'epsilon', NaN, 'children', zeros(1, 0), ...
                'value', numerators(valueIndex) / denominators(valueIndex), ...
                'numerator', numerators(valueIndex), ...
                'denominator', denominators(valueIndex));
            nodes(end + 1, 1) = node;
            nodeIndex = numel(nodes);
            ordered = group;
            count = 0;
            return;
        end
        epsilon = inf;
        pivot = group(1);
        for position = 2:numel(group)
            epsilon = min(epsilon, rationalDifferenceValuation( ...
                numerators(group(position)), denominators(group(position)), ...
                numerators(pivot), denominators(pivot)));
        end
        bits = zeros(size(group));
        for position = 1:numel(group)
            bits(position) = rationalBit(numerators(group(position)), ...
                denominators(group(position)), epsilon);
        end
        zero = group(bits == 0);
        one = group(bits == 1);
        if isempty(zero) || isempty(one)
            error('padic:radixSplit', 'First differing digit failed to split a node.');
        end
        [left, leftOrder, leftCount] = buildNode(zero, offset);
        [right, rightOrder, rightCount] = buildNode(one, offset + numel(zero));
        ordered = [leftOrder, rightOrder];
        node = struct('isLeaf', false, 'start', offset, ...
            'stop', offset + numel(group), ...
            'representative', numerators(pivot) / denominators(pivot), ...
            'epsilon', epsilon, 'children', [left, right], 'value', 0, ...
            'numerator', numerators(pivot), 'denominator', denominators(pivot));
        nodes(end + 1, 1) = node;
        nodeIndex = numel(nodes);
        count = 1 + leftCount + rightCount;
    end

    function validateNode(nodeIndex)
        node = nodes(nodeIndex);
        if node.isLeaf
            if node.stop ~= node.start + 1 || ...
                    node.numerator ~= orderedNumerators(node.start) || ...
                    node.denominator ~= orderedDenominators(node.start)
                error('padic:leafRange', 'Invalid radix leaf range.');
            end
            return;
        end
        expected = rationalResidue(orderedNumerators(node.start), ...
            orderedDenominators(node.start), node.epsilon);
        for position = node.start + 1:node.stop - 1
            actual = rationalResidue(orderedNumerators(position), ...
                orderedDenominators(position), node.epsilon);
            if actual ~= expected
                error('padic:prefixMismatch', ...
                    'Descendants disagree below epsilon.');
            end
        end
        left = nodes(node.children(1));
        right = nodes(node.children(2));
        childBits = [rationalBit(orderedNumerators(left.start), ...
            orderedDenominators(left.start), node.epsilon), ...
            rationalBit(orderedNumerators(right.start), ...
            orderedDenominators(right.start), node.epsilon)];
        if ~isequal(sort(childBits), [0, 1])
            error('padic:childBits', 'Children do not separate at epsilon.');
        end
        validateNode(node.children(1));
        validateNode(node.children(2));
    end
end

function value = rationalDifferenceValuation(n1, d1, n2, d2)
differenceNumerator = n1 * d2 - n2 * d1;
value = integerValuation(differenceNumerator);
end

function bit = rationalBit(numerator, denominator, depth)
residue = rationalResidue(numerator, denominator, depth + 1);
bit = mod(floor(residue / 2^depth), 2);
end

function residue = rationalResidue(numerator, denominator, digits)
if digits == 0
    residue = 0;
    return;
elseif digits > 50
    error('padic:integerFastPathUnavailable', ...
        'A collision depth exceeds the exact integer fast-path bound.');
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
