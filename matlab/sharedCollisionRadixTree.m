function tree = sharedCollisionRadixTree(lookup, selectedQIds, validate)
%SHAREDCOLLISIONRADIXTREE Extract one family's trie from the global trie.
%   Leaves not selected by qID are removed, unary nodes are suppressed, and
%   the branching depths of retained nodes are inherited unchanged. The
%   global trie and lookup remain read-only.

if nargin < 3, validate = false; end
if ~isstruct(lookup) || ~isscalar(lookup) || ...
        ~isfield(lookup, 'globalRadixTree') || isempty(lookup.globalRadixTree)
    error('padic:sharedCollisionLookup', ...
        'The lookup does not contain a global collision radix tree.');
end
selectedQIds = double(selectedQIds(:));
qCount = numel(lookup.numerators);
if isempty(selectedQIds) || any(~isfinite(selectedQIds)) || ...
        any(selectedQIds ~= fix(selectedQIds)) || ...
        any(selectedQIds < 1) || any(selectedQIds > qCount) || ...
        numel(unique(selectedQIds)) ~= numel(selectedQIds)
    error('padic:sharedCollisionIds', ...
        'selectedQIds must be distinct valid qIDs.');
end

globalTree = lookup.globalRadixTree;
selectedById = false(qCount, 1);
selectedById(selectedQIds) = true;
selectedInOrder = selectedById(globalTree.qIds);
prefixCount = [0; cumsum(selectedInOrder)];
nodeTemplate = globalTree.nodes(1);
nodes = repmat(nodeTemplate, 0, 1);
leafCount = numel(selectedQIds);
orderedNumerators = zeros(leafCount, 1);
orderedDenominators = zeros(leafCount, 1);
orderedQIds = zeros(leafCount, 1);
[root, retainedLeaves, internalCount] = retain(globalTree.root, 1);
if retainedLeaves ~= leafCount || root == 0
    error('padic:sharedCollisionTree', ...
        'The extracted collision tree lost selected qIDs.');
end
tree = struct('values', orderedNumerators ./ orderedDenominators, ...
    'numerators', orderedNumerators, 'denominators', orderedDenominators, ...
    'qIds', orderedQIds, 'root', root, 'nodes', nodes, ...
    'internalCount', internalCount, ...
    'vertexCount', leafCount + internalCount, ...
    'arithmeticMode', 'integer_shared');
if validate
    direct = integerCollisionRadixTree(struct( ...
        'numerators', lookup.numerators(selectedQIds), ...
        'denominators', lookup.denominators(selectedQIds)), true);
    assertEquivalentStructure(tree, direct);
end

    function count = selectedCount(node)
        count = prefixCount(node.stop) - prefixCount(node.start);
    end

    function [nodeIndex, leaves, internals] = retain(globalIndex, offset)
        globalNode = globalTree.nodes(globalIndex);
        if selectedCount(globalNode) == 0
            nodeIndex = 0;
            leaves = 0;
            internals = 0;
            return;
        end
        if globalNode.isLeaf
            qId = globalTree.qIds(globalNode.start);
            node = globalNode;
            node.start = offset;
            node.stop = offset + 1;
            nodes(end + 1, 1) = node;
            nodeIndex = numel(nodes);
            orderedNumerators(offset) = node.numerator;
            orderedDenominators(offset) = node.denominator;
            orderedQIds(offset) = qId;
            leaves = 1;
            internals = 0;
            return;
        end

        [left, leftLeaves, leftInternals] = retain( ...
            globalNode.children(1), offset);
        [right, rightLeaves, rightInternals] = retain( ...
            globalNode.children(2), offset + leftLeaves);
        leaves = leftLeaves + rightLeaves;
        internals = leftInternals + rightInternals;
        if left == 0
            nodeIndex = right;
            return;
        elseif right == 0
            nodeIndex = left;
            return;
        end
        node = globalNode;
        node.start = offset;
        node.stop = offset + leaves;
        node.children = [left, right];
        nodes(end + 1, 1) = node;
        nodeIndex = numel(nodes);
        internals = internals + 1;
    end
end

function assertEquivalentStructure(extracted, direct)
if ~isequal(extracted.numerators, direct.numerators) || ...
        ~isequal(extracted.denominators, direct.denominators) || ...
        extracted.internalCount ~= direct.internalCount
    error('padic:sharedCollisionTree', ...
        'The extracted and directly built collision trees differ.');
end
compareNodes(extracted.root, direct.root);

    function compareNodes(extractedIndex, directIndex)
        left = extracted.nodes(extractedIndex);
        right = direct.nodes(directIndex);
        if left.isLeaf ~= right.isLeaf || left.start ~= right.start || ...
                left.stop ~= right.stop || ~isequaln(left.epsilon, right.epsilon)
            error('padic:sharedCollisionTree', ...
                'The extracted and directly built collision nodes differ.');
        end
        if ~left.isLeaf
            compareNodes(left.children(1), right.children(1));
            compareNodes(left.children(2), right.children(2));
        end
    end
end
