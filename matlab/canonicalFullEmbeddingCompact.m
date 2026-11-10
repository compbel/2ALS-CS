function coordinates = canonicalFullEmbeddingCompact(data)
%CANONICALFULLEMBEDDINGCOMPACT Canonical embedding from compact tree arrays.
%   DATA is the structure returned by topologyDataFromCoordinates or
%   compactTopologyFromDigraph. The calculation uses only numeric arrays and
%   does not construct or traverse a MATLAB digraph.

required = {'sources', 'targets', 'taxon', 'key'};
if ~isstruct(data) || ~isscalar(data) || ~all(isfield(data, required))
    error('padic:compactTopology', ...
        'A compact topology must contain sources, targets, taxon, and key.');
end
taxon = data.taxon(:);
vertexCount = numel(taxon);
n = nnz(taxon > 0);
sources = data.sources(:);
targets = data.targets(:);
if n < 2 || vertexCount ~= 2*n - 1 || ...
        any(~isfinite(taxon)) || any(taxon ~= floor(taxon)) || any(taxon < 0) || ...
        ~isequal(sort(taxon(taxon > 0)).', 1:n) || ...
        numel(sources) ~= numel(targets) || numel(sources) ~= vertexCount - 1 || ...
        any(~isfinite(sources)) || any(~isfinite(targets)) || ...
        any(sources ~= floor(sources)) || any(targets ~= floor(targets)) || ...
        any(sources < 1 | sources > vertexCount) || ...
        any(targets < 1 | targets > vertexCount) || any(sources == targets)
    error('padic:compactTopology', 'Compact topology edge arrays are invalid.');
end

children = zeros(vertexCount, 2);
childCount = zeros(vertexCount, 1);
parent = zeros(vertexCount, 1);
for edgeIndex = 1:numel(sources)
    source = sources(edgeIndex);
    target = targets(edgeIndex);
    childCount(source) = childCount(source) + 1;
    if childCount(source) > 2 || parent(target) ~= 0
        error('padic:compactTopology', 'Compact topology is not a rooted binary tree.');
    end
    children(source, childCount(source)) = target;
    parent(target) = source;
end
root = find(parent == 0);
if numel(root) ~= 1 || any(childCount(taxon == 0) ~= 2) || ...
        any(childCount(taxon > 0) ~= 0)
    error('padic:compactTopology', 'Compact topology is not a rooted binary tree.');
end
visited = false(vertexCount, 1);
stack = root;
while ~isempty(stack)
    stackNode = stack(end);
    stack(end) = [];
    if visited(stackNode)
        error('padic:compactTopology', 'Compact topology contains a cycle.');
    end
    visited(stackNode) = true;
    reachableChildren = children(stackNode, :);
    stack = [stack, reachableChildren(reachableChildren > 0)]; %#ok<AGROW>
end
if ~all(visited)
    error('padic:compactTopology', 'Compact topology is disconnected.');
end

minimumTaxon = inf(vertexCount, 1);
computeMinimum(root);
if n <= 50
    coordinates = zeros(n, 1);
    visit(root, 0, 0);
else
    coordinates = sym(zeros(n, 1));
    visit(root, sym(0), 0);
end

    function minimum = computeMinimum(treeNode)
        if taxon(treeNode) > 0
            minimum = taxon(treeNode);
        else
            minimum = min(computeMinimum(children(treeNode, 1)), ...
                computeMinimum(children(treeNode, 2)));
        end
        minimumTaxon(treeNode) = minimum;
    end

    function visit(treeNode, value, depth)
        if taxon(treeNode) > 0
            coordinates(taxon(treeNode)) = value;
            return;
        end
        branchChildren = children(treeNode, :);
        [~, order] = sort(minimumTaxon(branchChildren));
        branchChildren = branchChildren(order);
        visit(branchChildren(1), value, depth + 1);
        if n <= 50
            branchValue = value + 2^depth;
        else
            branchValue = value + sym(2)^depth;
        end
        visit(branchChildren(2), branchValue, depth + 1);
    end
end
