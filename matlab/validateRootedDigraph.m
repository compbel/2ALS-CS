function information = validateRootedDigraph(tree)
%VALIDATEROOTEDDIGRAPH Validate the rooted-tree digraph data model.
%   Edges point from parent to child. Leaves have positive entries in
%   Nodes.Taxon; internal vertices have Taxon zero.

if ~isa(tree, 'digraph')
    error('padic:treeType', 'A tree must be a MATLAB digraph object.');
end
vertexCount = numnodes(tree);
if vertexCount < 3 || ~isdag(tree)
    error('padic:treeShape', 'A rooted binary tree must be a nontrivial DAG.');
end
if ~ismember('Taxon', tree.Nodes.Properties.VariableNames)
    error('padic:taxonMetadata', 'tree.Nodes must contain a numeric Taxon column.');
end
taxon = tree.Nodes.Taxon;
if ~isnumeric(taxon) || numel(taxon) ~= vertexCount || any(~isfinite(taxon)) || ...
        any(taxon ~= floor(taxon)) || any(taxon < 0)
    error('padic:taxonMetadata', 'Nodes.Taxon must contain nonnegative integers.');
end
indegrees = indegree(tree);
outdegrees = outdegree(tree);
roots = find(indegrees == 0);
if numel(roots) ~= 1 || any(indegrees(setdiff(1:vertexCount, roots)) ~= 1)
    error('padic:treeParents', 'The tree must have one root and one parent per nonroot vertex.');
end
leaves = find(outdegrees == 0);
internal = find(outdegrees ~= 0);
if any(outdegrees(internal) ~= 2) || any(taxon(internal) ~= 0)
    error('padic:binaryTree', 'Every internal vertex must have two children and Taxon zero.');
end
leafTaxa = sort(taxon(leaves).');
n = numel(leaves);
if ~isequal(leafTaxa, 1:n) || vertexCount ~= 2*n - 1 || numedges(tree) ~= vertexCount - 1
    error('padic:leafLabels', 'Leaf taxa must be exactly 1..n in a full rooted binary tree.');
end
reachable = bfsearch(tree, roots(1));
if numel(reachable) ~= vertexCount
    error('padic:disconnectedTree', 'Every vertex must be reachable from the root.');
end
information = struct('root', roots(1), 'leaves', leaves(:), ...
    'internal', internal(:), 'n', n, 'vertexCount', vertexCount);
end
