function key = topologyKey(tree, validate)
%TOPOLOGYKEY Canonical logical clade matrix for a rooted labeled topology.
%   Each row represents one nonroot internal clade and each column one taxon.
%   Rows are sorted lexicographically. The n=2 tree has a 0-by-2 key.

if nargin < 2
    validate = false;
end
if validate
    validateRootedDigraph(tree);
end
taxon = tree.Nodes.Taxon;
n = nnz(taxon > 0);
root = find(indegree(tree) == 0, 1);
internal = find(outdegree(tree) > 0);
nodes = setdiff(internal, root, 'stable');
clades = false(numnodes(tree), n);
order = flip(toposort(tree));
for node = reshape(order, 1, [])
    if taxon(node) > 0
        clades(node, taxon(node)) = true;
    else
        children = successors(tree, node);
        clades(node, :) = any(clades(children, :), 1);
    end
end
key = clades(nodes, :);
if ~isempty(key)
    key = logical(sortrows(double(key)));
end
end
