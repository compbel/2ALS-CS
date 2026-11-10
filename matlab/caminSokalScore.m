function score = caminSokalScore(tree, alignment)
%CAMINSOKALSCORE Irreversible parsimony by one bottom-up tree traversal.
%   All alignment columns are evaluated together. The input tree is assumed
%   to satisfy the rooted-digraph contract enforced by the search boundary.

taxon = tree.Nodes.Taxon;
edges = tree.Edges.EndNodes;
children = zeros(numnodes(tree), 2);
counts = zeros(numnodes(tree), 1);
parent = zeros(numnodes(tree), 1);
for edgeIndex = 1:size(edges, 1)
    source = edges(edgeIndex, 1);
    target = edges(edgeIndex, 2);
    counts(source) = counts(source) + 1;
    children(source, counts(source)) = target;
    parent(target) = source;
end
root = find(parent == 0, 1);
score = caminSokalScoreArrays(children, root, taxon, alignment);
end
