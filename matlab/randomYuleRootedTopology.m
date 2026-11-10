function tree = randomYuleRootedTopology(n)
%RANDOMYULEROOTEDTOPOLOGY Generate a Yule-Harding rooted binary topology.
%   Starting from a two-leaf tree, each new leaf is added by choosing one
%   extant leaf uniformly and splitting its pendant edge. Taxon labels are
%   randomly permuted after the unlabeled Yule shape has been generated.

if ~isscalar(n) || n < 2 || n ~= floor(n)
    error('padic:taxonCount', 'n must be an integer at least 2.');
end

tree = digraph([3; 3], [1; 2], [], 3);
tree.Nodes.Taxon = [1; 2; 0];
for taxon = 3:n
    leaves = find(outdegree(tree) == 0);
    chosenLeaf = leaves(randi(numel(leaves)));
    parent = predecessors(tree, chosenLeaf);
    oldVertexCount = numnodes(tree);
    newLeaf = oldVertexCount + 1;
    newInternal = oldVertexCount + 2;
    tree = addnode(tree, 2);
    tree.Nodes.Taxon(newLeaf) = taxon;
    tree.Nodes.Taxon(newInternal) = 0;
    tree = rmedge(tree, parent, chosenLeaf);
    tree = addedge(tree, [parent; newInternal; newInternal], ...
        [newInternal; chosenLeaf; newLeaf]);
end

permutation = randperm(n);
leafNodes = find(outdegree(tree) == 0);
for index = 1:numel(leafNodes)
    oldTaxon = tree.Nodes.Taxon(leafNodes(index));
    tree.Nodes.Taxon(leafNodes(index)) = permutation(oldTaxon);
end
tree = canonicalizeTopology(tree);
end
