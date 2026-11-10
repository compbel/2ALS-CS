function tree = randomRootedTopology(n)
%RANDOMROOTEDTOPOLOGY Generate a random rooted binary labeled digraph.

if ~isscalar(n) || n < 2 || n ~= floor(n)
    error('padic:taxonCount', 'n must be an integer at least 2.');
end
% Start with leaves 1,2 and root 3. New internal vertices and leaves receive
% temporary node IDs; Taxon stores the biological leaf label.
tree = digraph([3; 3], [1; 2], [], 3);
tree.Nodes.Taxon = [1; 2; 0];
for taxon = 3:n
    oldVertexCount = numnodes(tree);
    oldRoot = find(indegree(tree) == 0);
    edges = tree.Edges.EndNodes;
    location = randi(size(edges, 1) + 1) - 1;
    newLeaf = oldVertexCount + 1;
    newInternal = oldVertexCount + 2;
    tree = addnode(tree, 2);
    tree.Nodes.Taxon(newLeaf) = taxon;
    tree.Nodes.Taxon(newInternal) = 0;
    if location == 0
        tree = addedge(tree, [newInternal; newInternal], [oldRoot; newLeaf]);
    else
        parent = edges(location, 1);
        child = edges(location, 2);
        tree = rmedge(tree, parent, child);
        tree = addedge(tree, [parent; newInternal; newInternal], ...
            [newInternal; child; newLeaf]);
    end
end
permutation = randperm(n);
leafNodes = find(outdegree(tree) == 0);
for index = 1:numel(leafNodes)
    oldTaxon = tree.Nodes.Taxon(leafNodes(index));
    tree.Nodes.Taxon(leafNodes(index)) = permutation(oldTaxon);
end
tree = canonicalizeTopology(tree);
end
