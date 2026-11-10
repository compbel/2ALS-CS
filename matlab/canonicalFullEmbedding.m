function coordinates = canonicalFullEmbedding(tree, validate)
%CANONICALFULLEMBEDDING Canonical exact 2-adic embedding of a digraph tree.

if nargin < 2
    validate = true;
end
if validate
    tree = canonicalizeTopology(tree);
end
taxon = tree.Nodes.Taxon;
n = nnz(taxon > 0);
root = find(indegree(tree) == 0, 1);
useNumeric = n <= 50;
if useNumeric
    coordinates = zeros(n, 1);
else
    coordinates = sym(zeros(n, 1));
end
minimumTaxon = inf(numnodes(tree), 1);
order = flip(toposort(tree));
for node = reshape(order, 1, [])
    if taxon(node) > 0
        minimumTaxon(node) = taxon(node);
    else
        minimumTaxon(node) = min(minimumTaxon(successors(tree, node)));
    end
end
if useNumeric
    visit(root, 0, 0);
else
    visit(root, sym(0), 0);
end
if validate && ~isequal(topologyKey(topologyFromCoordinates(coordinates)), topologyKey(tree))
    error('padic:embeddingRoundTrip', 'Canonical embedding failed to reconstruct topology.');
end

    function visit(node, value, depth)
        children = successors(tree, node).';
        if isempty(children)
            coordinates(tree.Nodes.Taxon(node)) = value;
            return;
        end
        minima = minimumTaxon(children);
        [~, order] = sort(minima);
        children = children(order);
        visit(children(1), value, depth + 1);
        if useNumeric
            branchValue = value + 2^depth;
        else
            branchValue = value + sym(2)^depth;
        end
        visit(children(2), branchValue, depth + 1);
    end
end
