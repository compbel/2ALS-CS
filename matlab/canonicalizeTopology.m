function canonical = canonicalizeTopology(tree, validate)
%CANONICALIZETOPOLOGY Rebuild a rooted binary digraph deterministically.
%   Leaves become vertices 1..n in taxon order. Internal vertices are
%   numbered in postorder, with children ordered by minimum descendant taxon.

if nargin < 2
    validate = true;
end
if validate
    information = validateRootedDigraph(tree);
    rootNode = information.root;
    n = information.n;
else
    rootNode = find(indegree(tree) == 0, 1);
    n = nnz(tree.Nodes.Taxon > 0);
end
sources = zeros(2*n - 2, 1);
targets = zeros(2*n - 2, 1);
edgeIndex = 0;
nextInternal = n;
[root, ~] = rebuild(rootNode);
canonical = digraph(sources, targets, [], 2*n - 1);
canonical.Nodes.Taxon = [transpose(1:n); zeros(n - 1, 1)];
if root ~= 2*n - 1
    error('padic:canonicalRoot', 'Canonical postorder numbering failed.');
end

    function [newNode, minimumTaxon] = rebuild(oldNode)
        children = successors(tree, oldNode).';
        if isempty(children)
            newNode = tree.Nodes.Taxon(oldNode);
            minimumTaxon = newNode;
            return;
        end
        childNodes = zeros(1, 2);
        childMinimum = zeros(1, 2);
        for childIndex = 1:2
            [childNodes(childIndex), childMinimum(childIndex)] = rebuild(children(childIndex));
        end
        [childMinimum, order] = sort(childMinimum);
        childNodes = childNodes(order);
        nextInternal = nextInternal + 1;
        newNode = nextInternal;
        minimumTaxon = childMinimum(1);
        for childIndex = 1:2
            edgeIndex = edgeIndex + 1;
            sources(edgeIndex) = newNode;
            targets(edgeIndex) = childNodes(childIndex);
        end
    end
end
