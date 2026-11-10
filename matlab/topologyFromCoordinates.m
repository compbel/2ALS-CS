function tree = topologyFromCoordinates(coordinates)
%TOPOLOGYFROMCOORDINATES Exact 2-adic hierarchical tree as a digraph.

data = topologyDataFromCoordinates(coordinates);
tree = digraph(data.sources, data.targets, [], numel(data.taxon));
tree.Nodes.Taxon = data.taxon;
end
