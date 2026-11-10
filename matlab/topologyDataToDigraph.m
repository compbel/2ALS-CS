function topology = topologyDataToDigraph(data)
%TOPOLOGYDATATODIGRAPH Materialize a compact topology for returned output.

nodeTable = table(data.taxon(:), 'VariableNames', {'Taxon'});
topology = digraph(data.sources(:), data.targets(:), [], nodeTable);
end
