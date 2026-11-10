function text = topologyKeyText(treeOrKey)
%TOPOLOGYKEYTEXT Compact non-Newick serialization for maps and ordering.

if isa(treeOrKey, 'digraph')
    key = topologyKey(treeOrKey);
else
    key = treeOrKey;
    if ~islogical(key) || ~ismatrix(key)
        error('padic:keyType', 'A topology key must be a logical matrix.');
    end
end
n = size(key, 2);
bits = char('0' + reshape(key.', 1, []));
text = sprintf('%d:%s', n, bits);
end
