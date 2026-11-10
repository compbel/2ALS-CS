function key = topologyKeyFromValuations(valuations)
%TOPOLOGYKEYFROMVALUATIONS Canonical clade key without constructing a tree.
%   VALUATIONS is the symmetric v2 distance-level matrix with Inf diagonal.

if ~isnumeric(valuations) || size(valuations, 1) ~= size(valuations, 2) || ...
        any(any(valuations ~= transpose(valuations)))
    error('padic:valuationMatrix', ...
        'valuations must be a symmetric square numeric matrix.');
end
n = size(valuations, 1);
clades = false(n*(n - 1)/2, n);
count = 0;
for left = 1:n - 1
    for right = left + 1:n
        count = count + 1;
        clades(count, :) = valuations(left, :) >= valuations(left, right);
    end
end
clades = unique(clades(1:count, :), 'rows');
clades(sum(clades, 2) == n, :) = [];
key = logical(sortrows(double(clades)));
if size(key, 1) ~= max(n - 2, 0)
    error('padic:valuationTopology', ...
        'Valuation matrix did not induce a full binary topology.');
end
end
