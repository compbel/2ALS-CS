function score = caminSokalScoreArrays(children, root, taxon, alignment)
%CAMINSOKALSCOREARRAYS Irreversible parsimony on compact tree arrays.

alignment = double(alignment);
siteCount = size(alignment, 2);
vertexCount = numel(taxon);
zeroCost = inf(vertexCount, siteCount);
oneCost = inf(vertexCount, siteCount);
preorder = zeros(1, vertexCount);
count = 0;
stack = root;
while ~isempty(stack)
    node = stack(end);
    stack(end) = [];
    count = count + 1;
    preorder(count) = node;
    child = children(node, :);
    stack = [stack, child(child > 0)]; %#ok<AGROW>
end
for node = fliplr(preorder(1:count))
    if taxon(node) > 0
        states = alignment(taxon(node), :);
        zeroCost(node, :) = 0;
        zeroCost(node, states == 1) = inf;
        oneCost(node, :) = 0;
        oneCost(node, states == 0) = inf;
    else
        child = children(node, :);
        zeroCost(node, :) = sum(min(zeroCost(child, :), ...
            oneCost(child, :) + 1), 1);
        oneCost(node, :) = sum(oneCost(child, :), 1);
    end
end
score = sum(min(zeroCost(root, :), oneCost(root, :) + 1));
end
