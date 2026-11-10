function score = caminSokalScoreKey(key, alignment)
%CAMINSOKALSCOREKEY Irreversible binary score directly from internal clades.
%   For each character, begin with one gain per 1-leaf. Every all-1 internal
%   clade merges two gain components and reduces the score by one. KEY stores
%   all nonroot internal clades; the full root clade is handled separately.

alignment = double(alignment);
zeroMask = double(alignment == 0);
if isempty(key)
    containedInternal = zeros(0, size(alignment, 2));
else
    containedInternal = double(key) * zeroMask == 0;
end
rootAllOne = all(alignment == 1, 1);
siteScores = sum(alignment, 1) - sum(containedInternal, 1) - rootAllOne;
score = sum(siteScores);
end
