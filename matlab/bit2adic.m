function bit = bit2adic(q, k)
%BIT2ADIC Return zero-based 2-adic digit k of q in Z_2.

if ~isscalar(k) || k < 0 || k ~= floor(k)
    error('padic:badBit', 'k must be a nonnegative integer.');
end
bit = double(mod(floor(residue2adic(q, k + 1) / sym(2)^k), sym(2)));
end
