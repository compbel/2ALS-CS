function r = residue2adic(q, K)
%RESIDUE2ADIC Return q modulo 2^K for an exact q in Z_2.

if ~isscalar(K) || K < 0 || K ~= floor(K)
    error('padic:badDepth', 'K must be a nonnegative integer.');
end
q = sym(q);
[u, v] = numden(q);
uNumeric = double(u);
vNumeric = double(v);
if isfinite(vNumeric) && abs(vNumeric) < 2^50
    denominatorOdd = mod(vNumeric, 2) == 1;
else
    denominatorOdd = isAlways(mod(v, sym(2)) == 1);
end
if ~denominatorOdd
    error('padic:notIntegral', 'The rational parameter is not in Z_2.');
end
if K == 0
    r = sym(0);
    return;
end
if K <= 24 && isfinite(uNumeric) && isfinite(vNumeric) && ...
        abs(uNumeric) < 2^50 && abs(vNumeric) < 2^50
    modulusNumeric = 2^K;
    inverseNumeric = inverseModOddNumeric(vNumeric, modulusNumeric);
    r = sym(mod(mod(uNumeric, modulusNumeric) * inverseNumeric, modulusNumeric));
    return;
end
modulus = sym(2)^K;
inverse = inverseModOdd(v, modulus);
r = mod(u * inverse, modulus);
end

function inverse = inverseModOddNumeric(a, modulus)
oldR = mod(a, modulus);
r = modulus;
oldS = 1;
s = 0;
while r ~= 0
    quotient = floor(oldR / r);
    [oldR, r] = deal(r, oldR - quotient * r);
    [oldS, s] = deal(s, oldS - quotient * s);
end
if oldR ~= 1
    error('padic:noInverse', 'Denominator is not invertible modulo 2^K.');
end
inverse = mod(oldS, modulus);
end

function inverse = inverseModOdd(a, modulus)
% Extended Euclid, written explicitly to keep every operation exact.
oldR = mod(a, modulus);
r = modulus;
oldS = sym(1);
s = sym(0);
while ~isAlways(r == 0)
    quotient = floor(oldR / r);
    nextR = oldR - quotient * r;
    oldR = r;
    r = nextR;
    nextS = oldS - quotient * s;
    oldS = s;
    s = nextS;
end
if ~isAlways(oldR == 1)
    error('padic:noInverse', 'Denominator is not invertible modulo 2^K.');
end
inverse = mod(oldS, modulus);
end
