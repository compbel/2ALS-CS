function k = v2Exact(x)
%V2EXACT Exact 2-adic valuation of a nonzero symbolic integer or rational.
%   For x = u/v, returns v2(u)-v2(v). v2Exact(0) is intentionally undefined.

x = sym(x);
[u, v] = numden(x);
if double(u) == 0
    error('padic:v2Zero', 'v2(0) is infinity and is not accepted here.');
end
k = trailingTwos(u) - trailingTwos(v);
end

function k = trailingTwos(x)
x = abs(x);
numeric = double(x);
if isfinite(numeric) && numeric == fix(numeric) && numeric < 2^50
    k = 0;
    while mod(numeric, 2) == 0
        numeric = numeric / 2;
        k = k + 1;
    end
    return;
end
k = 0;
while isAlways(mod(x, sym(2)) == 0)
    x = x / 2;
    k = k + 1;
end
end
