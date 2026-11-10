function [values, usable] = exactIntegerDoubles(input)
%EXACTINTEGERDOUBLES Convert small exact symbolic integers without rounding.
%   The conservative 2^50 bound leaves room for exact subtraction and other
%   integer operations in the numeric fast path. Larger values stay symbolic.

input = sym(input);
[~, denominator] = numden(input);
values = double(input);
usable = all(double(denominator(:)) == 1) && ...
    all(isfinite(values(:))) && all(values(:) == fix(values(:))) && ...
    all(abs(values(:)) < 2^50);
end
