function mode = normalizeColumnPairFamilyMode(value)
%NORMALIZECOLUMNPAIRFAMILYMODE Validate NEVER/LOCAL/ALWAYS scheduling mode.
if ischar(value) && isrow(value)
    mode = upper(value);
elseif isstring(value) && isscalar(value)
    mode = upper(char(value));
else
    error('padic:columnPairMode', ...
        'useColumnPairFamilies must be NEVER, LOCAL, or ALWAYS.');
end
if ~any(strcmp(mode, {'NEVER', 'LOCAL', 'ALWAYS'}))
    error('padic:columnPairMode', ...
        'useColumnPairFamilies must be NEVER, LOCAL, or ALWAYS.');
end
end
