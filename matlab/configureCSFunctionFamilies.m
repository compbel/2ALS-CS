function [primaryFamilies, escapeFamilies, singleFamilies, pairFamilies, mode] = ...
        configureCSFunctionFamilies(instance, familyBuilder, useColumnPairFamilies)
%CONFIGURECSFUNCTIONFAMILIES Select the Camin-Sokal column-family schedule.
%   [PRIMARY, ESCAPE, SINGLES, PAIRS, MODE] =
%   CONFIGURECSFUNCTIONFAMILIES(INSTANCE, BUILDER, USECOLUMNPAIRFAMILIES)
%   supports three case-insensitive modes:
%       'NEVER'  - PRIMARY contains only single-column families and PAIRS
%                  is not constructed.
%       'LOCAL'  - PRIMARY contains single-column families and ESCAPE
%                  contains all allelic-state column-pair families.
%       'ALWAYS' - PRIMARY contains both single-column and
%                  column-pair families; ESCAPE is empty.
%
%   BUILDER is called once for the instance. In NEVER mode it is requested
%   to return only its first output, avoiding construction of pair families.

if ~isa(familyBuilder, 'function_handle')
    error('padic:familyBuilder', 'familyBuilder must be a function handle.');
end
mode = normalizeColumnPairFamilyMode(useColumnPairFamilies);

escapeFamilies = struct.empty(0, 1);
pairFamilies = struct.empty(0, 1);
switch mode
    case 'NEVER'
        singleFamilies = familyBuilder(instance);
        primaryFamilies = singleFamilies;
    case 'LOCAL'
        [singleFamilies, pairFamilies] = familyBuilder(instance);
        primaryFamilies = singleFamilies;
        escapeFamilies = pairFamilies;
    case 'ALWAYS'
        [singleFamilies, pairFamilies] = familyBuilder(instance);
        primaryFamilies = concatenateFamilyStructs(singleFamilies, pairFamilies);
end
end

function combined = concatenateFamilyStructs(first, second)
% Preserve family-specific replay metadata while making fields compatible.
if isempty(first)
    combined = second(:);
    return;
elseif isempty(second)
    combined = first(:);
    return;
end
if ~isstruct(first) || ~isstruct(second)
    error('padic:familyFormat', ...
        'The family builder outputs must be structure arrays.');
end
fieldOrder = [fieldnames(first); setdiff(fieldnames(second), ...
    fieldnames(first), 'stable')];
first = addMissingFields(first(:), fieldOrder);
second = addMissingFields(second(:), fieldOrder);
combined = [orderfields(first, fieldOrder); orderfields(second, fieldOrder)];
end

function families = addMissingFields(families, fieldOrder)
for index = 1:numel(fieldOrder)
    name = fieldOrder{index};
    if ~isfield(families, name)
        [families.(name)] = deal([]);
    end
end
end
