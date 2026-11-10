function catalog = loadPairAllelicStateTreeCatalog()
%LOADPAIRALLELICSTATETREECATALOG Load the cached allelic-state hierarchy table.
%   The validated 184-row catalog is bundled with the release.

persistent cachedCatalog
if isempty(cachedCatalog)
    functionDirectory = fileparts(mfilename('fullpath'));
    cacheFile = fullfile(functionDirectory, ...
        'pair_allelic_state_trees_kappa_all_orientations.mat');
    if ~isfile(cacheFile)
        error('padic:allelicStateCacheMissing', ...
            'The bundled allelic-state MAT catalog is missing.');
    end
    loaded = load(cacheFile, 'catalog');
    requiredFields = {'stateOrder', 'observedJointStates', 'stateMask', ...
        'stateCount', 'orientationIndex', 'orderedRootedTree', 'kappa', ...
        'd_H', 'rowCount'};
    if ~isfield(loaded, 'catalog') || ...
            ~all(isfield(loaded.catalog, requiredFields)) || ...
            loaded.catalog.rowCount ~= 184 || ...
            ~isequal(size(loaded.catalog.kappa), [184, 4])
        error('padic:allelicStateCacheInvalid', ...
            'The allelic-state MAT cache is invalid or incomplete.');
    end
    cachedCatalog = loaded.catalog;
end
catalog = cachedCatalog;
end
