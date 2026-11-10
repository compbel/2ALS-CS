%RUN_EXAMPLE Infer a rooted tree for the bundled binary alignment.
exampleDirectory = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(exampleDirectory);
addpath(fullfile(repositoryRoot, 'matlab'));

alignment = readmatrix(fullfile(exampleDirectory, 'alignment.csv'));
rng(20261110, 'twister');
result = twoALSCS(alignment, ...
    'PairMode', 'ALWAYS', ...
    'ParallelFamilies', false);

fprintf('Initial score: %g\n', result.initialScore);
fprintf('Final score:   %g\n', result.finalScore);
fprintf('Accepted moves: %d\n', result.acceptedMoveCount);
fprintf('Termination: %s\n', result.terminationReason);
disp(result.finalTree.Edges);
