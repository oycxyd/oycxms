function [scores, labels, sampledFiles, sampledRows] = qc_pca_aligned_batch(options)
%QC_PCA_ALIGNED_BATCH Sample spectra from aligned HDF5 files for PCA QC.
%
% Loads one HDF5 file at a time through h5toMat, randomly samples spectra,
% stores only sampled spectra, and clears the full dataset before reading
% the next file.
%
% Replace the PCA placeholder near the end with your PCA visualisation
% function, or pass one via options.pcaFunction.
%
% Example:
% qc_pca_aligned_batch( ...
%     path="C:\Data\Aligned", ...
%     pattern="*_aligned.h5", ...
%     spectraPerFile=100, ...
%     pcaFunction=@your_pca_function);

arguments
    options.path (1,1) string = string(pwd)
    options.pattern (1,1) string = "*_aligned.h5"
    options.spectraPerFile (1,1) double {mustBeInteger,mustBePositive} = 100
    options.maxTotalSpectra (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.randomSeed = 42
    options.normalise (1,1) logical = true
    options.scale (1,1) logical = true
    options.verbose (1,1) logical = true
end

if ~isempty(options.randomSeed)
    rng(options.randomSeed, "twister");
end

fileInfo = dir(fullfile(options.path, options.pattern));

if isempty(fileInfo)
    error("qc_pca_aligned_batch:NoFiles", ...
        "No files matching %s were found in %s.", ...
        options.pattern, options.path);
end

sampledFiles = arrayfun(@(f) fullfile(f.folder, f.name), ...
    fileInfo, "UniformOutput", false);

nFiles = numel(sampledFiles);
sampledRows = cell(1, nFiles);

% Store each file's sampled spectra temporarily as cells, then concatenate.
sampledSpectra = cell(1, nFiles);
sampledLabels = cell(1, nFiles);

nFeatures = [];

for f = 1:nFiles
    filePath = sampledFiles{f};

    if options.verbose
        fprintf("Loading file %d of %d: %s\n", ...
            f, nFiles, filePath);
    end

    % Your normal loader. Assumes:
    % dataset = nSpectra x nMzFeatures
    % mz      = 1 x nMzFeatures or nMzFeatures x 1
    [dataset, ~, mz] = h5toMat(filePath);

    nSpectra = size(dataset, 1);
    thisNFeatures = size(dataset, 2);

    if isempty(nFeatures)
        nFeatures = thisNFeatures;
    elseif thisNFeatures ~= nFeatures
        clear dataset mz
        error("qc_pca_aligned_batch:FeatureCountMismatch", ...
            ["File %s has %d features, but prior aligned files have %d. " ...
             "PCA requires a common aligned m/z feature axis."], ...
            filePath, thisNFeatures, nFeatures);
    end

    if numel(mz) ~= thisNFeatures
        clear dataset mz
        error("qc_pca_aligned_batch:MZLengthMismatch", ...
            ["File %s has %d m/z entries but %d intensity columns. " ...
             "Check the h5toMat output orientation."], ...
            filePath, numel(mz), thisNFeatures);
    end

    nTake = min(options.spectraPerFile, nSpectra);
    selected = sort(randperm(nSpectra, nTake));

    sampledRows{f} = selected;
    sampledSpectra{f} = single(dataset(selected, :));

    [~, fileLabel] = fileparts(filePath);
    sampledLabels{f} = repmat(string(fileLabel), nTake, 1);

    if options.verbose
        fprintf("  Sampled %d of %d spectra.\n", nTake, nSpectra);
    end

    % Important: release the full-file data before loading the next file.
    clear dataset mz selected
end

X = vertcat(sampledSpectra{:});
labels = vertcat(sampledLabels{:});

% The temporary cells can be discarded once X has been assembled.
clear sampledSpectra sampledLabels

% Optional overall cap after balanced per-file sampling.
if options.maxTotalSpectra > 0 && size(X,1) > options.maxTotalSpectra
    keep = randperm(size(X,1), options.maxTotalSpectra);
    X = X(keep, :);
    labels = labels(keep);
end

% Remove invalid / empty sampled spectra.
validRows = all(isfinite(X), 2) & any(X ~= 0, 2);

if options.verbose && any(~validRows)
    fprintf("Removing %d invalid or all-zero sampled spectra.\n", ...
        nnz(~validRows));
end

X = X(validRows, :);
labels = labels(validRows);

if size(X,1) < 3
    error("qc_pca_aligned_batch:TooFewSpectra", ...
        "Fewer than three valid spectra remain after filtering.");
end

% Optional feature-wise standardization for PCA.
if options.normalise
    X = TIC_norm(X);
    if options.scale
        X = log_trans(X);
    end
else
    X = double(X);
end

% -------------------------------------------------------------------------
% PCA visualisation
% 
% -------------------------------------------------------------------------
[pc,scores,v]=pca(X);
figure( ...
    "Name", "PCA QC: aligned spectra by file", ...
    "Color", "w", ...
    "Position", [100, 100, 1350, 750]);

hold on

fileLabels = unique(labels, "stable");
nFiles = numel(fileLabels);

% Deterministic, broadly separated hues.
hue = mod((0:nFiles-1)' * 0.61803398875, 1);
colours = hsv2rgb([ ...
    hue, ...
    0.70 * ones(nFiles,1), ...
    0.85 * ones(nFiles,1)]);

markers = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};
markerForFile = markers(mod(0:nFiles-1, numel(markers)) + 1);

for f = 1:nFiles
    mask = labels == fileLabels(f);

    scatter( ...
        scores(mask, 1), ...
        scores(mask, 2), ...
        24, ...
        colours(f, :), ...
        markerForFile{f}, ...
        "filled", ...
        "MarkerFaceAlpha", 0.55, ...
        "MarkerEdgeAlpha", 0.40, ...
        "DisplayName", fileLabels(f));
end

xlabel("PC 1 score");
ylabel("PC 2 score");
title("PCA QC of randomly sampled aligned spectra");
grid on
box on

legend( ...
    "Location", "eastoutside", ...
    "Interpreter", "none", ...
    "FontSize", 8);

hold off

end