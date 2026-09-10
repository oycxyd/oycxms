function [ppmthresh, metric, featurescale, MIscore, details] = migsearch_batch(options)
%MIGSEARCH_BATCH Select a ppm threshold for batch m/z recalibration.
%
% Reads one input file at a time, retains only its m/z vector and mean
% spectrum, then evaluates candidate VLM thresholds. Full data matrices are
% cleared immediately after each mean spectrum is calculated.
%
% This version is intended for centroid or unequal-length m/z data. The
% feature count used in FEATURESCALE is the number of common m/z features
% returned by CWT2CMZ at each tested ppm threshold, matching the original
% MIGSEARCH logic.
%
% Required supporting functions on the MATLAB path:
%   h5toMat, msreffind, cwt2cmz, adaptmatch, micontcont
%
% MSREFFIND must support:
%   msreffind('workspace', wsmzs=mzs, threshold=ppm, minfreq=fraction)
%
% Example:
% [ppmthresh, metric, featurescale, MIscore, details] = migsearch_batch( ...
%     path="C:\\Data\\Batch", ...
%     mode="raw", ...
%     thresholds=[5 10 20 30 40 50 60 70 80 100 150 300], ...
%     mzlow=50, mzhigh=1200, freq=0.4, minVLMfreq=0.5);

arguments
    options.path (1,1) string = string(pwd)
    options.mode (1,1) string = "raw"
    options.thresholds (1,:) double {mustBePositive, mustBeFinite} = ...
        [5 10 20 30 40 50 60 70 80 100 150 300]
    options.mzlow (1,1) double {mustBeFinite} = 50
    options.mzhigh (1,1) double {mustBeFinite} = 1200
    options.freq (1,1) double {mustBeGreaterThanOrEqual(options.freq,0), ...
                                mustBeLessThanOrEqual(options.freq,1)} = 0.4
    options.minVLMfreq (1,1) double {mustBeGreaterThanOrEqual(options.minVLMfreq,0), ...
                                      mustBeLessThanOrEqual(options.minVLMfreq,1)} = 0.5
    options.miBins (1,1) double {mustBeInteger, mustBePositive} = 5
    options.verbose (1,1) logical = true
end

if options.mzlow >= options.mzhigh
    error('migsearch_batch:InvalidMZRange', ...
        'mzlow must be less than mzhigh.');
end

[filenames, h5files] = localFindFiles(options.path, options.mode);
nFiles = numel(h5files);

if nFiles < 2
    error('migsearch_batch:TooFewFiles', ...
        'At least two files are needed for a batch threshold search.');
end

% -------------------------------------------------------------------------
% Load only m/z vectors and mean spectra. Do not retain batch datacubes.
% -------------------------------------------------------------------------
mzs = cell(1, nFiles);
meanSpecs = cell(1, nFiles);

for f = 1:nFiles
    if options.verbose
        fprintf('Loading mean spectrum: file %d of %d\n', f, nFiles);
        fprintf('  %s\n', filenames{f});
    end

    try
        [dataset, ~, mz] = h5toMat(h5files{f});
    catch ME
        error('migsearch_batch:ReadFailed', ...
            'Could not read %s: %s', h5files{f}, ME.message);
    end

    mz = double(mz(:));

    if isempty(mz) || any(~isfinite(mz))
        clear dataset mz
        error('migsearch_batch:InvalidMZ', ...
            'File %d has an empty or non-finite m/z vector.', f);
    end

    if isvector(dataset)
        meanSpec = double(dataset(:))';
    else
        meanSpec = mean(double(dataset), 1);
    end

    if numel(mz) ~= numel(meanSpec)
        datasetSize = size(dataset);
        error('migsearch_batch:DimensionMismatch', ...
            ['File %d has %d m/z values but %d mean-intensity values. ' ...
             'Datacube size: [%s].'], ...
            f, numel(mz), numel(meanSpec), num2str(datasetSize));
    end

    % CWT2CMZ and interpolation/matching are safest with ascending m/z.
    if any(diff(mz) < 0)
        [mz, order] = sort(mz);
        meanSpec = meanSpec(order);
    end

    mzs{f} = mz;
    meanSpecs{f} = meanSpec(:)';

    % Explicitly release the potentially large single-file datacube.
    clear dataset mz meanSpec order
end

mzLengths = cellfun(@numel, mzs);
meanLength = median(mzLengths);
[~, ind1] = min(abs(mzLengths - meanLength));

thresholds = options.thresholds(:)';
nThresholds = numel(thresholds);
minFreqCount = max(1, ceil(nFiles * options.freq));

% One MI comparison for every non-reference file and every threshold.
metric = NaN(nFiles - 1, nThresholds);

% One common-feature count per threshold. This is deliberately independent
% of interpolation-grid size and is obtained from cwt2cmz.
nfeatures = NaN(1, nThresholds);
numVLM = zeros(1, nThresholds);
ppmDisIni = NaN(1, nThresholds);
validThreshold = false(1, nThresholds);

for t = 1:nThresholds
    threshold = thresholds(t);

    if options.verbose
        fprintf('Testing VLM threshold %.6g ppm (%d of %d)\n', ...
            threshold, t, nThresholds);
    end

    % MSREFFIND only receives the small m/z-vector cell array.
    try
        [mzsRecal, VLM, VLMind, ppmDisIni(t)] = msreffind( ...
            'workspace', ...
            wsmzs=mzs, ...
            threshold=threshold, ...
            minfreq=options.minVLMfreq);
    catch ME
        warning('migsearch_batch:VLMFailure', ...
            'Skipping %.6g ppm because msreffind failed: %s', ...
            threshold, ME.message);
        continue
    end

    % A strict threshold can legitimately give no VLMs.
    if isempty(VLM) || isempty(VLMind) || ~iscell(mzsRecal) || ...
            numel(mzsRecal) ~= nFiles
        if options.verbose
            fprintf('  No usable virtual lock masses; skipped.\n');
        end
        continue
    end

    % Check that recalibration retained the one-to-one intensity/mz layout.
    recalValid = true;
    for f = 1:nFiles
        if numel(mzsRecal{f}) ~= numel(meanSpecs{f})
            warning('migsearch_batch:LengthChanged', ...
                ['Skipping %.6g ppm: recalibrated m/z length differs from ' ...
                 'mean-spectrum length in file %d (%d versus %d).'], ...
                threshold, f, numel(mzsRecal{f}), numel(meanSpecs{f}));
            recalValid = false;
            break
        end
    end

    if ~recalValid
        continue
    end

    % CWT2CMZ supplies the common centroid/features at this ppm threshold.
    try
        mzRecalCommon = cwt2cmz( ...
            mzsRecal, ...
            minFreq=minFreqCount, ...
            mzRange=[options.mzlow, options.mzhigh], ...
            ppm=threshold);
    catch ME
        warning('migsearch_batch:CWT2CMZFailure', ...
            'Skipping %.6g ppm because cwt2cmz failed: %s', ...
            threshold, ME.message);
        continue
    end

    mzRecalCommon = double(mzRecalCommon(:));
    mzRecalCommon = mzRecalCommon(isfinite(mzRecalCommon));
    mzRecalCommon = unique(mzRecalCommon, 'sorted');

    nfeatures(t) = numel(mzRecalCommon);

    if nfeatures(t) < options.miBins
        if options.verbose
            fprintf('  Only %d common features; skipped.\n', nfeatures(t));
        end
        continue
    end

    % ---------------------------------------------------------------------
    % Convert mean spectra to a shared feature-vector representation.
    % This is the same principle as the original MIGSEARCH: adapt common
    % CWT m/z features, then overwrite known VLM locations exactly.
    % ---------------------------------------------------------------------
    featureSpecs = ones(nFiles, nfeatures(t));

    % VLM(:,end) is the common/mean VLM location in the revised msreffind.
    % If the older msreffind layout is used, this still uses the final
    % common-reference column as in the original code.
    vlmCommon = double(VLM(:, end));
    [~, vlmFeatureInd] = ismember(vlmCommon, mzRecalCommon);

    for f = 1:nFiles
        currentMz = double(mzsRecal{f}(:));
        currentSpec = double(meanSpecs{f}(:));

        % Adapt common features to the current file's recalibrated axis.
        try
            [indd1, indd2] = adaptmatch(currentMz, mzRecalCommon', ppmDisIni(t));
        catch ME
            warning('migsearch_batch:AdaptMatchFailure', ...
                ['Threshold %.6g ppm, file %d: adaptmatch failed: %s'], ...
                threshold, f, ME.message);
            continue
        end

        % Guard against malformed indices returned by a helper.
        validMatch = indd1 >= 1 & indd1 <= numel(currentSpec) & ...
                     indd2 >= 1 & indd2 <= nfeatures(t);
        indd1 = indd1(validMatch);
        indd2 = indd2(validMatch);

        featureSpecs(f, indd2) = currentSpec(indd1);

        % Preserve intensities at the explicitly identified VLM indices.
        if size(VLMind, 2) >= f
            vlmSourceInd = VLMind(:, f);
            validVLM = vlmSourceInd >= 1 & ...
                       vlmSourceInd <= numel(currentSpec) & ...
                       vlmFeatureInd >= 1 & ...
                       vlmFeatureInd <= nfeatures(t);

            featureSpecs(f, vlmFeatureInd(validVLM)) = ...
                currentSpec(vlmSourceInd(validVLM));
        end
    end

    refSpec = featureSpecs(ind1, :);
    comparisonRow = 0;

    for f = 1:nFiles
        if f == ind1
            continue
        end

        comparisonRow = comparisonRow + 1;
        otherSpec = featureSpecs(f, :);

        % The original function uses MI(ref, other) normalized by MI(ref,ref).
        try
            selfMI = mi_cont_cont(refSpec, refSpec, options.miBins);
            otherMI = mi_cont_cont(refSpec, otherSpec, options.miBins);

            if isfinite(selfMI) && selfMI > 0 && isfinite(otherMI)
                metric(comparisonRow, t) = otherMI / selfMI;
            end
        catch ME
            warning('migsearch_batch:MIComputationFailed', ...
                'Threshold %.6g ppm, file %d: MI failed: %s', ...
                threshold, f, ME.message);
        end
    end

    validThreshold(t) = any(isfinite(metric(:, t)));

    if options.verbose
        fprintf('  Shared VLMs: %d; common CWT features: %d; mean MI: %.6f\n', ...
            size(VLM, 1), nfeatures(t), mean(metric(:, t), 'omitnan'));
    end
end

% -------------------------------------------------------------------------
% Reproduce the intended original feature scaling:
% compare the CWT common-feature count at each threshold with the median
% m/z-vector length. This is appropriate for centroid data, where vector
% length is itself a meaningful feature-count scale.
% -------------------------------------------------------------------------
featurescale = abs(nfeatures - meanLength);
meanMetric = mean(metric, 1, 'omitnan');
% 
% finiteScale = isfinite(featurescale) & validThreshold;
% if any(finiteScale)
%     maxScale = max(featurescale(finiteScale));
%     if maxScale > 0
%         scalePenalty = NaN(1, nThresholds);
%         scalePenalty(finiteScale) = featurescale(finiteScale) ./ maxScale;
%     else
%         scalePenalty = NaN(1, nThresholds);
%         scalePenalty(finiteScale) = 0;
%     end
% else
%     scalePenalty = NaN(1, nThresholds);
% end

% MIscore = meanMetric .* (1 - scalePenalty);

MIscore = meanMetric;
MIscore(~validThreshold) = NaN;

if ~any(isfinite(MIscore))
    error('migsearch_batch:NoValidThreshold', ...
        ['No candidate threshold produced usable shared VLMs, common CWT ' ...
         'features, and MI scores. Increase thresholds or reduce minVLMfreq/freq.']);
end

validInd = find(isfinite(MIscore));
[~, localBest] = max(MIscore(validInd));
bestInd = validInd(localBest);
ppmthresh = thresholds(bestInd);

if options.verbose
    fprintf('Selected ppm threshold: %.6g\n', ppmthresh);
    fprintf('  Mean MI: %.6f\n', meanMetric(bestInd));
    fprintf('  Common CWT features: %d\n', nfeatures(bestInd));
end

% No full datacubes are returned or retained.
details = struct();
details.filenames = filenames;
details.h5files = h5files;
details.mzs = mzs;
details.mean_specs = meanSpecs;
details.mz_lengths = mzLengths;
details.mean_length = meanLength;
details.reference_file = ind1;
details.reference_filename = filenames{ind1};
details.thresholds = thresholds;
details.nfeatures = nfeatures;
details.num_vlm = numVLM;
details.ppm_dis_ini = ppmDisIni;
details.valid_threshold = validThreshold;
details.mean_metric = meanMetric;
% details.scale_penalty = scalePenalty;
details.selected_index = bestInd;
end

function [ind1,ind2] = adaptmatch(mz, mz_recal,ini)
    for n = 1:length(mz_recal)
        [dif(n),ind1(n)] = min(abs(mz-mz_recal(n)));
    end
    ppm_dis = ini*sqrt(mz_recal);
    margin = 1.5;
    matched = (dif<= ppm_dis*margin);
    ind1 = ind1(matched);[ind1,iu] = unique(ind1);
    ind2 = (1:n);ind2 = ind2(matched);ind2 = ind2(iu);
end

function [filenames, h5files] = localFindFiles(basePath, mode)
%LOCALFINDFILES Return display names and HDF5 files passed to h5toMat.

basePath = char(basePath);
mode = lower(string(mode));

switch mode
    case "raw"
        listing = dir(fullfile(basePath, '*.raw'));
        filenames = arrayfun(@(f) fullfile(f.folder, f.name), ...
            listing, 'UniformOutput', false);
        h5files = cellfun(@(p) fullfile(p, 'datacube.h5'), filenames, ...
            'UniformOutput', false);

    case "recal"
        rawListing = dir(fullfile(basePath, '*.raw'));
        recalListing = dir(fullfile(basePath, '*_recal.h5'));

        rawNames = arrayfun(@(f) fullfile(f.folder, f.name), ...
            rawListing, 'UniformOutput', false);
        recalNames = arrayfun(@(f) fullfile(f.folder, f.name), ...
            recalListing, 'UniformOutput', false);

        filenames = [rawNames, recalNames];
        h5files = cell(1, numel(filenames));

        for k = 1:numel(rawNames)
            h5files{k} = fullfile(rawNames{k}, 'datacube_recal.h5');
        end

        for k = 1:numel(recalNames)
            h5files{numel(rawNames) + k} = recalNames{k};
        end

    case {"mz5", ".mz5"}
        listing = dir(fullfile(basePath, '*.mz5'));
        filenames = arrayfun(@(f) fullfile(f.folder, f.name), ...
            listing, 'UniformOutput', false);
        h5files = cell(1, numel(filenames));

        for k = 1:numel(filenames)
            [folder, stem] = fileparts(filenames{k});
            h5files{k} = fullfile(folder, [stem, '.h5']);
        end

    case {"h5", ".h5"}
        listing = dir(fullfile(basePath, '*.h5'));
        filenames = arrayfun(@(f) fullfile(f.folder, f.name), ...
            listing, 'UniformOutput', false);
        h5files = filenames;

    otherwise
        error('migsearch_batch:UnknownMode', ...
            'Unknown mode "%s". Use raw, recal, mz5, or h5.', mode);
end

if isempty(h5files)
    error('migsearch_batch:NoFiles', ...
        'No files found in %s for mode "%s".', basePath, mode);
end

missing = h5files(~cellfun(@isfile, h5files));
if ~isempty(missing)
    error('migsearch_batch:MissingFile', ...
        'Expected HDF5 input does not exist:\n%s', strjoin(missing, '\n'));
end
end
