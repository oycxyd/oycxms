function [mzs_output, VLM, VLM_ind, ppm_dis_ini] = msreffind(mode, options)
%MSREFFIND Find shared virtual lock masses across a batch of spectra.
%
% Inputs
% ------
% mode:
%   'workspace' : use options.wsdatasets and options.wsmzs
%   'raw'       : search options.path for *.raw/datacube.h5
%   'recal'     : search current folder for *.raw and *_recal.h5
%   'mz5'       : search options.path for *.mz5, using matching .h5 files
%
% Options
% -------
% threshold : maximum matching difference in ppm
% minfreq   : minimum fraction of files in which a reference must match
%             Default: 0.5
% wsdatasets, wsmzs : required for workspace mode
%
% Outputs
% -------
% mzs_output : recalibrated m/z vectors, in the original input order
% VLM        : [per-file matched m/z columns, common-reference column]
% VLM_ind    : indices of each VLM in each original m/z vector
% ppm_dis_ini: mean fitted coefficient for delta(m/z) ~ k*sqrt(m/z)
%
% If no common lock masses survive the selected threshold:
%   VLM = zeros(0, nFiles+1)
%   VLM_ind = zeros(0, nFiles)
%   ppm_dis_ini = NaN
%   mzs_output = input mzs, unchanged

arguments
    mode (1,:) char
    options.path (1,:) char = pwd
    options.threshold (1,1) double {mustBePositive, mustBeFinite} = 100
    options.minfreq (1,1) double {mustBeGreaterThanOrEqual(options.minfreq,0), ...
                                   mustBeLessThanOrEqual(options.minfreq,1)} = 0.9
    options.wsmzs cell = {}
end

% -------------------------------------------------------------------------
% Load data according to mode
% -------------------------------------------------------------------------
datasets = {};
mzs = {};

switch lower(mode)
    case 'workspace'
        mzs = options.wsmzs;

    case 'raw'
        base_path = options.path;
        files = dir(fullfile(base_path, '*.raw'));

        for f = 1:numel(files)
            file_path = fullfile(files(f).folder, files(f).name);
            [datasets{f}, ~, mzs{f}] = h5toMat( ...
                fullfile(file_path, 'datacube.h5'));
        end

    case 'recal'
        base_path = options.path;
        raw_files = dir(fullfile(base_path, '*.raw'));
        recal_files = dir(fullfile(base_path, '*_recal.h5'));
        files = [raw_files; recal_files];

        for f = 1:numel(files)
            file_path = fullfile(files(f).folder, files(f).name);
            [~, ~, ext] = fileparts(file_path);

            if strcmpi(ext, '.raw')
                [datasets{f}, ~, mzs{f}] = h5toMat( ...
                    fullfile(file_path, 'datacube_recal.h5'));
            else
                [datasets{f}, ~, mzs{f}] = h5toMat(file_path);
            end
        end

    case {'mz5', '.mz5'}
        base_path = options.path;
        files = dir(fullfile(base_path, '*.mz5'));

        for f = 1:numel(files)
            [~, file_stem] = fileparts(files(f).name);
            h5_file = fullfile(base_path, [file_stem, '.h5']);
            [datasets{f}, ~, mzs{f}] = h5toMat(h5_file);
        end

    otherwise
        error('msreffind:UnknownMode', ...
            'Unsupported mode "%s". Use workspace, raw, recal, or mz5.', mode);
end

% -------------------------------------------------------------------------
% Validate inputs
% -------------------------------------------------------------------------
nFiles = numel(mzs);

if nFiles == 0
    error('msreffind:NoFiles', ...
        'No spectra/mz vectors were provided or found for mode "%s".', mode);
end

for f = 1:nFiles
    if isempty(mzs{f})
        error('msreffind:EmptyMZ', ...
            'mzs{%d} is empty. Every file must contain a nonempty m/z vector.', f);
    end

    mzs{f} = double(mzs{f}(:));

    if any(~isfinite(mzs{f}))
        error('msreffind:InvalidMZ', ...
            'mzs{%d} contains NaN or Inf values.', f);
    end

    if any(diff(mzs{f}) < 0)
        [mzs{f}, order] = sort(mzs{f});

        if ~isempty(datasets) && ~isempty(datasets{f})
            datasets{f} = datasets{f}(:, order);
        end
    end
end

% -------------------------------------------------------------------------
% Sort internally: use longest m/z axis as candidate reference source.
% Keep I so outputs are restored to original order at the end.
% -------------------------------------------------------------------------
mz_lengths = cellfun(@numel, mzs);
[~, I] = sort(mz_lengths, 'descend');

mzs_sorted = mzs(I);
reference_mz = mzs_sorted{1};
nCandidates = numel(reference_mz);

% Row m = one candidate VLM.
% Column f = its match in file f.
% NaN means "not matched in that file".
matches = NaN(nCandidates, nFiles);
matches(:, 1) = reference_mz;

% -------------------------------------------------------------------------
% Match candidate reference m/z values to each file
% -------------------------------------------------------------------------
for f = 2:nFiles
    this_mz = mzs_sorted{f};

    for m = 1:nCandidates
        reference_value = reference_mz(m);

        [abs_diff, nearest_ind] = min(abs(this_mz - reference_value));
        ppm_diff = abs_diff / reference_value * 1e6;

        if ppm_diff <= options.threshold
            matches(m, f) = this_mz(nearest_ind);
        end
    end
end

% -------------------------------------------------------------------------
% Keep VLM candidates found in at least minfreq of files.
% This handles low ppm thresholds safely.
% -------------------------------------------------------------------------
match_frequency = mean(isfinite(matches), 2);
keep = match_frequency >= options.minfreq;

matches = matches(keep, :);

if isempty(matches)
    % Preserve expected dimensions even when no VLMs exist.
    mzs_output = mzs;
    VLM = zeros(0, nFiles + 1);
    VLM_ind = zeros(0, nFiles);
    ppm_dis_ini = NaN;

    warning('msreffind:NoVLM', ...
        ['No virtual lock masses found at %.6g ppm across at least %.0f%% ' ...
         'of %d files. Returning unchanged m/z vectors.'], ...
        options.threshold, 100 * options.minfreq, nFiles);
    return
end

% -------------------------------------------------------------------------
% Remove duplicate candidates.
% Use the second file if available, otherwise the reference file.
% -------------------------------------------------------------------------
if nFiles >= 2
    duplicate_key = matches(:, 2);
else
    duplicate_key = matches(:, 1);
end

% Candidate rows not seen in the selected duplicate key should not survive
% uniqueness selection as NaN entries.
valid_key = isfinite(duplicate_key);

matches = matches(valid_key, :);

if isempty(matches)
    mzs_output = mzs;
    VLM = zeros(0, nFiles + 1);
    VLM_ind = zeros(0, nFiles);
    ppm_dis_ini = NaN;

    warning('msreffind:NoUniqueVLM', ...
        'No unique virtual lock masses remain at %.6g ppm.', options.threshold);
    return
end

[~, ia] = unique(duplicate_key(valid_key), 'stable');
matches = matches(ia, :);

% Calculate a common target lock mass for every retained row.
common_reference = mean(matches, 2, 'omitnan');

% -------------------------------------------------------------------------
% Recalibrate each m/z vector.
% VLM_ind has columns in sorted internal order initially.
% -------------------------------------------------------------------------
nVLM = size(matches, 1);
mzs_recal_sorted = mzs_sorted;
VLM_ind_sorted = zeros(nVLM, nFiles);
ini_sorted = NaN(nFiles, 1);

for f = 1:nFiles
    original_mz = mzs_sorted{f};
    recal_mz = original_mz;

    valid_match = isfinite(matches(:, f));

    if ~any(valid_match)
        % No references in this file: leave it unchanged.
        mzs_recal_sorted{f} = recal_mz;
        continue
    end

    matched_values = matches(valid_match, f);
    [is_present, location] = ismember(matched_values, original_mz);

    valid_locations = location(is_present);
    valid_targets = common_reference(valid_match);
    valid_targets = valid_targets(is_present);

    VLM_ind_sorted(valid_match, f) = location;

    if isempty(valid_locations)
        mzs_recal_sorted{f} = recal_mz;
        continue
    end

    % Substitute matched values with the shared reference values.
    recal_mz(valid_locations) = valid_targets;

    % Estimate delta(m/z) = k*sqrt(m/z) + b from the matched VLMs.
    % Need at least two points to estimate both k and b.
    displacement = recal_mz(valid_locations) - original_mz(valid_locations);
    usable = isfinite(displacement) & isfinite(original_mz(valid_locations));

    if nnz(usable) >= 2
        xs = sqrt(original_mz(valid_locations(usable)));
        A = [xs(:), ones(nnz(usable), 1)];
        p = A \ displacement(usable);

        ini_sorted(f) = p(1);
    elseif nnz(usable) == 1
        % With one point only, fit a zero-intercept k*sqrt(m/z) model.
        x0 = sqrt(original_mz(valid_locations(usable)));
        d0 = displacement(usable);

        if x0 ~= 0
            ini_sorted(f) = d0 / x0;
        end
    end

    mzs_recal_sorted{f} = recal_mz';
end

% -------------------------------------------------------------------------
% Restore original file order.
% I maps sorted positions -> original positions.
% -------------------------------------------------------------------------
mzs_output = cell(1, nFiles);
mzs_output(I) = mzs_recal_sorted;

VLM_ind = zeros(nVLM, nFiles);
VLM_ind(:, I) = VLM_ind_sorted;

VLM = [matches, common_reference];
ppm_dis_ini = mean(ini_sorted, 'omitnan');

end