function data_select = loadMetadataScans(filename, metadata)
%LOADMETADATASCANS Load scan ranges specified in a metadata cell array.
%
%   data_select = loadMetadataScans(filename, metadata)
%
% Inputs
%   filename     Full path or filename of the raw/mz5 data file.
%   metadata     Cell array containing:
%                column 2: filename identifier to match
%                column 3: first scan number
%                column 4: last scan number
%   use_metadata Logical or numeric flag. When false/0, no scans are read.
%
% Output
%   data_select  N-by-1 cell array. Each cell contains a 2-by-M array:
%                row 1: m/z values
%                row 2: intensity values
%
% Dependencies
%   readraw2spec
%   checkSindex

data_select = {};

[~, baseName, ext] = fileparts(filename);
ext = lower(ext);

% Match either the full filename/path or just the basename stored in metadata.
metadataNames = string(metadata(:, 2));
matchingRows = find( ...
    contains(metadataNames, string(filename)) | ...
    contains(metadataNames, baseName));

if isempty(matchingRows)
    warning("loadMetadataScans:NoMetadataMatch", ...
        "No metadata rows matched '%s'.", filename);
    return
end

% Load and validate the spectrum index once for mz5 files.
if strcmp(ext, ".mz5")
    spectrumIndex = double(h5read(filename, "/SpectrumIndex"));
    spectrumIndex = checkSindex(spectrumIndex);
elseif ~strcmp(ext, ".raw")
    error("loadMetadataScans:UnsupportedFileType", ...
        "Unsupported file extension '%s'. Expected .raw or .mz5.", ext);
end

for row = matchingRows(:)'

    firstScan = metadata{row, 3};
    lastScan  = metadata{row, 4};

    if iscell(firstScan)
        firstScan = firstScan{1};
    end

    if iscell(lastScan)
        lastScan = lastScan{1};
    end

    firstScan = double(firstScan);
    lastScan = double(lastScan);

    % validateattributes(firstScan, {"numeric"}, ...
    %     {"scalar", "finite", "integer", "positive"}, ...
    %     mfilename, "metadata start scan");
    % 
    % validateattributes(lastScan, {"numeric"}, ...
    %     {"scalar", "finite", "integer", "positive"}, ...
    %     mfilename, "metadata end scan");

    if lastScan < firstScan
        warning("loadMetadataScans:InvalidRange", ...
            "Skipping metadata row %d because end scan (%d) is before start scan (%d).", ...
            row, lastScan, firstScan);
        continue
    end

    for scanNumber = firstScan:lastScan

        if strcmp(ext, ".raw")
            [mz, intensity] = readraw2spec(filename, scanNumber);

        else
            [mz, intensity] = readMz5Scan( ...
                filename, spectrumIndex, scanNumber);
        end

        % Ensure both are row vectors before concatenating.
        mz = mz(:)';
        intensity = intensity(:)';

        if numel(mz) ~= numel(intensity)
            warning("loadMetadataScans:LengthMismatch", ...
                ["Skipping scan %d from '%s': m/z length (%d) does not " ...
                 "match intensity length (%d)."], ...
                scanNumber, filename, numel(mz), numel(intensity));
            continue
        end

        data_select{end + 1, 1} = [mz; intensity];
    end
end

end


function [mz, intensity] = readMz5Scan(filename, spectrumIndex, scanNumber)
%READMZ5SCAN Read a single spectrum from an mz5 HDF5 file.

nScans = numel(spectrumIndex);

if scanNumber > nScans
    error("loadMetadataScans:ScanOutOfRange", ...
        "Requested scan %d, but SpectrumIndex contains only %d scans.", ...
        scanNumber, nScans);
end

if scanNumber == 1
    startIndex = 1;
    nPoints = spectrumIndex(1);
else
    startIndex = spectrumIndex(scanNumber - 1) + 1;
    nPoints = spectrumIndex(scanNumber) - spectrumIndex(scanNumber - 1);
end

if nPoints < 1
    error("loadMetadataScans:EmptyScan", ...
        "Scan %d contains no spectrum.", scanNumber);
end

intensity = h5read( ...
    filename, ...
    "/SpectrumIntensity", ...
    double(startIndex), ...
    double(nPoints));

mzDelta = h5read( ...
    filename, ...
    "/SpectrumMZ", ...
    double(startIndex), ...
    double(nPoints));

% mz5 stores the first m/z value followed by successive m/z increments.
mz = cumsum(double(mzDelta));

end