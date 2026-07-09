function [raw_specs, dims, chromoT, chromoI, missing_mask] = mz5toMat_row(fileList)

if isstring(fileList)
    fileList = cellstr(fileList);
end

if isstruct(fileList)
    fileList = {fileList.name};
end

nFiles = length(fileList);

raw_specs_rows = cell(nFiles, 1);
chromoT = [];
chromoI = [];
nSpec_per_file = zeros(nFiles,1);

for f = 1:nFiles

    filename = fileList{f};

    chromoT_line = h5read(filename,'/ChomatogramTime');
    chromoI_line = h5read(filename,'/ChromatogramIntensity');
    Sindex = h5read(filename,'/SpectrumIndex');

    chromoT = [chromoT; chromoT_line(:)];
    chromoI = [chromoI; chromoI_line(:)];

    nSpec = length(Sindex);
    nSpec_per_file(f) = nSpec;

    raw_specs_line = cell(1, nSpec);

    parfor n = 1:nSpec

        if n == 1
            start = 1;
            count = double(Sindex(n));
        else
            start = double(Sindex(n-1)) + 1;
            count = double(Sindex(n)) - double(Sindex(n-1));
        end

        spec = h5read(filename,'/SpectrumIntensity',start,count);
        mz   = h5read(filename,'/SpectrumMZ',start,count);

        mz = cumsum(mz);

        raw_specs_line{n} = [mz(:)'; spec(:)'];
    end

    raw_specs_rows{f} = raw_specs_line;

end

disp('nSpec per file:')
disp(nSpec_per_file)
disp('Unique nSpec:')
disp(unique(nSpec_per_file))

nRows = nFiles;
nCols = max(nSpec_per_file);
dims = [nRows, nCols];

raw_specs_matrix = cell(nRows, nCols);
missing_mask = true(nRows, nCols);

for f = 1:nFiles
    nSpec = nSpec_per_file(f);
    raw_specs_matrix(f, 1:nSpec) = raw_specs_rows{f};
    missing_mask(f, 1:nSpec) = false;
end

for r = 1:nRows
    for c = 1:nCols
        if isempty(raw_specs_matrix{r,c})
            raw_specs_matrix{r,c} = fill_missing_spectrum_median(raw_specs_matrix, r, c);
        end
    end
end

raw_specs = reshape(raw_specs_matrix, 1, []);

end


function spec_filled = fill_missing_spectrum_median(raw_specs_matrix, r, c)

% First use 4-neighbour
offsets = [
    -1  0
     1  0
     0 -1
     0  1
];

neighbour_specs = collect_neighbours(raw_specs_matrix, r, c, offsets);

% If too few valid neighbours, use 8-neighbour
if numel(neighbour_specs) < 2
    offsets = [
        -1 -1
        -1  0
        -1  1
         0 -1
         0  1
         1 -1
         1  0
         1  1
    ];
    neighbour_specs = collect_neighbours(raw_specs_matrix, r, c, offsets);
end

% If still no neighbour, leave empty
if isempty(neighbour_specs)
    spec_filled = [];
    return
end

mz_ref = neighbour_specs{1}(1,:);
int_mat = zeros(numel(neighbour_specs), numel(mz_ref));

for k = 1:numel(neighbour_specs)

    mz_k = neighbour_specs{k}(1,:);
    int_k = neighbour_specs{k}(2,:);

    if numel(mz_k) == numel(mz_ref) && max(abs(mz_k - mz_ref)) < 1e-6
        int_mat(k,:) = int_k;
    else
        int_mat(k,:) = interp1(mz_k, int_k, mz_ref, 'linear', 0);
    end

end

int_median = median(int_mat, 1);

spec_filled = [mz_ref; int_median];

end


function neighbour_specs = collect_neighbours(raw_specs_matrix, r, c, offsets)

[nRows, nCols] = size(raw_specs_matrix);
neighbour_specs = {};

for i = 1:size(offsets,1)

    rr = r + offsets(i,1);
    cc = c + offsets(i,2);

    if rr < 1 || rr > nRows || cc < 1 || cc > nCols
        continue
    end

    if ~isempty(raw_specs_matrix{rr,cc})
        neighbour_specs{end+1} = raw_specs_matrix{rr,cc};
    end

end

end