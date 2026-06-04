function [raw_specs, dims, chromoT, chromoI] = mz5toMat_row(fileList)

if isstring(fileList)
    fileList = cellstr(fileList);
end

if isstruct(fileList)
    fileList = {fileList.name};
end

nFiles = length(fileList);

raw_specs = {};
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
            spec = h5read(filename,'/SpectrumIntensity',1,double(Sindex(n)));
            mz   = h5read(filename,'/SpectrumMZ',1,double(Sindex(n)));
        else
            start = double(Sindex(n-1)) + 1;
            count = double(Sindex(n)) - double(Sindex(n-1));

            spec = h5read(filename,'/SpectrumIntensity',start,count);
            mz   = h5read(filename,'/SpectrumMZ',start,count);
        end

        for i = 2:length(mz)
            mz(i) = mz(i-1) + mz(i);
        end

        raw_specs_line{n} = [mz'; spec'];
    end

    raw_specs = [raw_specs, raw_specs_line];

end

if numel(unique(nSpec_per_file)) ~= 1
    error('failed to align pixel and generate dims');
end

%dims = [nFiles, nSpec_per_file(1)];
dims = [nSpec_per_file(1),nFiles];
end