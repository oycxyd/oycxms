function [raw_specs, dims] = raw2mat(filename)
% filename = '2019_09_01_porklivertest_100umpixel_neg Analyte 1.raw';

    [~,p2] = watersPackages(filename);
    numS = calllib('MassLynxRaw','getScansInFunction',p2,1);

    raw_specs = {};xy = {};
    for n = 1:numS
        [mz, spectrum,xy{n}] = readraw2spec(filename,n);
        raw_specs{n} = cat(1,mz,spectrum);
    end

    xy = cell2mat(xy);
    pixel = abs(xy(1)-xy(3));
    x_len = round(abs(xy(end-1)-xy(1))/pixel)+1;
    y_len = round(abs(xy(end)-xy(2))/pixel)+1;

    dims = [x_len,y_len];
    
    if length(specs) ~= x_len.*y_len
        disp('Warning! The pixel dimension according to the metadata does not match the file size.')
    end
end