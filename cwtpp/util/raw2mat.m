function [raw_specs, dims] = raw2mat(filename)
% load Waters raw data to Work
    [~,p2] = watersPackages(filename);
    numS = calllib('MassLynxRaw','getScansInFunction',p2,1);

    raw_specs = {};xy = {};
    counter1 = parfor_wait(numS, 'Waitbar', true);
    parfor n = 1:numS
        counter1.Send;
        [mz, spectrum,xy{n}] = readraw2spec(filename,n);
        raw_specs{n} = cat(1,mz,spectrum);
    end
    counter1.Destroy
    
    xy = cell2mat(xy);
    pixel = abs(xy(1)-xy(3));
    x_len = round(abs(xy(end-1)-xy(1))/pixel)+1;
    y_len = round(abs(xy(end)-xy(2))/pixel)+1;

    dims = [x_len,y_len];
    
    if length(raw_specs) ~= x_len*y_len
        disp('Warning! The pixel dimension according to the metadata does not match the file size.')
    end
end