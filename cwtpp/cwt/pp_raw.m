function [list_of_peaks, dims] = pp_raw(filename, options)
     arguments
        filename
        options.T0 {mustBeNumeric} = 1000;
        options.ite {mustBeNumeric} = 100;
        options.Nparticles {mustBeNumeric} = 1/2;
        options.Imin {mustBeNumeric} = 100;
        options.G {mustBeNumeric} = 3;
        options.L {mustBeNumeric} = 2;
    end        
    [~,p2] = watersPackages(filename);
    numS = calllib('MassLynxRaw','getScansInFunction',p2,1);
    list_of_peaks = cell(numS,1);
    
    counter1 = parfor_wait(numS, 'Waitbar', false);
    for n = 1:numS
        counter1.Send;
        [mz,spectrum,xy] = readraw2spec(filename,n);
        dims_x(n) = xy(1);
        dims_y(n) = xy(2);
        list_of_peaks{n} = ridgedetect(mz,spectrum,T0 = options.T0, ite = options.ite, ...
           Nparticles = options.Nparticles, Imin=options.Imin, G= options.G, L=options.L);       
    end

    try
        dims = get2Dcoord((dims_x),(dims_y));
        dims = size(dims);
        dims = fliplr(dims);
    catch
        warning('Cannot obtain XY dims from raw file; check if data is an image!');
        [list_of_peaks, dims] = guessroi(list_of_peaks,dims_x);
    end
    counter1.Destroy
end