function [list_of_peaks] = pp_imzml(raw_specs, dims,options)
     arguments
        raw_specs
        dims
        options.T0 {mustBeNumeric} = 1000;
        options.ite {mustBeNumeric} = 100;
        options.Nparticles {mustBeNumeric} = 1/2;
        options.Imin {mustBeNumeric} = 100;
        options.G {mustBeNumeric} = 3;
        options.L {mustBeNumeric} = 2;
    end  
    numS = dims(1)*dims(2);

    list_of_peaks = cell(numS,1);
    
    counter1 = parfor_wait(numS, 'Waitbar', false);
    parfor n = 1:numS
        counter1.Send;
        spectrum = raw_specs{n}(2,:);
        mz = raw_specs{n}(1,:);
        list_of_peaks{n} = ridgedetect(mz,spectrum,T0 = options.T0, ite = options.ite, ...
           Nparticles = options.Nparticles, Imin=options.Imin, G= options.G, L=options.L);
    end
    counter1.Destroy
end