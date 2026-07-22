function [list_of_peaks] = pp_workspace(raw_specs, options)
    arguments
        raw_specs
        options.T0 {mustBeNumeric} = 1000;
        options.ite {mustBeNumeric} = 100;
        options.Nparticles {mustBeNumeric} = 1/2;
        options.Imin {mustBeNumeric} = 100;
        options.G {mustBeNumeric} = 3;
        options.L {mustBeNumeric} = 2;
    end
    numS = length(raw_specs);
    
    list_of_peaks = cell(numS,1);
    
    counter1 = parfor_wait(numS, 'Waitbar', false);
    parfor n = 1:numS
        counter1.Send;
        spectrum = raw_specs{n}(2,:);
        mz = raw_specs{n}(1,:);
        list_of_peaks{n} = ridgedetect(mz,spectrum,T0 = options.T0, ite = options.ite, ...
           Nparticles = options.Nparticles, Imin=options.Imin, G= options.G, L=options.L);
    end
    % if exist('dims','var') == 0
    %     dims = [0,0];
    % end
    counter1.Destroy
end