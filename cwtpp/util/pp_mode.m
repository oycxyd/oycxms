%% new preprocessing workflow
%% switch function for different file formats/modes

function [list_of_peaks,dims] = pp_mode(filename, options)
    arguments
        filename
        options.mode string = '.raw';
        options.T0 {mustBeNumeric} = 1000;
        options.ite {mustBeNumeric} = 100;
        options.Nparticles {mustBeNumeric} = 1/2;
        options.Imin {mustBeNumeric} = 100;
        options.G {mustBeNumeric} = 3;
        options.L {mustBeNumeric} = 2;
    end
    

end
