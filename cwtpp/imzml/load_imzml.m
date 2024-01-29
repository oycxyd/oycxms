%% load imzml imaging files
function [raw_specs,dims,imzML] = load_imzml(varargin)
    if isempty(varargin)
        [filename,path] = uigetfile('*.imzml');
        filename = fullfile(path,filename);
    else
        filename = fullfile(pwd,varargin{1});
    end

    try
        javaclasspath('D:\BOX\Box Sync\RA\codes\Imaging\ambient MSI\oycxms\cwtpp\imzml\imzML converter\imzMLConverter.jar');
        imzML = imzMLConverter.ImzMLHandler.parseimzML(filename);
%         dims = [imzML.getWidth(),imzML.getHeight()];
        dims = [imzML.getHeight(),imzML.getWidth()];

        raw_specs = {};
        i = 1;
        for p = 1:dims(1)
            for q = 1:dims(2)
                spec = imzML.getSpectrum(p,q);
                if isempty(spec)
                    mz_i = zeros(10,1);
                    int_i = zeros(10,1);
                else
                    mz_i = spec.getmzArray();
                    int_i = spec.getIntensityArray();
                end
                raw_specs{i} = cat(1,mz_i',int_i');
                i = i+1;
            end
        end
    catch
        warning('sth went wrong!check you selected a .imzml file');
    end
end