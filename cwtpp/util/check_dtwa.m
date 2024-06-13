function [mean_specs] = check_dtwa(varargin)
    if isempty(varargin)
        dname = uigetdir();
        cd (dname);
        filenames=[dir('*.raw');dir('*.imzml')];
    end
    mean_specs = [];
    parfor n = 1:length(filenames)
        disp(['loading mean spec from file ',num2str(n),'/',num2str(length(filenames))])
        filename = filenames(n).name;
        if n == 1
            [data,~,mz]=h5toMat([filename,'\datacube_aligned.h5']);
            mean_specs(n,:) = mean(data);
        else
            [data]=h5toMat([filename,'\datacube_aligned.h5']);
            mean_specs(n,:) = mean(data);
        end
    end
    figure,imagesc(mean_specs,'XData',mz);xlabel('m/z')
end