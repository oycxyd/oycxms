function [list_of_peaks,dims, mode] =  cwtpp(filename,options)

arguments
    filename = 'default'
    options.Imin = 100
    options.use_metadata = 0
    options.dir = ''
end

% tic
if strcmp(filename,'default')
    [path] = uigetdir();
    cd (path);
    files=[dir('*.raw');dir('*.imzml');dir('*.mz5')];
    filename=files(1).name;
end
try
    [~,file,ext] = fileparts(filename);
    if strcmp(ext, '.raw')
        mode = '.raw';
    elseif strcmp(ext, '.imzML')
        mode = '.imzml';
    else
        mode = '.mz5';
    end
catch
    mode = 'workspace';
end

[list_of_peaks,dims] = pp_mode(filename,mode=mode,Imin=options.Imin);
if options.use_metadata == 1
    mode = '.raw';
    save([options.dir,'\cwtpeaks'],'list_of_peaks','dims')
else
    if strcmp(mode,'workspace') == 0
        if strcmp(mode,'.raw')
            save([filename,'\cwtpeaks'],'list_of_peaks','dims')
        else

            save([file,'_cwtpeaks'],'list_of_peaks','dims')
        end
    end
end
% toc
