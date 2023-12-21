function [list_of_peaks,dims, mode] =  cwtpp(varargin)
tic

if isempty(varargin)
    [path] = uigetdir();
    cd (path);
    files=[dir('*.raw');dir('*.imzml')];
    filename=files(1).name;
    [~,~,ext] = fileparts(filename);
    if strcmp(ext, '.raw')
        mode = '.raw';
    else
        mode = '.imzml';
    end
elseif ischar(varargin{1})
    filename = varargin{1};
    [~,~,ext] = fileparts(filename);
    if strcmp(ext, '.raw')
        mode = '.raw';
    else
        mode = '.imzml';
    end
else
    filename = varargin{1};
    mode = 'workspace';
end

if nargin > 1
    Imin = varargin{2};
else
    Imin = 100;
end

if nargin > 2
    use_metadata = varargin{3};
else
    use_metadata = 0;
end
    
[list_of_peaks,dims] = pp_mode(filename,mode,Imin);
if use_metadata == 1
    mode = '.raw';
else
    if strcmp(mode,'workspace') == 0
        if strcmp(mode,'.raw')
            save([filename,'\cwtpeaks'],'list_of_peaks','dims')
        else
            save(['cwtpeaks'],'list_of_peaks','dims')
        end
    end
end


toc
