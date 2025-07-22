function [raw_specs,dims,chromoT,chromoI] = mz5toMat(varargin)

% arguments
    % filename string
    % options.mod = []
    % options.param = 1
% end

if isempty(varargin)
    [file,dir] = uigetfile('*.mz5');
    filename = [dir,'\',file];
else
    filename = varargin{1};
end

chromoT=h5read(filename,['/ChomatogramTime']);
chromoI=h5read(filename,['/ChromatogramIntensity']);
Sindex = h5read(filename,['/SpectrumIndex']);

raw_specs = {};
tic
parfor n = 1:length(Sindex)
    % disp(n)
    if n == 1
        spec = h5read(filename,['/SpectrumIntensity'],1,double(Sindex(n)));
        mz = h5read(filename,['/SpectrumMZ'],1,double(Sindex(n)));
    else
        try
            spec = h5read(filename,['/SpectrumIntensity'],double(Sindex(n-1))+1,double(Sindex(n)-Sindex(n-1)));
            mz = h5read(filename,['/SpectrumMZ'],double(Sindex(n-1))+1,double(Sindex(n)-Sindex(n-1)));
        catch
            spec = h5read(filename,['/SpectrumIntensity'],double(Sindex(n-1))+1,1);
            mz = h5read(filename,['/SpectrumMZ'],double(Sindex(n-1))+1,1);
        end
    end
    for i=2:length(mz)
        mz(i) = mz(i-1)+mz(i);
    end
    raw_specs{n} = cat(1,mz',spec');
end
toc


[~,locs] = findpeaks(diff(chromoT),'MinPeakProminence',mean(diff(chromoT)));
if floor(length(Sindex)/unique(diff(locs))) ~= length(Sindex)/unique(diff(locs))
    disp('non-discrete image dimension detected.Maybe not an image.')
    dims = [0,0];
else
    dims = [unique(diff(locs)),length(Sindex)/unique(diff(locs))];
end