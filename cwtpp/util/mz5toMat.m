function [raw_specs,dims,chromoT,chromoI,features] = mz5toMat(filename, options)

arguments
    filename string = []
    options.mode = []
    % options.param = 1
end

if isempty(filename)
    [file,dir] = uigetfile('*.mz5');
    filename = [dir,'\',file];
end


if strcmp(options.mode,'TQ')
    chromoT=h5read(filename,['/ChomatogramTime']);
    Sindex = double(h5read(filename,['/ChromatogramIndex']));
    chromoI=double(h5read(filename,['/ChromatogramIntensity']));
    Sintensity = chromoI;
    Sintensity= reshape(Sintensity,[Sindex(1),length(Sindex)]);
    mz = h5read(filename,['/ChromatogramList']).id;
    raw_specs = Sintensity(:,2:end); % remove the TIC line
    features = mz(2:end); % remove the TIC line
else
    chromoT=h5read(filename,['/ChomatogramTime']);
    chromoI=h5read(filename,['/ChromatogramIntensity']);
    Sindex = double(h5read(filename,['/SpectrumIndex']));
    if find(round(Sindex/1e9,4)==round(2^32/1e9,4)) % to correct for 32-bit overflow
        disp('overflow in Sindex detected. Correcting...')
        f_ind = find(round(Sindex/1e9,4)==round(2^32/1e9,4));
        Sindex(f_ind+1:end) = Sindex(f_ind+1:end)+2^32;
    end
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
end

[~,locs] = findpeaks(diff(chromoT),'MinPeakProminence',mean(diff(chromoT)));
if length(unique(diff(locs)))>1
    disp('non-discrete image dimension detected.Maybe not an image.')
    [raw_specs, dims] = guessroi(raw_specs,chromoT);
    disp(['we are guessing that the image dimension is: ', num2str(dims(1)),' by ', num2str(dims(2))])
else
    dims = [unique(diff(locs)),length(Sindex)/unique(diff(locs))];
end

%% function to trim a line at the bottom of image due to missing pixels
function [raw_specs_new, est_dims] = guessroi(raw_specs,chromoT)
    [~,locs] = findpeaks(diff(chromoT),'MinPeakProminence',mean(diff(chromoT)));
    x_dims = min(unique(diff(locs)));
    y_dims = floor(size(raw_specs,1)/x_dims);
    est_dims = [x_dims,y_dims];
    raw_specs_new = raw_specs(1:est_dims(1)*est_dims(2),:);
end

end