function image = raw_plotter(ion,filename)
    [~,~,ext] = fileparts(filename);

    if strcmp(ext,'.raw')
        [~,p2] = watersPackages(filename);
        numS = calllib('MassLynxRaw','getScansInFunction',p2,1);
        image = zeros(numS, 1);
        disp('loading from raw file')
        parfor n = 1:numS
            [mz,spectrum,xy] = readraw2spec(filename,n);
            dims_x(n) = xy(1);
            dims_y(n) = xy(2);
            try
                [~,ind] = min(abs(mz-ion));
                image(n) = spectrum(ind);
            catch
                image(n) = sum(spectrum(:));
            end
        end
        dims = get2Dcoord((dims_x),(dims_y));
        dims = size(dims);
        dims = fliplr(dims);
    elseif strcmp(ext,'.mz5')
        disp('loading from mz5 file')
        Sindex = double(h5read(filename,['/SpectrumIndex']));
        image = zeros(length(Sindex), 1);
        if find(round(Sindex/1e9,4)==round(2^32/1e9,4)) % to correct for 32-bit overflow
        f_ind = find(abs(round(Sindex/1e9,4)-round(2^32/1e9,4))<1e-4);
        f_ind((diff(f_ind)<10)) = [];
        disp(['overflow in Sindex detected. Correcting ',num2str(length(f_ind)),' discontinuities...'])
            for i = 1:length(f_ind)
                Sindex(f_ind(i)+1:end) = Sindex(f_ind(i)+1:end)+2^32;
            end
            % figure,plot(Sindex)
        end
        chromoT=h5read(filename,['/ChomatogramTime']);
    
        parfor n = 1:length(Sindex)
            if n == 1
                spectrum = h5read(filename,['/SpectrumIntensity'],1,double(Sindex(n)));
                mz = h5read(filename,['/SpectrumMZ'],1,double(Sindex(n)));
            else
                try
                    spectrum = h5read(filename,['/SpectrumIntensity'],double(Sindex(n-1))+1,double(Sindex(n)-Sindex(n-1)));
                    mz = h5read(filename,['/SpectrumMZ'],double(Sindex(n-1))+1,double(Sindex(n)-Sindex(n-1)));
                catch
                    spectrum = h5read(filename,['/SpectrumIntensity'],double(Sindex(n-1))+1,1);
                    mz = h5read(filename,['/SpectrumMZ'],double(Sindex(n-1))+1,1);
                end
            end
            for i=2:length(mz)
                mz(i) = mz(i-1)+mz(i);
            end
            try
                [~,ind] = min(abs(mz-ion));
                image(n) = spectrum(ind);
            catch
                image(n) = sum(spectrum(:));
            end
        end
    
        [~,locs] = findpeaks(diff(chromoT),'MinPeakProminence',mean(diff(chromoT)));
        if length(unique(diff(locs)))>1 | mod(length(Sindex),unique(diff(locs)))>0
            disp('non-discrete image dimension detected.Maybe not an image.')
            [image, dims] = guessroi(image,chromoT);
            disp(['we are guessing that the image dimension is: ', num2str(dims(1)),' by ', num2str(dims(2))])
        else
            dims = [unique(diff(locs)),length(Sindex)/unique(diff(locs))];
        end
    end
    image = reshape(image,dims);
    figure,imagesc(image);axis image;colormap('magma')
end

