 function [list_of_peaks, dims] = pp_mz5(filename, options)  
         arguments
            filename
            options.T0 {mustBeNumeric} = 1000;
            options.ite {mustBeNumeric} = 100;
            options.Nparticles {mustBeNumeric} = 1/2;
            options.Imin {mustBeNumeric} = 100;
            options.G {mustBeNumeric} = 3;
            options.L {mustBeNumeric} = 2;
        end       
        % Read and correct Sindex once
        Sindex = double(h5read(filename, '/SpectrumIndex'));
        Sindex = checkSindex(Sindex);
        chromoT = h5read(filename, '/ChomatogramTime');

        numS = length(Sindex);
        list_of_peaks = cell(numS,1);
        
        counter1 = parfor_wait(numS, 'Waitbar', false);
        parfor n = 1:numS
            % disp(n)
            counter1.Send;
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
                list_of_peaks{n} = ridgedetect(mz,spectrum,T0 = options.T0, ite = options.ite, ...
                Nparticles = options.Nparticles, Imin=options.Imin, G= options.G, L=options.L);
            catch
                warning(['error on spectrum# ',num2str(n)])
            end
        end
        
        [~,locs] = findpeaks(diff(chromoT),'MinPeakProminence',mean(diff(chromoT)));
        if length(unique(diff(locs)))>1 | mod(length(Sindex),unique(diff(locs)))>0
            disp('non-discrete image dimension detected.Maybe not an image.')
            [list_of_peaks, dims] = guessroi(list_of_peaks,chromoT);
            disp(['we are guessing that the image dimension is: ', num2str(dims(1)),' by ', num2str(dims(2))])
        else
            dims = [unique(diff(locs)),length(Sindex)/unique(diff(locs))];
        end
        counter1.Destroy
    end