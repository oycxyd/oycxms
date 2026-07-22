function [list_of_peaks] = ridgedetect(mz,spectrum,options)
    arguments
        mz
        spectrum
        options.T0 {mustBeNumeric} = 1000;
        options.ite {mustBeNumeric} = 100;
        options.Nparticles {mustBeNumeric} = 1/2;
        options.Imin {mustBeNumeric} = 100;
        options.G {mustBeNumeric} = 3;
        options.L {mustBeNumeric} = 2;
    end        
if length(spectrum)< 5
    list_of_peaks = [0; 0];
else
%         scales = 1:2:32;
    scales = 1:2:16;
    cwttest = cwtft(spectrum,'scales',scales,'wavelet','mexh');% need to figure out scale for CWT
    wcoefs = cwttest.cfs;

%% crazy-climber/other max. finding
    [rows,cols] = size(wcoefs);
    % if isempty(options.Nparticles)
    %     disp('default Nparticles')
    %     Nparticles = floor(cols * rows/2); % number recommended by Zheng et al.
    % else
    %     Nparticles = options.Nparticles;
    % end
    Nparticles = floor(cols * rows*options.Nparticles);
    Tf = 0;
    T = [options.T0:(Tf-options.T0)/(options.ite-1):Tf];
    
    % its_per_stage = rows*4;
    % T0 = 500;
    % ite = 100;
    % T = 500*[ones(1, its_per_stage)*1 ones(1, its_per_stage)*.1 ones(1, its_per_stage)*.01 ones(1, its_per_stage)*.001]; 
    % not sure what this does, need to read about simulated annealing
    % T = [500 400 300 200 100];
    % T = [1000 800 600 400 200];
    if options.ite > 0
        ridgetable = CrazyClimber((wcoefs), T, Nparticles);
        % figure,imagesc(ridgetable);colorMap = jet(256);colormap(colorMap); colorbar;
        ridgetable = flipud(-ridgetable);
        try
            ridgetable = round(spectrum.*(ridgetable));
        catch
            ridgetable = round(spectrum'.*(ridgetable));
        end
    else
        ridgetable = flipud(wcoefs);
    end
    % ridgetable(ridgetable<0)=0;

        % find local maxima for all levels
    maxima = {};
    for i = 1:length(scales)
        try
    %     [pks,locs] = findpeaks(ridgetable(i,:),'MinPeakHeight',Imin,'MinPeakDistance',3);
            [pks,locs] = findpeaks(ridgetable(i,:),'MinPeakProminence',options.Imin,'MinPeakDistance',2);
            maxima{i} = cat(1,pks,locs);
        catch
            pks = 0; locs = 0;
            maxima{i} = cat(1,pks,locs);
        end
    end

    %% detect peaks from ridges
    % initiliase ridges
    ridges = {};
    ridges_ini = maxima{1};
    if isempty(ridges_ini)
        disp('no peak detected.')
        try
            list_of_peaks = [mz; zeros(1,length(spectrum))];
        catch
            list_of_peaks = [mz'; zeros(1,length(spectrum))];
        end
    else
        for i = 1:size(ridges_ini,2)
            ridges{i} = ridges_ini(:,i);
        end

        % search listing
        for i=2:length(scales)-1
            indices = [];
            maxima_i = maxima{i};
            for j=1:length(ridges)
                ridge = ridges{j};
                [val,ind]=min(abs(maxima_i(2,:)-ridge(2,end)));
                if abs(val)<=3
                    ridges{j} = cat(2,ridge,maxima_i(:,ind));
                    indices = [indices maxima_i(2,ind)];
                else
                    ridges{j} = cat(2,ridge,zeros(2,1));
                end
            end
            [~,indices] = setdiff(maxima_i(2,:),indices);
            new_ridges = maxima_i(:,indices);
            if length(new_ridges) == 2
                ridges{end+1} = (new_ridges);
            else
                for p = 1:length(new_ridges)
                    ridges{end+1} = (new_ridges(:,p));
                end
            end
        %     G = G+1;
        end

        % filtering & generate peak list
        G = options.G; %maximum gap between levels
        L = options.L; %minimun length of ridge
        peaks = [];
        for i = 1:length(ridges)
            ridge = ridges{i};
            if length(ridge(ridge == 0))>= G*2 || length(ridge) <= L
                ridges{i} = [];
            else
                [~,max_ind] = max(abs(ridge(1,:)));
                peak_ind = ridge(2,max_ind);
                try
                    while spectrum(peak_ind)<spectrum(peak_ind+1)
                        peak_ind = peak_ind+1;
                    end
                    while spectrum(peak_ind)<spectrum(peak_ind-1)
                        peak_ind = peak_ind-1;
                    end
                catch
                    continue
                end
                peaks(:,i) = [mz(peak_ind); spectrum(peak_ind)];
            end
        end
        if isempty(peaks)
            try
                list_of_peaks = [mz; zeros(1,length(spectrum))];
            catch
                list_of_peaks = [mz'; zeros(1,length(spectrum))];
            end
        else
            [~,ia] = unique(peaks(1,:));
            peaks = peaks(:,ia);
            peaks(:,1)=[];
            list_of_peaks = peaks;
        end
    end
end