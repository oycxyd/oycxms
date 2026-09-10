function [cwtpp_params] = infogsearch(filename,nsamples,meta_switch)   
    if nargin < 2
        nsamples = 5;
        meta_switch = 0;
    end
    %% estimate baseline using spectral entropy
    try
        Sindex = double(h5read(filename,['/SpectrumIndex']));
        Sindex = checkSindex(Sindex);
        [sample_ind] = find(diff(Sindex)==round(median(diff(Sindex))));
    catch
        med = median(cellfun(@length, filename));
        [sample_ind] = find((cellfun(@length, filename))==round(med));
        meta_switch = 1;
    end
    if length(sample_ind) < nsamples
        nsamples = length(sample_ind);
    end

    disp('optimising hyperparameters...')
    tic
    for m = 1:nsamples
        ind = sample_ind(m);
        if meta_switch
            mz = filename{m}(1,:);
            spectrum = filename{m}(2,:);
        else
            ind = ind+1;
            spectrum = h5read(filename,['/SpectrumIntensity'],double(Sindex(ind(1)-1))+1,double(Sindex(ind(1))-Sindex(ind(1)-1)));
            mz = h5read(filename,['/SpectrumMZ'],double(Sindex(ind(1)-1))+1,double(Sindex(ind(1))-Sindex(ind(1)-1)));
            mz = mz';spectrum = spectrum';
            for i=2:length(mz)
                mz(i) = mz(i-1)+mz(i);
            end
        end
        % figure,plot(mz,spectrum)
        
        % tic
        b = 0:max(spectrum)/500:max(spectrum);
        % b = 0:500;
        S = [];
        for n = 1:length(b)
            sub = spectrum-b(n);sub(sub<0)=0;
            S(n) = Sent(sub);
        end
        % figure,plot((b),gradient(S))
        % figure,plot((b(2:end)),diff(gradient(S)))
        cutoff = find(abs(diff(gradient(S)))<0.01);
        cutoff = b(cutoff(1)+1);
        Starget = Sent(spectrum(spectrum>cutoff));
        % toc

    %% grid search the rest and constrain to baseline-subtracted spectral info.
        Nparticles_range = [1/4, 1/2, 1];
        G_range = 1:5;
        L_range = 1:5;
        [A,B,C] = ndgrid(Nparticles_range,G_range,L_range);
        params = [A(:), B(:), C(:)];
        
        nComb   = size(params,1);
        obj  = zeros(nComb,1);
        cwttest{1} = cat(1,mz,spectrum);
        
        parfor n = 1:nComb
            % disp(n)
            p = params(n,:);
            [cwtpeaks]=cwtpp(cwttest,...
                G=(p(2)),L=(p(3)),...
                Nparticles=(p(1)),verbose=false);
            [~,idxInY] = ismember(cwtpeaks{1}(1,:),mz);
            obj(n) = (Sent(cwtpeaks{1}(2,:))-Starget)...
            +mean(mean(spectrum(idxInY)-cwtpeaks{1}(2,:)));
        end
        [~,opt] = min(abs(obj));
        opt_params(:,m) = cat(2,cutoff,params(opt,:));
    end
    cwtpp_params(1) = mean(opt_params(1,:));
    cwtpp_params(2) = median(opt_params(2,:));
    cwtpp_params(3) = median(opt_params(3,:));
    cwtpp_params(4) = median(opt_params(4,:));
    toc
end