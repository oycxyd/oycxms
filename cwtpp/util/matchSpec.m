function [cmz,aligned_peaks, range, shifts] = matchSpec(cwtpeaks, thresh, options)
% using cwt2cmz (James McKenzie) to find a common m/z axis 
% and match intensities from all scans/pixels
arguments
    cwtpeaks cell
    thresh (1,:) {mustBeNumeric,mustBeReal} = []
    options.freq (1,:) {mustBeNumeric,mustBeReal} = 0.05
end
%% load in detected peaks & initialise parameters
% mzs = cellfun(@(x) x(1,:), cwtpeaks, 'UniformOutput', false);
% datasets = cellfun(@(x) x(2,:), cwtpeaks, 'UniformOutput', false);
% tic
    if isempty(thresh)
        shifts = [];
        % quick look at number of common peaks that are found on >5% of data within
        % various ppm windows (5-800)
        range = [5 10 20 30 50 100 150 200 300 500 600 700 800];
        for i = 1:length(range)
            [pks] = cwt2cmz(cwtpeaks,'ppm',range(i),'mzRange',[50 1200], ...
                'minFreq',round(length(cwtpeaks)*0.05));
            if isempty(pks)
                pks = 0;
            end
            shifts(i) = length(pks);
        end
        % figure,plot(range,shifts)
        
        % define ppm threshold as that at which max. of peaks is detected (by
        % cwtpp)
        [thresh] = knee_pt(shifts,range);
    end
% [~,ind] = min(shifts - max(cellfun(@length, cwtpeaks)));[thresh] = 50*ind;
disp(['matching peaks using an estimated ppm of ',num2str(thresh)])
[cmz] = cwt2cmz(cwtpeaks,'ppm',thresh,'mzRange',[50 1200], ...
    'minFreq',round(length(cwtpeaks)*options.freq));
if isempty(cmz) cmz = 0; end
%% match peaks from all pixels
aligned_peaks = [];
M = length(cwtpeaks);N = length(cmz);
multi = {};
parfor m = 1:M
    multi_check = 0;
    for n = 1:N
        dum = cwtpeaks{m};
        [inds]=(abs(dum(1,:)-cmz(n)))/cmz(n)*1e6 <= thresh;
        if isempty(find(inds,1)) == 0
            if (length(find(inds)))>1
                multi_check = cat(2,multi_check,cmz(n));
                % [~,rank] = min(abs(dum(1,inds)-cmz(n)));
                % inds = inds(rank);
            end
            aligned_peaks(m,n) = sum(dum(2,inds));
        else
            aligned_peaks(m,n) = 0;
        end
    end
    if length(multi_check)>1
        multi{m} = multi_check;
    end
end