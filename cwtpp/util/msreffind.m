function [mzs_recal] = msreffind(varargin)
    datasets = varargin{1};
    for i = 1:length(datasets)
       lens(i) = size(datasets{i},2);
    end
    [~,I] = sort(lens,'descend');
    datasets = datasets(I);
    mzs = mzs(I);

    mz_recal = mzs{1};
    mzs_recal = mzs;
    ms_references = [];
    [~,locs] = findpeaks(mean(datasets{1}),'MinPeakProminence',mean(mean(datasets{1})));
    mz_recal_p = mz_recal(locs);
    ms_references(:,1) = mz_recal_p;
    % ms_references(:,1) = mz_recal;

    if nargin < 2
        threshold = 200;
    end

    mzs_recal = {};
    for n = 2:length(mzs)
        disp(n)
        mz_raw = mzs{n};
        [~,locs] = findpeaks(mean(datasets{n}),'MinPeakProminence',mean(mean(datasets{n})));
        mz_raw_p = mz_raw(locs);
    %     [mz_new_p] = MSrecal(mz_raw_p,references, 1000);

    %     mz_raw(locs) = mz_new_p;
    %     [mz_new] = MSrecal(mz_raw,references, 1000);
    %     mzs_recal{n} = mz_new;
        for m = 1:length(ms_references)
                [diff, ind] = min( abs(mz_raw_p-ms_references(m,1)) );
    %             [diff, ind] = min( abs(mz_raw-ms_references(m,1)) );
                ppm = diff/ms_references(m,1)*10^6;
                if ppm <= threshold
                    ms_references(m,n) = ind;
                else
                    ms_references(m,:)=0;
                end
        end     
    end
    ms_references(ms_references(:,1)==0,:)=[];

    for i = 2:length(mzs_recal)
        dum = mzs_recal{i};
        dum(ms_references(:,i)) = ms_references(:,1);
        mzs_recal{i} = dum;
        clear dum
    end
end