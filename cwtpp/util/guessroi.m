function [raw_specs_new, est_dims] = guessroi(raw_specs,chromoT)
%% function to guesstimate image dimensions and to correct those of imcomplete/glitched runs
    [~,locs] = findpeaks(diff(chromoT),'MinPeakProminence',mean(diff(chromoT)));
    x_dims = min(unique(diff(locs)));
    y_dims = floor(size(raw_specs,1)/x_dims);
    est_dims = [x_dims,y_dims];
    raw_specs_new = raw_specs(1:est_dims(1)*est_dims(2),:);
end