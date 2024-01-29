mzs_recal = {};
ppms = [];
threshold = 300;%ppm

for i = 1:length(mzs)
% for i = 1:5
    disp(i)
    [mzs_recal{i},ppm] = MSrecal(mzs{i},references, threshold);
    ppms(i) = mean(ppm);
end

select = (ppms<=threshold);
mzs_recal = mzs_recal(select);
