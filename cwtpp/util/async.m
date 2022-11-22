%% Asymmetric synchronisation after DTW

function x1 = async (spec_i,i1,i2)
    s = unique(i2);
    x1 = zeros(1,length(s));
    for i = 1:length(s)
        dum = spec_i(i1);
        x_mean = mean(dum(i2==s(i)));
        x1(1,i) = x_mean;
    end
end
