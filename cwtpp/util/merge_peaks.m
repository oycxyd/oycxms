function [mz_f,data_f] = merge_peaks(mz,data)
test = diff(mz);
test = cat(2,1,test);
indc = (test>0.01);

mz_f = [];data_f = [];
counter1 = parfor_wait(length(indc), 'Waitbar', true);
for n = 1:length(indc)
    counter1.Send;
    if indc(n)
        mz_f = cat(1,mz_f,mz(n));
        data_f = cat(2,data_f,data(:,n));
    else
        mz_f(end) = (mz_f(end)+mz(n))/2;
        data_f(:,end) = data_f(:,end)+ data(:,n);
    end
end
counter1.Destroy
