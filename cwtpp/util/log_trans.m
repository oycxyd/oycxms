function [Data_log] = log_trans(data)

tmp_log = data(data ~= 0);
logOS = nanmedian(tmp_log);
Data_log = log10(data+logOS);

end