function [filt,thresh] = Sent_filter(spectrum,term)

% mz = mzs{1};
% spectrum = mean(datasets{1});
% figure,stem(mz,spectrum,'Marker','none')

% b = [0,50,100,200,500,1000,5000,1e4,1e5,1e6,1e7,1e8];
b = 0:10:5000;
ent = [];
for n = 1:length(b)
    dum = (spectrum-b(n));dum(dum<0|dum==0)=1;ent(n,:) = dum;
    S(n)=Sent(dum);
end
% figure,plot(log10(b),S)

% [~,thresh] = min(S);
% figure
% [FitResults1,GOF1,~,coefs]=peakfit([log10(b(thresh-2:thresh+2))' ...
    % S(thresh-2:thresh+2)'], ...
    % 0,0,1,28,2,1,0,0);
% base = 10^(-coefs(2)/(2*coefs(1)));
% base = 10^((4*coefs(1)*coefs(3)-coefs(2)^2)/(4*coefs(1)));
if nargin < 2
    term = 0.001;
end
base = find(diff(gradient(S))<term);
thresh = b(base(1)+1);

filt = (spectrum<thresh);
% spectrum_f = spectrum;spectrum_f(filt) = [];
% mz_f = mz;mz_f(filt) = [];
% figure,stem(mz_f,spectrum_f,'Marker','none')