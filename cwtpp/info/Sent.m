function [s_entropy] = Sent(spectrum)
    spectrum=spectrum/sum(spectrum);
    spectrum = spectrum(spectrum~=0);
    s_entropy=-(sum(spectrum.*(log2(spectrum))));
end