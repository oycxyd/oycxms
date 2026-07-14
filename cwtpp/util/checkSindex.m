function Sindex_fixed = checkSindex(Sindex)
    if find(round(Sindex/1e9,4)==round(2^32/1e9,4)) % to correct for 32-bit overflow
        f_ind = find(abs(round(Sindex/1e9,4)-round(2^32/1e9,4))<1e-4);
        f_ind((diff(f_ind)<10)) = [];
        disp(['overflow in Sindex detected. Correcting ',num2str(length(f_ind)),' discontinuities...'])
        Sindex_fixed = Sindex;
        for i = 1:length(f_ind)
            Sindex_fixed(f_ind(i)+1:end) = Sindex_fixed(f_ind(i)+1:end)+2^32;
        end
        % figure,plot(Sindex)
    else
        Sindex_fixed = Sindex;
    end
end