function TIC_image = TIC(raw_specs, newdims)
TIC_image = [];
    tic
    counter1 = waitbar(0,'Summing...');
    if isa(raw_specs,'cell') == 1
        for j=1:length(raw_specs)
            waitbar(j/length(raw_specs), counter1);
            dum = double(cell2mat(raw_specs(j)));
            TIC_image(j) = sum(dum(2,:));
        end
    else
        dims = size(raw_specs);
        for j=1:dims(1)
            waitbar(j/length(raw_specs), counter1);
            TIC_image(j) = sum(raw_specs(j,:));
        end
    end
    close(counter1);
    TIC_image = reshape(TIC_image, newdims(1), newdims(2));
    figure;imagesc(fliplr(rot90(TIC_image, -1)));colormap(turbo(256));colorbar;
    toc
end