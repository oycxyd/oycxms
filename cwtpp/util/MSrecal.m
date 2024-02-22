function [new_mz,ppms2,coefs] = MSrecal(mz_raw,references, threshold)

    % threshold = 100;%in ppm

    % y = references;
    stop = 0;
    x = [];
    y = [];
    ppms1 = [];
    if length(references) < 2
        [diff] = min( abs(mz_raw-references(1)) );
        ppms1 = diff/references(1)*10^6;
        if ppms1 <= threshold
            new_mz = mz_raw+diff;
            ppms2 = ppms1;
            disp(['the average ppm after correction is ', num2str((ppms2))])
        else
            disp(['the lock mass chosen was not found within ', num2str(threshold),' ppm!'])
        end
    else
        while stop == 0
            for i = 1:length(references)
                [diff, ind] = min( abs(mz_raw-references(i)) );
            %     [diff, ind] = min( abs(new_mz-references(i)) );
                ppms1(i) = diff/references(i)*10^6;
                if ppms1(i) <= threshold
                    x = cat(1,x,mz_raw(ind));
                        y = cat(1,y,references(i));
                end
            end
    
                % ppms1 = ppms1';
            disp(['the average ppm before correction is ', num2str(mean(ppms1))])
    
            n = 3;
            [~, ~,~,coefs]=peakfit([x y],0,0,1,28,n,1,0,0,0,0);
    
            % [dist] = dtw(x,references);
            % mz_new = async (x,i1,i2);
            % 
            % [cow_params] = optim_cow(x,[1 3 1 1],[1 3 50 0.15],y);
            % [~,mz_new] = cow(y,x,1,0.1);
    
            for p = 1:n+1
                if p == 1
                    new_mz = coefs(end);
                else
                    new_mz = new_mz + coefs(end-p+1)*mz_raw.^(p-1);
                end
            end
    
            ppms2 = [];
            for i = 1:length(references)
            %     [diff, ind] = min( abs(mz_raw-references(i)) );
                [diff, ind] = min( abs(new_mz-references(i)) );
                ppms2(i) = diff/references(i)*10^6;
            end
            disp(['the average ppm after correction is ', num2str(mean(ppms2))])
    
            if mean(ppms2) < mean(ppms1)
                mz_raw = new_mz;
            else
                stop = 1;
                new_mz = mz_raw;
            end
        end
    end
end