%% Function to calculate Mander's coefficients (co-occurence)
function [M1, M2, MOC]=manders(image1, image2)
    image1 = image1(:);
    image2 = image2(:);
    xi = [];
    yi = [];

    for i = 1:length(image1)
        if image2(i) > 0
            xi(i) = image1(i);
        elseif image2(i) == 0
            xi(i) = 0;
        end
        if image1(i) > 0
            yi(i) = image2(i);
        elseif image1(i) == 0
            yi(i) = 0;
        end     
    end
    M1 = sum(xi)/sum(image1);
    M2 = sum(yi)/sum(image2);
    MOC = sum(image1.*image2)/(sqrt(sum(image1.^2))*sqrt(sum(image2.^2)));
end