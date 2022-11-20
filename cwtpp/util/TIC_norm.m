function [specs_norm] = TIC_norm(specs_input)
    counter1 = parfor_wait(size(specs_input,1), 'Waitbar', true);
    parfor i = 1:size(specs_input,1)
        counter1.Send;
        if sum(specs_input(i,:))~= 0
            specs_norm(i,:) = specs_input(i,:)/sum(specs_input(i,:));
        else
            specs_norm(i,:) = zeros(1,length(specs_input(i,:)));
        end
    end
    counter1.Destroy
end