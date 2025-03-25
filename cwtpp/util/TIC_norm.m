function [specs_norm] = TIC_norm(specs_input)
% Original also highly bloated.

% This is quicker by typically >1000 fold (typically ~ 3 seconds), less
% when the wait bar is turned off
% tic;
specs_norm = bsxfun(@rdivide,specs_input,sum(specs_input,2));
specs_norm(isnan(specs_norm)|isinf(specs_norm))=0;
% toc;

% Benchmark cmparison below
% tic;
% [a] = original(specs_input);
% toc;
% 
% 
% df = bsxfun(@minus,specs_norm,a);
% df = bsxfun(@rdivide,df,specs_input) * 100;
% 
% disp(['Largest difference is ' sprintf('%0.6f',max(df(:))) '%']);


end


function [specs_norm] = original(specs_input)

%counter1 = parfor_wait(size(specs_input,1), 'Waitbar', true);
parfor i = 1:size(specs_input,1)
    %counter1.Send;
    if sum(specs_input(i,:))~= 0
        specs_norm(i,:) = specs_input(i,:)/sum(specs_input(i,:));
    else
        specs_norm(i,:) = zeros(1,length(specs_input(i,:)));
    end
end

%counter1.Destroy

end