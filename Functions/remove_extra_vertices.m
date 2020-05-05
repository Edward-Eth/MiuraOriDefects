function [fcn_output] = remove_extra_vertices(fcn_input)
%%
flag = zeros(1,length(fcn_input));
for ii = 1:size(fcn_input,1)-1  
    for jj = ii+1:size(fcn_input,1)    
        if abs(norm(fcn_input(ii,:) - fcn_input(jj,:))) < 1e-6          
            flag(jj) = 1;           
        end
    end
end
fcn_output = fcn_input(flag == 0,:);