function [Pbin]=PijrNZbin(Event,trueRGM)

% This function computes the exceedance probability P(y_ij?Y_ir ) for each 
% ground motion map j, site i, and return period r, where Y_ir is the 
% ground shaking from the “true” hazard curve for return period r at site i, 
% and y_ij is the ground shaking at site i caused by ground motion map j.

% Probabilities are all 0 or 1 since there's no uncertainty in ground
% motion for a map

sitenum=length(Event(1).GM(:,1));   % Number of sites
count=length(Event);                % Number of ground motion maps

for i=1:count
    for j=1:sitenum
        for r=1:length(trueRGM(1,:))
            if Event(1,i).GM(j,1)>trueRGM(j,r)
                % If ground motion y_ij>Y_ir, assign 1
                Pbin(i).gm(j,r)=1
            else % If ground motion y_ij<=Y_ir, assign 1
                Pbin(i).gm(j,r)=0
            end
        end
    end
end

end