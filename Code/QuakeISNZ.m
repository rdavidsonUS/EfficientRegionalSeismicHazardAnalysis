function M=QuakeISNZ(faults,backgrounds,distMchar)

% This function uses Monte Carlo simulation with importance sampling of the 
% magnitudes to simulate a set of earthquake scenarios. 
% For background sources, sample magnitude following Gutenberg-Richter
%     relationship
% For faults, 
%     distMchar=1 means a truncated magnitude distribution function is used
%     distMchar=0 means a delta function is used (just one possible magnitude, Mchar)
% For all sources, a bin size of 0.1 is assumed for sampling from the 
%    magnitude distributions.

count=1;

% Sample earthquakes from fault sources
for f=1:length(faults)
    if distMchar ==0
        %Delta function,only one earthquake per fault
        M(count,1)=count;
        M(count,2)=f;
        M(count,3)=faults(f).Mchar;
        M(count,4)=1;
        M(count,5)=1/faults(f).Tyear;
        count=count+1;
    else
        abin=faults(f).Mchar-0.2:0.1:faults(f).Mchar+0.2;
        for j=1:length(abin)-1
            M(count,1)=count;
            M(count,2)=f;
            M(count,3)=abin(j)+0.05;
            M(count,4)=1;
            F=@(x) exp(-(x-faults(f).Mchar).^2/(2*0.06^2))/(0.06*sqrt(2*pi));
            M(count,5)=(1-exp(-1/faults(f).Tyear))*quadl(F,abin(j),abin(j+1));
            count=count+1;
        end      
    end   
    %disp(['faultnumber',num2str(f)]);
end

% Sample earthquakes from background sources
for i=1:length(backgrounds)
    ara(i)=backgrounds(i).avalue;
    arb(i)=backgrounds(i).bvalue;
    Mmina(i)=backgrounds(i).Mmin;
    Mmaxa(i)=backgrounds(i).Mmax;
    abin=Mmina(i):0.1:Mmaxa(i); % divide the magnitude range of i into partitions by 0.1
    for j=1:length(abin)-1
        M(count,1)=count;           % Scenario ID number
        M(count,2)=i;               % ID of source
        M(count,3)=abin(j)+0.05;    % Magnitude of eq
        M(count,4)=0;               % Indicator of background (0) or fault (1) source
        M(count,5)= exp(-(10^(ara(i)-arb(i)*abin(j+1))))-exp(-(10^(ara(i)-arb(i)*abin(j))));    % Annual occurrence probability
        count=count+1;
    end
    %disp(['background seismicity',num2str(i)]);
end

end