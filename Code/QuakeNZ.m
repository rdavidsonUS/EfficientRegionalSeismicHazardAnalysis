function [M]=QuakeNZ(years,faults,backgrounds,distMchar)

% This function uses conventional Monte Carlo simulation (MCS) to simulate 
% a set of earthquake scenarios for a user-specified number of years.
% For background sources, sample magnitude following Gutenberg-Richter
%     relationship
% For faults, 
%     distMchar=1 means a truncated magnitude distribution function is used
%     distMchar=0 means a delta function is used (just one possible magnitude, Mchar)

count=1;

for yr=1:years          % Step through years
    % Sample earthquakes from fault sources
    for f=1:length(faults)
        randquake=rand(1);   % Generate random number to determing if EQ occurs and possible magnitude
        if randquake<=(1/faults(f).Tyear)  % Determine if EQ occurs in year yr
           if distMchar ==0    % Delta function --> only Mchar happens
               M(count,1)=yr;               % Year of eq 
               M(count,2)=f;                % ID of source
               M(count,3)=faults(f).Mchar;  % Magnitude of eq
               M(count,4)=1;                % Indicator of background (0) or fault (1) source
               M(count,5)=1/years;          % Annual occurrence probability
               count=count+1;
           else                % Truncated normal function for magnitude
               M(count,1)=yr;
               M(count,2)=f;
               M(count,3)=max(faults(f).Mchar-0.2,min(faults(f).Mchar+randn(1)*0.06,faults(f).Mchar+0.2));
               M(count,4)=1;
               M(count,5)=1/years;
               count=count+1;
            end      
        end   
    end
    
    % Sample earthquakes from background sources
    for i=1:length(backgrounds)
        ara(i)=backgrounds(i).avalue;
        arb(i)=backgrounds(i).bvalue;
        Mmina(i)=backgrounds(i).Mmin;
        Mmaxa(i)=backgrounds(i).Mmax;
        
        % Calculate min and max range for N, and P for each background source 
        Nmina(i)=10^(ara(i)-arb(i)*Mmina(i));
        Nmaxa(i)=10^(ara(i)-arb(i)*Mmaxa(i));
        Pmina(i)=1-exp(-Nmina(i));
        Pmaxa(i)=1-exp(-Nmaxa(i));

        randquake2=rand(1);
        if (randquake2>Pmaxa(i))&&(randquake2<Pmina(i))     % If true, then eq happens
            Nevent=-log(1-randquake2);
            M(count,1)=yr;
            M(count,2)=i;
            M(count,3)=(ara(i)-log10(Nevent))/arb(i);
            M(count,4)=0;
            M(count,5)=1/years;
            count=count+1;
        end
    end  
end