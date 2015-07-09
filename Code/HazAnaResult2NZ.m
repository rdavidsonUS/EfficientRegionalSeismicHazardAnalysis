function [HCD,RGM]=HazAnaResult2NZ(Event,returnperiods,trueRGM)    

% This function calculates hazard curves for a set of ground motion maps,
% It is for the case when comparing results that are ground motion maps (for earthquake scenarios, use HazAnaResult1NZ.m)

% Unlike in HazAnaResult1NZ.m, which integrates over the ground
%   motion distribution to compute the exceedence probability, this
%   function does not.

% Only have to input trueRGM in case the reduced set hazard curve is flat at one
% of the input return periods. In that case, the function sets the ground motion 
% for that return period to be closest to trueRGM value.

% HCD is a matrix containing hazard curve info for each site 
% RGM is matrix containing ground motion intensity value at each return period

sitenum=length(Event(1,1).GM(:,1));     % Number of sites
count=length(Event);                    % Number of ground motion maps
% Define ground motion intensities at which to compute hazard curve
k=0.001:0.001:3;
HCD=zeros(sitenum,length(k));

for i=1:sitenum
    for j=1:count
        site(i,j)=Event(1,j).GM(i,1);
    end
end

for i=1:sitenum
    n=1;
    for gm=0.001:0.001:3
        a=find(site(i,:)>gm);
        for j=1:length(a)
            HCD(i,n)=HCD(i,n)+Event(a(j)).AP;
        end
        n=n+1;
    end
%    disp(['sitenum',num2str(i)]);
end

g=1./returnperiods;        % Annual exceedence probabilities = 1/return periods

for j=1:length(g)
    for i=1:sitenum
        a=find(HCD(i,:)>g(j));  % Points on hazard curve above return period g(j)
        b=find(HCD(i,:)<g(j));  % Points on hazard curve below return period g(j)
        c=find(HCD(i,:)>0.9999999*g(j)&HCD(i,:)<1.0000001*g(j));    % Points on hazard curve equal to return period g(j)
        
        if numel(b)==0      % If no points below return period, use rightmost point on hazard curve
           RGM(i,j)=gm(a(length(a)));
        elseif numel(c)>0 && length(trueRGM)>1          % If hazard curve is flat at return period g(j), use point closest to trueRGM
            % If doing as part of 'MCS' case, trueRGM=999 and you can't do this
            if k(c(1))<=trueRGM(i,j) 
                if k(c(length(c)))>=trueRGM(i,j)
                    RGM(i,j)=trueRGM(i,j);
                else
                    RGM(i,j)=k(c(length(c)));
                end
            else
                RGM(i,j)=k(c(1));
            end
        else   % In most cases, interpolate hazard curve to obtain ground motion at return period g(j)
            RGM(i,j)=k(a(length(a)))+(k(b(1))-k(a(length(a))))*(HCD(i,a(length(a)))-g(j))/(HCD(i,a(length(a)))-HCD(i,b(1)));
        end
%        disp(['sitenum',num2str(i)]);
    end
end

end