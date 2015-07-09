function [HCD,RGM]=HazAnaResult1NZ(Event,returnperiods,trueRGM)  
% This function calculates hazard curves for a set of earthquake scenarios,
% It is for the case when comparing results that are earthquake
% scenarios (for ground motion maps, use HazAnaResult2NZ.m)
% It is assumed that the input Event was generated using a run of SeisEvent1 
%  with number of eqs same as number of eqs in M (i.e., lnmean
%  and sigmas have already been determined) 

% Unlike in HazAnaResult2NZ.m, this function integrates over the ground
%   motion distribution to compute the exceedence probability
% It assumes ground motion has a normal distribution truncated at +/- 3 s.d.

% Only have to input trueRGM in case the reduced set hazard curve is flat at one
% of the input return periods. In that case, the function sets the ground motion 
% for that return period to be closest to trueRGM value.

% HCD is a matrix containing hazard curve info for each site 
% RGM is matrix containing ground motion intensity value at each return period

sitenum=length(Event(1,1).GM(:,1));     % Number of sites
count=length(Event);                    % Number of earthquake scenarios
% Define ground motion intensities at which to compute hazard curve
gm=0.001:0.001:3;             

HCD=zeros(sitenum,length(gm));

for i=1:sitenum
    for j=1:count   % for each eq j
        lnmean(i,j)=Event(1,j).GM(i,2); 
        residual(i,j)=sqrt((Event(1,j).GM(i,3))^2+(Event(1,j).GM(i,4))^2);
    end
end

for i=1:sitenum
    for j=1:count   % for each eq j
        % Define truncated ground motion intensity distribution
        pd=makedist('Normal','mu',lnmean(i,j),'sigma',residual(i,j));
        t=truncate(pd,lnmean(i,j)-3*residual(i,j),lnmean(i,j)+3*residual(i,j));
        for k=1:length(gm)
            site(i,k).eq(j)=cdf(t,log(gm(k)),'upper');
            wsite(i,k).eq(j)=site(i,k).eq(j)*Event(j).AP;
        end
%        disp(['site',num2str(i),'eq',num2str(j)]);
    end
    for k=1:length(gm)
        HCD(i,k)=sum(wsite(i,k).eq(:));
    end
end

g=1./returnperiods;     % Annual exceedence probabilities = 1/return periods

for j=1:length(g)
    for i=1:sitenum
        a=find(HCD(i,:)>g(j));      % Points on hazard curve above return period g(j)
        b=find(HCD(i,:)<g(j));      % Points on hazard curve below return period g(j)
        c=find(HCD(i,:)>0.9999999*g(j)&HCD(i,:)<1.0000001*g(j));    % Points on hazard curve equal to return period g(j)
        
        if numel(b)==0      % If no points below return period, use rightmost point on hazard curve 
           RGM(i,j)=gm(a(length(a)));
        elseif numel(c)>0 && length(trueRGM)>1    % If hazard curve is flat at return period g(j), use point closest to trueRGM
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
        else  % In most cases, interpolate hazard curve to obtain ground motion at return period g(j)
            RGM(i,j)=gm(a(length(a)))+(gm(b(1))-gm(a(length(a))))*(HCD(i,a(length(a)))-g(j))/(HCD(i,a(length(a)))-HCD(i,b(1)));
        end
%        disp(['sitenum',num2str(i)]);
    end
end

end