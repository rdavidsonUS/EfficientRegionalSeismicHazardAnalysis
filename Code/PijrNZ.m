function P=PijrNZ(M,sites,trueRGM,Event)

% This function computes the exceedance probability P(y_ij?Y_ir ) for each 
% earthquake scenario j, site i, and return period r, where Y_ir is the 
% ground shaking from the “true” hazard curve for return period r at site i, 
% and y_ij is the ground shaking at site i caused by earthquake scenario j.

% Event is result of running SeisEvent1NZ.m with num=number of eqs in M 
%  Just to get mean and sigma's

sitenum=length(sites);                  % Number of sites
eqnum=length(M);                        % Number of earthquake scenarios
lnmean=zeros(sitenum,eqnum,'double');
residual=zeros(sitenum,eqnum,'double');

for i=1:sitenum
    for j=1:eqnum
        lnmean(i,j)=Event(1,j).GM(i,2);    % Mean value of ground motion at site i caused by eq j
        residual(i,j)=sqrt((Event(1,j).GM(i,3))^2+(Event(1,j).GM(i,4))^2);  % standard divation of residual of ground motion at site i caused by eq j
    end
    disp(['sitenum',num2str(i)]);
end

for i=1:sitenum
    for j=1:eqnum
        for r=1:length(trueRGM(1,:))    
            % Integrate over ground motion distribution, truncating at -/+ 3sd
            pd=makedist('Normal','mu',lnmean(i,j),'sigma',residual(i,j));
            t=truncate(pd,lnmean(i,j)-3*residual(i,j),lnmean(i,j)+3*residual(i,j));
            P(j).gm(i,r)=cdf(t,log(trueRGM(i,r)),'upper');
        end
    end
    disp(['sitenum',num2str(i)]);
end

end