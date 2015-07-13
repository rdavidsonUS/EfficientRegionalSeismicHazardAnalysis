function [Event]=SeisEvent1NZ(inputfolder,M,faults,backgrounds,sites,CorrCoef,Tsec,num)

% This function is the same as SeisEventNZ.m but can include spatial 
% correlation among residuals (according to Jayaram and Baker 2009) 
% and simulates multiple ground motion maps for each earthquake scenario.
% Ground motion distribution is truncated at +/- 3 standard deviations

% CorrCoef=1 means include spatial correlation among the residuals
%         =0 means do NOT include spatial correlation among the residuals
% num = Total number of ground motion maps generated from all the eq scenarios
%       (They are distributed equally across earthquake scenarios). 

% Tsec = Period in seconds of ground motion intensity metric 
%        (Tsec=0 refers to PGA) 

b=40.7-15*Tsec;

randseed=1;
rand('state',randseed);
sitesnumber=length(sites);   % Number of sites
k=1;

if mod(num,length(M))~=0
    binnum1=floor(num/length(M))*ones(length(M),1);
    binnum2=zeros(length(M),1);
    y=randsample(length(M),mod(num,length(M)));
    for i=1:length(y)
        binnum2(y(i))=1;
    end
    binnum=binnum1+binnum2;
else
    binnum=num/length(M)*ones(length(M),1);
end

for i=1:length(M)
    for j=1:binnum(i)
        M(k,1)=k;                           % ID of EQ
        M(k,2)=M(i,2);                      % Generated EQ's source number
        M(k,3)=M(i,3)-0.05+(2*j-1)/(20*binnum(i));  % Generate EQ from selected bin chosenbin(i), the number of EQS is based on the value of binmark(cb(i)))
        M(k,4)=M(i,4);                      %  Indicator of background (0) or fault (1) source
        M(k,5)=M(i,5)/length(binnum(i));    %  Annual occurrence probability 
        k=k+1;
    end
end

eqnum=length(M);
[Event]=SeisEventNZ(inputfolder,M,faults,backgrounds,sites,Tsec);

% If NO spatial correlation among residuals
if CorrCoef == 0
    for i=1:eqnum
        rands=randn(1,2);   % Generate 2 random numbers, one for inter-event residual and one for intra-event residual
        rands(rands>3)=3;   % Truncate at -/+ 3 sd
        rands(rands<-3)=-3;
        for j=1:sitesnumber
            Event(i).GM(j,1)=exp(Event(i).GM(j,2)+rands(1,1)*Event(i).GM(j,3)+rands(1,2)*Event(i).GM(j,4));
        end
    end
end

% If spatial correlation among residuals
if CorrCoef == 1
    for i=1:sitesnumber
        for j=1:sitesnumber
            d(i,j)=ppdistanceNZ(sites(i).lat,sites(i).lon,0,sites(j).lat,sites(j).lon,0);
            R(i,j)=exp(-3*d(i,j)/b);  % Correlation coefficient from Jayaram and Baker (2009)
        end
    end
    
    MU=zeros(1,sitesnumber);
    for i=1:eqnum
        for j=1:sitesnumber
            c(j)=Event(i).GM(j,3);
        end    
        D=diag(c(j));           % Diagonal marix of intra-event standard deviation of each site for eq i
        S=D*R*D;                % Compute covariance matrix
        randintra=mvnrnd(MU,S); % Sample values from multivariate normal distribution
        for j=1:sitesnumber     % Truncate at -/+ 3 sd
            randintra(randintra(j)>(3*Event(i).GM(j,3)))=3*Event(i).GM(j,3);
            randintra(randintra(j)<(-3*Event(i).GM(j,3)))=-3*Event(i).GM(j,3);
        end 
        randinter=(randn(1)*Event(i).GM(1,4));
        randinter(randinter>(3*Event(i).GM(1,4)))=3*Event(i).GM(1,4);
        randinter(randinter<(-3*Event(i).GM(1,4)))=(-3*Event(i).GM(1,4));
        randinter=randinter*ones(1,sitesnumber);
        
        for j=1:sitesnumber
            Event(i).GM(j,1)=exp(Event(i).GM(j,2)+randinter(j)+randintra(j));
        end
    end
end

for i=1:eqnum
    Event(i).AP=M(i,5);         %  Annual occurrence probability
end
end