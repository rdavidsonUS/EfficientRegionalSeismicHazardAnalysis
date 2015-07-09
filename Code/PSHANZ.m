function [AEP,pshaRGM]=PSHANZ(inputfolder,backgrounds,faults,sites,Tsec,returnperiods)

% This function performs a probabilistic seismic hazard analysis (PSHA) to 
% compute the hazard curve at each specified site
% It assumes ground motion has a normal distribution truncated at +/- 3 s.d.

% AEP contains matrix of hazard curve information at each site
% pshaRGM contains ground motion intensity value at each return period 

% Tsec is period of specified ground motion intensity parameter (Tsec=0 means PGA)
% returnperiods is vector of returnperiods e.g. (250 500 1000 2500)

% Define values of ground motion at which to evaluate hazard curve
if Tsec==0  % If PGA, use y=0g to 2g at 0.05g intervals
    ybins=0.05:0.05:2;    
else        % If Sa, use y=0g to 5g at 0.05g intervals
    ybins=0.05:0.05:5;    
end

P=zeros(1,length(ybins));
AEP=zeros(length(ybins),(1+length(sites)));
pshaRGM=zeros(length(sites),length(returnperiods));

for i=1:length(sites)
    Lambday=zeros(1,length(ybins));
    for b=1:length(backgrounds)
        Mbin=backgrounds(b).Mmin:0.1:backgrounds(b).Mmax; % divide the magnitude range of i into partitions by 0.1
        for m=1:(length(Mbin)-1)
            GRa(b)=backgrounds(b).avalue;
            GRb(b)=backgrounds(b).bvalue;
            % Peq=Rate of eq from source b of magnitude m
            Peq=exp(-(10^(GRa(b)-GRb(b)*Mbin(m+1))))-exp(-(10^(GRa(b)-GRb(b)*Mbin(m))));
            [lnmean,vtra,vter,GM]=GroundMotionNZ(inputfolder,Tsec,(Mbin(m)+0.05),backgrounds(b).sensemov,backgrounds(b).Rrup(i),backgrounds(b).depth,backgrounds(b).depth,sites(i).soilt);
            residual=sqrt(vtra^2+vter^2);
            % Define truncated normal distribution
            pd=makedist('Normal','mu',lnmean,'sigma',residual);
            t=truncate(pd,lnmean-3*residual,lnmean+3*residual);
            for y=1:length(ybins)
                P(1,y)=cdf(t,log(ybins(y)),'upper');
            end
            % Annual rate at which y is exceeded =
            %   sum of annual rate due to eq from source b of mag m
            Lambday=Lambday+Peq*P;
        end
% disp(['Finished background source ',num2str(b)]);
    end
    for f=1:length(faults)
            Peq=1/faults(f).Tyear;
            [lnmean,vtra,vter,GM]=GroundMotionNZ(inputfolder,Tsec,faults(f).Mchar,faults(f).sensemov,faults(f).Rrup(i),faults(f).dbottom,faults(f).dtop,sites(i).soilt);
            residual=sqrt(vtra^2+vter^2);
            % Define truncated normal distribution
            pd=makedist('Normal','mu',lnmean,'sigma',residual);
            t=truncate(pd,lnmean-3*residual,lnmean+3*residual);
            for y=1:length(ybins)
                P(1,y)=cdf(t,log(ybins(y)),'upper');
            end
            % Annual rate at which y is exceeded =
            %   sum of annual rate due to eq from source b of mag m
            Lambday=Lambday+Peq*P;
    end
    % Compute annual exceedence probability for each site i and ground motion y
    AEP(:,1)=ybins';
    AEP(:,1+i)=(1-exp(-Lambday))';
% disp(['Finished site ',num2str(i)]);
end

% Add computation of RGM by linearly interpolating hazard curve 
xq=1./returnperiods;   % Define query points as 1/return periods
for i=1:length(sites)
    x=AEP(:,1+i);                % Define sample points x as AEPs 
    v=AEP(:,1);                   % Define values as ground motions y
    pshaRGM(i,:)=interp1(x',v',xq);  % Interpolated ground motions at specified return periods
end
end

