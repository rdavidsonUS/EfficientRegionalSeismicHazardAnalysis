function [Event]=SeisEventNZ(inputfolder,M,faults,backgrounds,sites,Tsec)

% This function computes the ground motion intensity value at every site i 
% for each specified earthquake scenario j. It does not consider spatial 
% correlation of ground motion residuals, and it calculates only one ground 
% motion intensity value for each earthquake scenario and each site. 

sitesnumber=length(sites);   % Number of sites

quakesize=size(M);
quakelength=quakesize(1);    % Number of earthquake scenarios
Event(1,quakelength)=struct('year',[],'SourceID',[],'rectype',[],'M',[],'GM',[]);

% Determine ground motions
for eq=1:quakelength
    Event(eq).year=single(M(eq,1));        % Year eq occurred if generated using conventional MCS or scenario ID if generated using importance sampling
    Event(eq).SourceID=single(M(eq,2));    % Source ID 
    Event(eq).rectype=single(M(eq,4));     % Indicator of background (0) or fault (1) source
    Event(eq).M=single(M(eq,3));           % Magnitude of eq
    for mo=1:sitesnumber
        if M(eq,4)==0
            [lnmean,vtra,vter,GM]=GroundMotionNZ(inputfolder,Tsec,M(eq,3),backgrounds(M(eq,2)).sensemov,backgrounds(M(eq,2)).Rrup(mo),backgrounds(M(eq,2)).depth,backgrounds(M(eq,2)).depth,sites(mo).soilt);
        else
            [lnmean,vtra,vter,GM]=GroundMotionNZ(inputfolder,Tsec,M(eq,3),faults(M(eq,2)).sensemov,faults(M(eq,2)).Rrup(mo),faults(M(eq,2)).dbottom,faults(M(eq,2)).dtop,sites(mo).soilt);
        end   
        Event(eq).GM(mo,1)=single(GM);      % Mean of ground motion 
        Event(eq).GM(mo,2)=single(lnmean);  % ln(mean) of ground motion
        Event(eq).GM(mo,3)=single(vtra);    % intra-event residual
        Event(eq).GM(mo,4)=single(vter);    % inter-event residual
    end
       
%    disp(['earthquakenumber=',num2str(eq)])
end  

end
