function [MHCEs]=sensNjrTONmcNZ(Neqsets,Mlargest,Plargest,TrueRGM,returnperiods,inputfolder,faults,backgrounds,sites,CorrCoef,Tsec)

% This function can be used to generate candidate ground motion maps for
% multiple sets of input earthquake scenarios (Neqsets) and multiple numbers of
% candidate ground motion maps (Nmc).

% It replaces Steps 3b (UpdateMNZ.m) and 4 (SeisEvent3NZ.m) when doing the
% analysis for just one reduced set of earthquake scenarios and one value
% of Nmc.

% The function sets a random seed so common random variables are used
% to reduce sampling variability in the comparison across runs 

% MHCEs is a matrix that stores the MHCE values for the different runs
% sensNjrTONmcOutput%d.mat is saved for each run and contains the resulting
%   RedM2,RedP2,and Event2 for each run
 
% Nmc is vector with the numbers of candidate ground motion maps to
%   generate (e.g., Nmc=[100 500 1000 2000 3000 4000 5000];)
% Mlargest and Plargers are the M and P structures for the largest set of 
% candidate earthquake scenarios considered (from EQReduceNZ.m)
% TrueRGM is true hazard curve info to compare to
% returnperiods is vector of returnperiods e.g. (250 500 1000 2500)
% CorrCoef=1 means include spatial correlation among the residuals
%         =0 means do NOT include spatial correlation among the residuals
% Tsec = Period in seconds of ground motion intensity metric (Tsec=0 refers to PGA) 

% Read in Pj's indicating which EQs were chosen in the 1st implementation of the optimization 
%  and what the new Pj's are 

Njc=zeros(Neqsets);
Njr=zeros(Neqsets);

for i=1:Neqsets                                 % Step through columns of Pjs.txt
    Njc(i)=dlmread('Pjs.txt','',[0 i 0 i]);   % Num candidate eq scenarios in opt 1 run
    Njr(i)=dlmread('Pjs.txt','',[1 i 1 i]);   % Num eq scenarios allowed in reduced set in opt 1 run (actual num may be smaller)
    Pjs(:,i)=dlmread('Pjs.txt','',[2 i (Njc(i)+2) i];    % Pj values output from each opt 1 run
end

Nruns=length(Njc)*length(Nmc);          % Num runs in this analysis to generate candidate GM maps

MHCEs=zeros(Nruns,5);                   % Row for each run
                                        % columns with Njc, Njr, Nmc, num maps/eq, and MHCE
s=rng;
count=1;
for j=1:length(Njc)     % Loop through runs from optimization 1 output
    % Make new smaller M and P structures to store eq scenario info and P(yij>=Yir)'s
    % Similar to what UpdateMNZ.m does for one run
    RedM2=Mlargest(Pjs(:,j)>0,:);
    RedM2(:,5)=Pjs(Pjs(:,j)>0,j);
    RedP2=Plargest(Pjs(:,j)>0);
    Njr_actual(j)=length(RedM2);        % Num eq scenarios in reduced set
    for k=1:length(Nmc)          % Loop through Nmc values
        if Njr_actual(j)<= Nmc(k)
            rng(s);              % Set random seed so common random variables are used across runs
            [Event2]=SeisEvent3NZ(inputfolder,RedP2,RedM2,faults,backgrounds,sites,CorrCoef,Tsec,TrueRGM,Nmc(k));
            [~,~,Means,~]=Summary2NZ('Given',returnperiods,Event2,999,TrueRGM);
            MHCEs(count,1)=Nmc(k);
            MHCEs(count,2)=Njc(j);
            MHCEs(count,3)=Njr(j);
            MHCEs(count,4)=Njr_actual(j);
            MHCEs(count,5)=Nmc(k)/Njr_actual(j);
            MHCEs(count,6)=Means(1,1);
            filenameout=sprintf('C:\Users\Davidson\Documents\New Zealand\GNS hazard project\From Natalia\sensNjrTONmcOutput%d.mat',count);
            save(filenameout,'MHCEs','RedM2','RedP2','Event2');
            count=count+1;    
        end
    end
end

figure;   % MHCE vs Num maps per eq scenario
scatter(MHCEs(:,5),MHCEs(:,6));
title('MHCE vs number of maps per earthquake scenario');
xlabel('Number of maps per earthquake scenario');
ylabel('Mean hazard curve error, MHCE');

end


