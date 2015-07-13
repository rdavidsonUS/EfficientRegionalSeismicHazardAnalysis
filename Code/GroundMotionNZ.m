function [lnGM,sigmaintra,tauinter,GM]=GroundMotionNZ(inputfolder,Tsec,Mag,ftype,Rrup,dbottom,dtop,soiltype)

% This function estimates the ground motion intensity in g at a single 
% specified site caused by a single specified earthquake scenario. 

% Outputs:
%   lnGM=ln(mean) in g
%   sigmaintra=intra-event residual
%   tauinter=inter-event residual
%   GM=mean ground motion intensity in g

% Indicate which McVerry coefficients to use by specifying the filename in line 59

% Location of coefficients in (Intensity parameter for each row in McVerry table)
% 1:PGA
% 2:PGA'
% 3:0.075s
% 4:0.10s
% 5:0.20s
% 6:0.30s
% 7:0.40s
% 8:0.50s
% 9:0.75s
%10:1.00s
%11:1.50s
%12:2.00s
%13:3.00s
if Tsec ==0
    refT =1;
elseif Tsec==0.075
    refT=3;
elseif Tsec==0.10
    refT=4;
elseif Tsec==0.2
    refT=5;
elseif Tsec==0.3
    refT=6;
elseif Tsec==0.4
    refT=7;
elseif Tsec==0.5
    refT=8;
elseif Tsec==0.75
    refT=9;
elseif Tsec==1
    refT=10;
elseif Tsec==1.5
    refT=11;
elseif Tsec==2
    refT=12;
elseif Tsec==3
    refT=13;
else
    refT=1;
    disp('Incorrect Tsec, this code will run  for PGA');
end

% Import coefficients in McVerry 2006 Table 5 or 6 
% (Change file name to switch table used)
fileID=fopen(sprintf('McVerryT5.csv',inputfolder));
tempcoef=textscan(fileID,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','Delimiter',',','Headerlines',1);
fclose(fileID);
fields={'C1','C3','C4','C5','C6','C8','C10','C11','C12','C13','C15','C17','C18','C19','C20','C24','C29','C30','C32','C33','C43','C46','SigmaM6','Sigslope','Tau','SigtotM6'};
% tempcoef=csvread(sprintf('McVerryT6.csv',inputfolder),1,0,[1,0,14,28]);
coef=cell2struct(tempcoef,fields,2);

%Define Rvol
if strcmp(ftype,'nv')==1
    Rvol = Rrup;
else
    Rvol=0;
end

%Define Hc Centroid depth
Hc=(dbottom+dtop)/2;

%Define CN,CR
if strcmp(ftype,'nv')==1 || strcmp(ftype,'nn')==1
    CN=-1; 
else
    CN=0;
end 
if strcmp(ftype,'rv')==1 
    CR = 1.0;
elseif strcmp(ftype,'sr')==1 || strcmp(ftype,'rs')
    CR=0.5;
else
    CR=0;
end

%Define SI, DS, deltaC, deltaD
if strcmp(ftype,'if')==1
    SI = 1;
else
    SI=0;
end 
DS=0;
if strcmp(ftype,'if')==0
    if Hc>=30
        DS=1; 
    end
end 
if strcmp(soiltype,'C')==1
    deltaC=1;
else
    deltaC=0;
end 
if strcmp(soiltype,'D')==1||strcmp(soiltype,'E')==1
    deltaD=1;
else
    deltaD=0;
end 

% Identify model
subduction =0;
if strcmp(ftype,'if')==1
    subduction =1;
elseif Hc>=30
    subduction =1;
end

if  subduction == 1
    % Subduction
    % Equation 2 to get ln(PGA_AB)
    lnPGAAB=coef.C11(1)+(coef.C12(1)+(coef.C15(1)-coef.C17(1))*coef.C19(1))*(Mag-6)...
        +coef.C13(1)*(10-Mag)^3+coef.C17(1)*log(Rrup+coef.C18(1)*exp(coef.C19(1)*Mag))...
        +coef.C20(1)*Hc+coef.C24(1)*SI+coef.C46(1)*Rvol*(1-DS);
    
    % Equation 2 to get ln(PGAprime_AB)
    lnPGAprimeAB=coef.C11(2)+(coef.C12(2)+(coef.C15(2)-coef.C17(2))*coef.C19(2))*(Mag-6)...
        +coef.C13(2)*(10-Mag)^3+coef.C17(2)*log(Rrup+coef.C18(2)*exp(coef.C19(2)*Mag))...
        +coef.C20(2)*Hc+coef.C24(2)*SI+coef.C46(2)*Rvol*(1-DS);
    
    if strcmp(soiltype,'D')==1 || strcmp(soiltype,'C')==1 || strcmp(soiltype,'E')
        if Tsec ==0
            % Equation 4 to get PGA_CD 
            % Using correction for Eqn. 4 where PGAprimeAM is replaced by SAprimeAB
            GM=exp(lnPGAAB+coef.C29(1)*deltaC+(coef.C30(1)*(log(exp(lnPGAAB)+0.03))+coef.C43(1))*deltaD);
        else
            % Equation 4 to get ln(PGAprime_CD) and ln(PGA_CD)
            lnPGAprimeCD=lnPGAprimeAB+coef.C29(2)*deltaC+(coef.C30(2)*(log(exp(lnPGAprimeAB)+0.03))+coef.C43(2))*deltaD;
            lnPGACD=lnPGAAB+coef.C29(1)*deltaC+(coef.C30(1)*(log(exp(lnPGAAB)+0.03))+coef.C43(1))*deltaD;
           
            % Equation 2 to get ln(SAprime_AB)
            lnSAprimeAB=coef.C11(refT)+(coef.C12(refT)+(coef.C15(refT)-coef.C17(refT))*coef.C19(refT))*(Mag-6)...
                +coef.C13(refT)*(10-Mag)^3+coef.C17(refT)*log(Rrup+coef.C18(refT)*exp(coef.C19(refT)*Mag))...
                +coef.C20(refT)*Hc+coef.C24(refT)*SI+coef.C46(refT)*Rvol*(1-DS);
            
            % Equation 4 to get ln(SAprime_CD)
            lnSAprimeCD=lnSAprimeAB+coef.C29(refT)*deltaC+(coef.C30(refT)*(log(exp(lnSAprimeAB)+0.03))+coef.C43(refT))*deltaD;

            % Equation 6 to get SA_CD
            GM=exp(lnSAprimeCD)*exp(lnPGACD)/exp(lnPGAprimeCD);
        end
    else
        % Soil types A/B
        if Tsec ==0
            % Get PGA_AB
            GM=exp(lnPGAAB);
        else
            % Equation 2 to get ln(SAprime_AB)
            lnSAprimeAB=coef.C11(refT)+(coef.C12(refT)+(coef.C15(refT)-coef.C17(refT))*coef.C19(refT))*(Mag-6)...
                +coef.C13(refT)*(10-Mag)^3+coef.C17(refT)*log(Rrup+coef.C18(refT)*exp(coef.C19(refT)*Mag))...
                +coef.C20(refT)*Hc+coef.C24(refT)*SI+coef.C46(refT)*Rvol*(1-DS);
            % Equation 6 to get SA_AB
            GM=exp(lnSAprimeAB)*exp(lnPGAAB)/exp(lnPGAprimeAB);
        end
    end
    % End of subduction
else
    % Crustal
    % No FHW term is included, as noted in McVerry et al. (2006, Section 6.1)
    % Equation 1 to get ln(PGAprime_AB) and ln(PGA_AB) 
    lnPGAAB=coef.C1(1)+coef.C4(1)*(Mag-6)+coef.C3(1)*((8.5-Mag)^2)+coef.C5(1)*Rrup...
        +(coef.C8(1)+coef.C6(1)*(Mag-6))*log(sqrt(Rrup^2+coef.C10(1)^2))+coef.C46(1)*Rvol...
        +coef.C32(1)*CN+coef.C33(1)*CR;
    
    lnPGAprimeAB=coef.C1(2)+coef.C4(2)*(Mag-6)+coef.C3(2)*((8.5-Mag)^2)+coef.C5(2)*Rrup...
        +(coef.C8(2)+coef.C6(2)*(Mag-6))*log(sqrt(Rrup^2+coef.C10(2)^2))+coef.C46(2)*Rvol...
        +coef.C32(2)*CN+coef.C33(2)*CR;
    
    if strcmp(soiltype,'D')==1 || strcmp(soiltype,'C')==1|| strcmp(soiltype,'E')==1
        if Tsec ==0
            % Equation 4 to get PGA_CD
            GM=exp(lnPGAAB+coef.C29(1)*deltaC+(coef.C30(1)*(log(exp(lnPGAAB)+0.03))+coef.C43(1))*deltaD);
        else
            % Equation 4 to get ln(PGAprime_CD) and ln(PGA_CD) 
            lnPGAprimeCD=lnPGAprimeAB+coef.C29(2)*deltaC+(coef.C30(2)*(log(exp(lnPGAprimeAB)+0.03))+coef.C43(2))*deltaD;
            lnPGACD=lnPGAAB+coef.C29(1)*deltaC+(coef.C30(1)*(log(exp(lnPGAAB)+0.03))+coef.C43(1))*deltaD;
           
            % Equation 1 to get ln(SAprime_AB)
            lnSAprimeAB=coef.C1(refT)+coef.C4(refT)*(Mag-6)+coef.C3(refT)*((8.5-Mag)^2)+coef.C5(refT)*Rrup...
                +(coef.C8(refT)+coef.C6(refT)*(Mag-6))*log(sqrt(Rrup^2+coef.C10(refT)^2))+coef.C46(refT)*Rvol...
                +coef.C32(refT)*CN+coef.C33(refT)*CR;

            % Equation 4 to get ln(SAprime_CD)
            lnSAprimeCD=lnSAprimeAB+coef.C29(refT)*deltaC+(coef.C30(refT)*(log(exp(lnSAprimeAB)+0.03))+coef.C43(refT))*deltaD;
            
            % Equation 6 to get SA_CD
            GM=exp(lnSAprimeCD)*exp(lnPGACD)/exp(lnPGAprimeCD);
        end
    else
        %Soil types A/B
        if Tsec ==0
            % Get PGA_AB
            GM=exp(lnPGAAB);
        else
            % Equation 1 to get ln(SAprime_AB)
            lnSAprimeAB=coef.C1(refT)+coef.C4(refT)*(Mag-6)+coef.C3(refT)*((8.5-Mag)^2)+coef.C5(refT)*Rrup...
                +(coef.C8(refT)+coef.C6(refT)*(Mag-6))*log(sqrt(Rrup^2+coef.C10(refT)^2))+coef.C46(refT)*Rvol...
                +coef.C32(refT)*CN+coef.C33(refT)*CR;
            % Equation 6 to get SA_AB
            GM=exp(lnSAprimeAB)*exp(lnPGAAB)/exp(lnPGAprimeAB);
        end
    end
    %End of crustal model
end

%Equations 8 and 9 in McVerry
tauinter=coef.Tau(refT);
if Mag<5
    sigmaintra=coef.SigmaM6(refT)-coef.Sigslope(refT);
elseif Mag<7
    sigmaintra=coef.SigmaM6(refT)+coef.Sigslope(refT)*(Mag-6);
else
    sigmaintra=coef.SigmaM6(refT)+coef.Sigslope(refT);
end
lnGM=log(GM);

end

