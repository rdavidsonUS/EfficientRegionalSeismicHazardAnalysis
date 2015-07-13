function [ Rrup ] = pfaultdistanceNZ( siteLAT,siteLON,ID,faults )
% Computes the Rrup distance (km) from a fault to a point using the method
% described in Kaklamanos et al. 2011.
% It first computes RJB and the source-to-site azimuth, then Rrup.

R=6371;     % Radius of the earth (km)
%Sine formula to find width of fault plane
width =(faults(ID).dbottom-faults(ID).dtop)/sin(faults(ID).dip*pi/180);
if faults(ID).dip==90
    widthSURFproj=0;
else
    widthSURFproj=width*cos(faults(ID).dip*pi/180);
end
% Create arrays for coordinates defining corners of surface projection of fault, 
% on strike (STRIKE) and not on strike (NOSTRIKE)
cooSTRIKElat=zeros(2,faults(ID).numseg);      
cooSTRIKElon=zeros(2,faults(ID).numseg);
cooNOSTRIKElat=zeros(2,faults(ID).numseg);
cooNOSTRIKElon=zeros(2,faults(ID).numseg);

cooSTRIKElat(1,1)=faults(ID).coostartLAT;
cooSTRIKElon(1,1)=faults(ID).coostartLON;
cooSTRIKElat(2,faults(ID).numseg)=faults(ID).cooendLAT;
cooSTRIKElon(2,faults(ID).numseg)=faults(ID).cooendLON;
if faults(ID).numseg>1
    for i = 2:faults(ID).numseg
        cooSTRIKElat(2,i-1)=faults(ID).segments(1,i-1);
        cooSTRIKElon(2,i-1)=faults(ID).segments(2,i-1); 
        cooSTRIKElat(1,i)=faults(ID).segments(1,i-1);
        cooSTRIKElon(1,i)=faults(ID).segments(2,i-1);
    end
end

% Find coordinates of corners of surface projection of fault that are not
% on strike
for i = 1:faults(ID).numseg
    % Check if fault is vertical (dip=90 deg) and if so, compute dist
    % (=Rrup for segment i)
    
    if faults(ID).dip==90

        % Distance between end points of segment
        faultlineD=ppdistanceNZ(cooSTRIKElat(1,i),cooSTRIKElon(1,i),0,cooSTRIKElat(2,i),cooSTRIKElon(2,i),0);
        % Distance between start point of segment and site
        startTOsiteD=ppdistanceNZ(cooSTRIKElat(1,i),cooSTRIKElon(1,i),0,siteLAT,siteLON,0);
        % Distance between end point of segment and site
        endTOsiteD=ppdistanceNZ(cooSTRIKElat(2,i),cooSTRIKElon(2,i),0,siteLAT,siteLON,0);
        % Cosine of angle opposite line connecting end point and site, using law of
        %   cosines in rectangular coordinates
        COSendTOsiteOpA=(faultlineD^2+startTOsiteD^2-endTOsiteD^2)/(2*faultlineD*startTOsiteD);
        % Vector of distance from start point to location of shortest
        %   distance from line of strike to site. The direction is given by
        %   minus or plus sign
        startTOhD=startTOsiteD*COSendTOsiteOpA;

        if startTOhD<0      % If site is before the start point
           RJB=startTOsiteD;
        elseif startTOhD>faultlineD     % If site is after the end point  
           RJB=endTOsiteD;
        else                % If site is between start and end points
           RJB=sqrt(startTOhD^2+startTOsiteD^2);
        end

        % Compute Rrup for segment i
        dist=sqrt(RJB^2+(faults(ID).dtop)^2);
    else
        % angle between lines from center of earth to strike and opposite side of surface projection
        b=widthSURFproj/R;
        % Get bearing angle of strike
        theta=atan2(sin((cooSTRIKElon(2,i)-cooSTRIKElon(1,i))*pi/180)*cos(cooSTRIKElat(2,i)*pi/180),cos(cooSTRIKElat(1,i)*pi/180)*sin(cooSTRIKElat(2,i)*pi/180)-sin(cooSTRIKElat(1,i)*pi/180)*cos(cooSTRIKElat(2,i)*pi/180)*cos((cooSTRIKElon(2,i)-cooSTRIKElon(1,i))*pi/180));
        % atan2 returns values between -pi and pi convert to 0 to 2pi
        if theta<0
            theta = 2*pi+theta;
        end
        % To get the bearing angle of the perpendicular plane +pi/2
        theta2=theta+pi/2;

        % If angle is greater than 2 pi reduce to first circle
        if theta2>2*pi
            theta2=theta2-2*pi;
        end

        %Using the bearing angle and distance get new coordinates
        cooNOSTRIKElat(1,i)=asin(sin(cooSTRIKElat(1,i)*pi/180)*cos(b)+cos(cooSTRIKElat(1,i)*pi/180)*sin(b)*cos(theta2))*180/pi;
        cooNOSTRIKElon(1,i)=cooSTRIKElon(1,i)+atan2(sin(theta2)*sin(b)*cos(cooSTRIKElat(1,i)*pi/180),cos(b)-sin(cooSTRIKElat(1,i)*pi/180)*sin(cooNOSTRIKElat(1,i))*pi/180)*180/pi;
        cooNOSTRIKElat(2,i)=asin(sin(cooSTRIKElat(2,i)*pi/180)*cos(b)+cos(cooSTRIKElat(2,i)*pi/180)*sin(b)*cos(theta2))*180/pi;
        cooNOSTRIKElon(2,i)=cooSTRIKElon(2,i)+atan2(sin(theta2)*sin(b)*cos(cooSTRIKElat(2,i)*pi/180),cos(b)-sin(cooSTRIKElat(2,i)*pi/180)*sin(cooNOSTRIKElat(2,i)*pi/180))*180/pi;

        % Define coordinates of 4 corners of surface projection of fault
        % (include 1st one again as 5th point to check 4th side of projection in loop)
        cornersLat=[cooSTRIKElat(1,i) cooSTRIKElat(2,i) cooNOSTRIKElat(1,i) cooNOSTRIKElat(2,i) cooSTRIKElat(1,i)];
        cornersLon=[cooSTRIKElon(1,i) cooSTRIKElon(2,i) cooNOSTRIKElon(1,i) cooNOSTRIKElon(2,i) cooSTRIKElon(1,i)];

        % First check if RJB=0
        if min(cornersLat)<=siteLAT && siteLAT<=max(cornersLat) && min(cornersLon)<=siteLON && siteLON<=max(cornersLon)
            RJB=0;
        else  % If RJB not zero, continue
            for k=1:4
                % Distance between end points of segment
                faultlineD=ppdistanceNZ(cornersLat(k),cornersLon(k),0,cornersLat(k+1),cornersLon(k+1),0);
                % Distance between start point of segment and site
                startTOsiteD=ppdistanceNZ(cornersLat(k),cornersLon(k),0,siteLAT,siteLON,0);
                % Distance between end point of segment and site
                endTOsiteD=ppdistanceNZ(cornersLat(k+1),cornersLon(k+1),0,siteLAT,siteLON,0);
                % Cosine of angle opposite line connecting end point and site, using law of
                %   cosines in rectangular coordinates
                COSendTOsiteOpA=(faultlineD^2+startTOsiteD^2-endTOsiteD^2)/(2*faultlineD*startTOsiteD);
                % Vector of distance from start point to location of shortest
                %   distance from line of strike to site. The direction is given by
                %   minus or plus sign
                startTOhD=startTOsiteD*COSendTOsiteOpA;

                if startTOhD<0      % If site is before the start point
                   RJBtemp=startTOsiteD;
                elseif startTOhD>faultlineD     % If site is after the end point  
                   RJBtemp=endTOsiteD;
                else                % If site is between start and end points
                   RJBtemp=sqrt(startTOhD^2+startTOsiteD^2);
                end
                if k==1 || RJBtemp<RJB
                    RJB=RJBtemp;
                end
            end
        end
        
        % Find source-to-site azimuth = angle between positive fault strike
        % direction and line connecting site to closest point on surface
        % projection of top edge of rupture, with clockwise angles assumed
        % positive (Kaklamanos et al. 2001, p.1224)

        % Distance between start and end point of segment
        faultlineD=ppdistanceNZ(cooSTRIKElat(1,i),cooSTRIKElon(1,i),faults(ID).dtop,cooSTRIKElat(2,i),cooSTRIKElon(2,i),faults(ID).dtop);
        % Distance between start point of segment and site
        startTOsiteD=ppdistanceNZ(cooSTRIKElat(1,i),cooSTRIKElon(1,i),faults(ID).dtop,siteLAT,siteLON,0);
        % Distance between end point of segment and site
        endTOsiteD=ppdistanceNZ(cooSTRIKElat(2,i),cooSTRIKElon(2,i),faults(ID).dtop,siteLAT,siteLON,0);
        % Cosine of angle opposite line connecting end point and site, using law of
        %   cosines in rectangular coordinates
        COSendTOsiteOpA=(faultlineD^2+startTOsiteD^2-endTOsiteD^2)/(2*faultlineD*startTOsiteD);
        % Vector of distance from start point to location of shortest
        %   distance from line of strike to site. The direction is given by
        %   minus or plus sign
        startTOhD=startTOsiteD*COSendTOsiteOpA;

        %Get absolute value of azimuth of site versus strike line. 
        if startTOhD<0      % If site is before the start point
           ABSazimuthF=pi-acos(startTOhD/startTOsiteD);
        elseif startTOhD>faultlineD     %If site is after the end point  
           ABSazimuthF=acos(startTOhD/startTOsiteD);
        else
           ABSazimuthF=pi/2;
        end

        % Identify location with respect to line through strike using same
        % logic but now side of surface projection is considered strike
        faultlineD2=ppdistanceNZ(cooSTRIKElat(2,i),cooSTRIKElon(2,i),0,cooNOSTRIKElat(2,i),cooNOSTRIKElon(2,i),0);
        startTOsiteD2=ppdistanceNZ(cooSTRIKElat(2,i),cooSTRIKElon(2,i),0,siteLAT,siteLON,0);
        endTOsiteD2=ppdistanceNZ(cooNOSTRIKElat(2,i),cooNOSTRIKElon(2,i),0,siteLAT,siteLON,0);
        COSendTOsiteOpA2=(faultlineD2^2+startTOsiteD2^2-endTOsiteD2^2)/(2*faultlineD2*startTOsiteD2);
        startTOhD2=startTOsiteD2*COSendTOsiteOpA2;
        if startTOhD2<0    % If site is before the start point - left of strike
           azimuthF=-ABSazimuthF;
        else
           azimuthF=ABSazimuthF;
        end

        % Compute Rx using equations in Table 3 in Kaklamanos et al. 2011
        if (0<=azimuthF && azimuthF<(pi/2)) || ((pi/2)<azimuthF && azimuthF<=pi)
            if RJB*abs(tan(azimuthF)) <= widthSURFproj
                Rx=RJB*abs(tan(azimuthF));
            else
                Rx=RJB*tan(azimuthF)*cos(azimuthF - asin(widthSURFproj*cos(azimuthF)/RJB));
            end
        elseif azimuthF==(pi/2)
            if RJB>0    
                Rx=RJB+widthSURFproj;
            else
                Rx=startTOsiteD*sin(COSendTOsiteOpA);
            end    
        else
            Rx=RJB*sin(azimuthF);
        end

        % Compute primeRrup (Kaklamanos et al. 2011 Eqns. 15-17)
        if Rx<faults(ID).dtop*tan(pi*faults(ID).dip/180)
            primeRrup=sqrt(Rx^2+faults(ID).dtop^2);
        elseif Rx<=faults(ID).dtop*tan(pi*faults(ID).dip/180)+width*sec(pi*faults(ID).dip/180)
            primeRrup=Rx*sin(pi*faults(ID).dip/180)+faults(ID).dtop*cos(pi*faults(ID).dip/180);
        else
            primeRrup=sqrt((Rx-width*cos(pi*faults(ID).dip/180))^2+(faults(ID).dtop+width*sin(pi*faults(ID).dip/180))^2);
        end

        % Compute Ry (Kaklamanos et al. 2011 Eqns. 18-20)
        if  azimuthF==pi/2||azimuthF==-pi/2
            Ry=0;
        elseif azimuthF==0||azimuthF==-pi||azimuthF==pi
            Ry=RJB;
        else
            Ry=abs(Rx*cot(azimuthF));
        end
        % Compute Rrup for segment i (Kaklamanos et al. 2011 Eqn. 14)
        dist=sqrt(primeRrup^2+Ry^2);
     end
   
     % Update Rrup for fault
     if i==1||Rrup>dist
        Rrup=dist;
     end
end
end

