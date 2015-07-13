function [ faults,backgrounds,sites, trueRGM ] = loadingNZ( inputfolder,maxD )

% This function imports the input files and stores the data in variables. 
% To save space, fault sources, background seismicity sources, and “true” 
% hazard data (if included) are only saved if they are within the 
% user-specified distance maxD of the study area boundaries. This procedure 
% is designed assuming a rectangular study area defined by the site data. 
% For efficiency, the Rrup distances are saved for each fault and background 
% source for later use. It also assumes Mmin=5.25 for all background sources 
% (set in line 205).

% inputfolder: Location of all input files 
% maxD: Maximum distance of sources to be included (in km) 
%       (if maxD<0, all sources are included).

% Create study area structure
studyarea.north=0;
studyarea.south=0;
studyarea.west=0;
studyarea.east=0;

% Loads site information and study area coordinates
txt1=sprintf('chch_pilot_sites.csv',inputfolder);
Openfile1=fopen(txt1);
fileID4=textscan(Openfile1,'%f %f %s','Delimiter',',','Headerlines',1);
fclose(Openfile1);
fields={'lon','lat','soilt'};
tempSITE=cell2struct(fileID4,fields,2);

for i=1:length(tempSITE.lat(:))
    sites(i).lat=tempSITE.lat(i);
    if i==1||sites(i).lat<studyarea.south
        studyarea.south=sites(i).lat;
    end
    if i==1||sites(i).lat>studyarea.north
        studyarea.north=sites(i).lat;
    end
    sites(i).lon=tempSITE.lon(i);
    if i==1||sites(i).lon>studyarea.east
        studyarea.east=sites(i).lon;
    end
    if i==1||sites(i).lon<studyarea.west
        studyarea.west=sites(i).lon;
    end
    sites(i).soilt=tempSITE.soilt(i);
end

% Create matrix of edge coordinates
count=1;
for i=1:length(sites)
    if sites(i).lat==studyarea.south || sites(i).lat==studyarea.north
        edgecoo(count,1)=sites(i).lat;
        edgecoo(count,2)=sites(i).lon;
        count = count+1;
    elseif sites(i).lon==studyarea.east || sites(i).lon==studyarea.west
        edgecoo(count,1)=sites(i).lat;
        edgecoo(count,2)=sites(i).lon;
        count=count+1;
    end
end

% Import given "true" hazard, eliminating points different than in region
% of study
fileID5=sprintf('%sNSHM_PGA.xyz',inputfolder);
tempGM=dlmread(fileID5);

for i=1:length(tempGM)
    for j=1:length(sites)
        if abs(sites(j).lat-tempGM(i,2))<0.001 && abs(sites(j).lon-tempGM(i,1))<0.001
            % Collects ground motion only for sites and gives a value of zero
            % to points with value -NaN
            if isnan(tempGM(i,3))==0
                trueRGM(j,1)=tempGM(i,3);
            else
                trueRGM(j,1)=0;
            end
            if isnan(tempGM(i,4))==0
                trueRGM(j,2)=tempGM(i,4);
            else
                trueRGM(j,2)=0;
            end
        	if isnan(tempGM(i,5))==0
                trueRGM(j,3)=tempGM(i,5);
            else
                trueRGM(j,3)=0;
            end
            if isnan(tempGM(i,6))==0
                trueRGM(j,4)=tempGM(i,6);
            else
                trueRGM(j,4)=0;
            end
        end
    end
end


% Import fault information
fileID1=fopen(sprintf('FUN1111.DAT',inputfolder));

finfo = dir(sprintf('FUN1111.DAT',inputfolder));
fsize = finfo.bytes;

if fsize>0
    % Matrix to load the fault segments coordinates. Modify the size if
    % more segments than in given faults
    
    % Skip the first three lines
    tline = fgetl(fileID1);
    tline = fgetl(fileID1);
    tline = fgetl(fileID1);
    count=1;
    while ~feof(fileID1)
        clear cooF cooF2;
        %	Fault name, sense of movement (rv=reverse, sr=strike slip reverse, 
        %nv=normal volcanic, ss=pure strike slip, sn=strike slip normal)
        faults(count).name= fscanf(fileID1, '%s', 1);
        if ~isempty(faults(count).name)
        faults(count).sensemov= fscanf(fileID1, '%s', 1);
        % Number of segments
        faults(count).numseg= fscanf(fileID1, '%dD', 1);
        % dbottom = depth to bottom, dtop = depth to top
        faults(count).dip=fscanf(fileID1, '%f', 1);
        faults(count).dipdir=fscanf(fileID1, '%f', 1);
        faults(count).dbottom=fscanf(fileID1, '%f', 1);
        faults(count).dtop=fscanf(fileID1, '%f', 1);
        % Coordinates are in degrees and decimal minutes -> convert to
        % decimal-degrees
        degree = fscanf(fileID1, '%f', 1);
        minutes=fscanf(fileID1, '%f', 1);
        faults(count).coostartLAT=-(degree+minutes/60);
        degree = fscanf(fileID1, '%f', 1);
        minutes=fscanf(fileID1, '%f', 1);
        faults(count).coostartLON=degree+minutes/60;
        degree = fscanf(fileID1, '%f', 1);
        minutes=fscanf(fileID1, '%f', 1);
        faults(count).cooendLAT=-(degree+minutes/60);
        degree = fscanf(fileID1, '%f', 1);
        minutes=fscanf(fileID1, '%f', 1);
        faults(count).cooendLON=(degree+minutes/60);
        faults(count).Mchar=fscanf(fileID1, '%f', 1);
        faults(count).Tyear=fscanf(fileID1, '%f', 1);
        % Coordinates per segment in temporary matrix
        cooF=fscanf(fileID1, '%f', [faults(count).numseg,8]);
        cooF2=zeros(length(cooF(:,1)),length(cooF(1,:)));
        numcoo=1;
        for j=1:length(cooF(1,:))
            for i=1:length(cooF(:,1))
                cooF2(1+floor((numcoo-1)/8),numcoo-(floor((numcoo-1)/8))*8)=cooF(i,j);
                numcoo=numcoo+1;
            end
        end
        % Pass coordinates of intermediate points to structure. 
        for i=1:(faults(count).numseg-1)
            faults(count).segments(1,i)=-(cooF2(i,5)+cooF2(i,6)/60); % LAT
            faults(count).segments(2,i)=cooF2(i,7)+cooF2(i,8)/60; % LON
        end
        endID=fscanf(fileID1, '%d', 1);
        if endID ~=-1
            disp('Error importing faultdata');
        end
        if maxD >0
            shortestdistance = 0;
            for i =1:length(edgecoo)
                % Check distance between fault and set of points on the edge
                % of study area
                distance=pfaultdistanceNZ( edgecoo(i,1),edgecoo(i,2),count,faults );
                if i==1 ||shortestdistance>distance
                    shortestdistance=distance;
                end
            end
            if shortestdistance<=maxD
                % If shortestdistance is within desired maxD keep fault and 
                % compute Rrup for all sites, otherwise remove
                for k=1:length(sites)
                    faults(count).Rrup(k)=pfaultdistanceNZ( sites(k).lat,sites(k).lon,count,faults );
                end
                count = count +1;
            else
                % Clean faults
                faults(count)=[];
            end
        else
            count = count +1;
        end
        else
            faults(count)=[];
        end
    end
end
fclose(fileID1);

% Read background seismicity information
fileID2=fopen(sprintf('NZBCK211.DAT',inputfolder));

finfo = dir(sprintf('NZBCK211.DAT',inputfolder));
fsize = finfo.bytes;

if fsize>0
    % Matrix to load the fault segments coordinates. Modify the size if
    % more segments than in given faults
    tempGR=zeros(1,6);
    count=1;
    for i=1:6
        tempGR(1,i)=fscanf(fileID2, '%f', 1);
    end
    notuse=fscanf(fileID2, '%f', 1);
    while ~feof(fileID2)
        N1= fscanf(fileID2, '%f', 1);
        if ~isempty(N1)
        N2= fscanf(fileID2, '%f', 1);
        N3=fscanf(fileID2, '%f', 1);
        backgrounds(count).bvalue=fscanf(fileID2, '%f', 1);
        backgrounds(count).Mmin=5.25;
        backgrounds(count).Mmax=fscanf(fileID2, '%f', 1);
        notuse=fscanf(fileID2, '%f', 1);
        backgrounds(count).sensemov=fscanf(fileID2, '%s', 1);
        if  strcmp(backgrounds(count).sensemov,'nv')==1||strcmp(backgrounds(count).sensemov,'nn')==1||strcmp(backgrounds(count).sensemov,'ss')==1||strcmp(backgrounds(count).sensemov,'sr')==1||strcmp(backgrounds(count).sensemov,'rv')==1||strcmp(backgrounds(count).sensemov,'if')==1||strcmp(backgrounds(count).sensemov,'ro')==1
            notuse=fscanf(fileID2, '%f', 1);
        else
            backgrounds(count).sensemov='ss';
        end
        notuse=fscanf(fileID2, '%f', 1);
        backgrounds(count).cooLAT=-fscanf(fileID2, '%f', 1);
        backgrounds(count).cooLON=fscanf(fileID2, '%f', 1);
        backgrounds(count).depth=fscanf(fileID2, '%f', 1);

        if maxD >0
            % Identify location of point with respect to study area
            % Estimate distance according to location
            if backgrounds(count).cooLAT<studyarea.south 
            % Three bottom quadrants. 
                if backgrounds(count).cooLON< studyarea.west
                % Quadrant west so distance is to SW corner
                    distance=ppdistanceNZ(backgrounds(count).cooLAT,backgrounds(count).cooLON,backgrounds(count).depth,studyarea.south,studyarea.west,0);
               elseif backgrounds(count).cooLON< studyarea.east
                % It is the middle quadrant to the south. Assume closest
                % distance in same LON as point
                    distance=ppdistanceNZ(backgrounds(count).cooLAT,backgrounds(count).cooLON,backgrounds(count).depth,studyarea.south,backgrounds(count).cooLON,0);
                else
                % SE quadrant. Distance to corner point
                    distance=ppdistanceNZ(backgrounds(count).cooLAT,backgrounds(count).cooLON,backgrounds(count).depth,studyarea.south,studyarea.east,0);
                end
            elseif backgrounds(count).cooLAT<= studyarea.north
            % Between South and North of area
                if backgrounds(count).cooLON< studyarea.west
                    % We compare to West edge but at LAT of point
                    distance=ppdistanceNZ(backgrounds(count).cooLAT,backgrounds(count).cooLON,backgrounds(count).depth,backgrounds(count).cooLAT,studyarea.west,0);
                elseif backgrounds(count).cooLON<= studyarea.east
                % Right in study area
                    distance = 0;
                else
                    distance=ppdistanceNZ(backgrounds(count).cooLAT,backgrounds(count).cooLON,backgrounds(count).depth,backgrounds(count).cooLAT,studyarea.east,0);
                end
            else
                if backgrounds(count).cooLON< studyarea.west
                    distance=ppdistanceNZ(backgrounds(count).cooLAT,backgrounds(count).cooLON,backgrounds(count).depth,studyarea.north,studyarea.west,0);
                elseif backgrounds(count).cooLON< studyarea.east
                    distance=ppdistanceNZ(backgrounds(count).cooLAT,backgrounds(count).cooLON,backgrounds(count).depth,studyarea.north,backgrounds(count).cooLON,0);
                else
                    distance=ppdistanceNZ(backgrounds(count).cooLAT,backgrounds(count).cooLON,backgrounds(count).depth,studyarea.north,studyarea.east,0);
                end
            end
            if distance<=maxD
                % If shortestdistance is within desired maxD, keep background 
                % source and compute Rrup for all sites, else remove
                tempvalue=tempGR(1,2)*(10^(-tempGR(1,1)*backgrounds(count).bvalue))+tempGR(1,4)*(10^(-tempGR(1,3)*backgrounds(count).bvalue))+tempGR(1,6)*(10^(-tempGR(1,5)*backgrounds(count).bvalue));
                backgrounds(count).avalue=log10((N1+N2+N3)/(tempvalue));
                for k=1:length(sites)
                    backgrounds(count).Rrup(k)=ppdistanceNZ( sites(k).lat,sites(k).lon,0,backgrounds(count).cooLAT,backgrounds(count).cooLON,backgrounds(count).depth );
                end

                count = count +1;
            else
                % Clean source
                backgrounds(count)=[];
            end
        else
            tempvalue=tempGR(1,2)*(10^(-tempGR(1,1)*backgrounds(count).bvalue))+tempGR(1,4)*(10^(-tempGR(1,3)*backgrounds(count).bvalue))+tempGR(1,6)*(10^(-tempGR(1,5)*backgrounds(count).bvalue));
            backgrounds(count).avalue=log10((N1+N2+N3)/(tempvalue));
            count = count +1;
        end
        end
    end
end
fclose(fileID2);

end

