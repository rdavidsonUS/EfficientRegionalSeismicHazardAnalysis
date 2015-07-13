function [RedM,RedP]=EQreduceNZ(M,P,contribpct)

% This function takes a set of earthquake scenarios in M and selects those
% that together make up 'contribpct' of the contribution
% "contribpct" is from 0 to 1. Due to the averaging, sometimes it goes over
% one, so to include all eqs, let contribpct be 2

% It also makes figures of the cumulative contribution vs number of eqs 
%    (a) for all eqs and (b) for only those selected eqs

sitenum=length(P(1).gm(:,1));       % Number of sites i
numofr=length(P(1).gm(1,:));        % Number of return periods r
contribution=zeros(1,length(M));    % Number of eq scenarios j

for i=1:sitenum
    for r=1:numofr
        for j=1:length(P)
            % Numerator in brackets in Han and Davidson (2012, Eq. 1)
            a(j)=P(1,j).gm(i,r)*M(j,5);     
        end
        suma=sum(a);    % Denominator in brackets in Han and Davidson (2012, Eq. 1)
        for j=1:length(P)
            % Bracketed quantity in Han and Davidson (2012, Eq. 1)
            a(j)=a(j)/suma;
        end
        contribution=contribution+a;
    end 
end
contribution=contribution./sitenum/numofr;  % Contribution Cj
M(:,6)=contribution.';
Msort=sortrows(M,-6);   % sort M based on contribution in descending order
% Compute cumulative contribution for each earthquake when sorted in M(:,7)
Msort(1,7)=Msort(1,6);
out=0;
maxj=length(M);
for j=2:length(M)
    Msort(j,7)=Msort(j-1,7)+Msort(j,6);
    if Msort(j,7)>contribpct && out==0
        maxj=j-1;               % maxj is the last EQ to keep in RedM
        out=1;
    end
end
RedM=Msort(1:maxj,:);      % Update M to include only selected eq scenarios 
RedP=P(1,RedM(:,1));       % Update P to include only selected eq scenarios

% Plot cumulative contribution vs number of eqs
numeqs=1:length(M);      % Vector counting number of eqs
figure;
plot(numeqs,Msort(:,7));
ylim([0 1]);
title('Cumulative contribution vs number of earthquakes');
ylabel('Cumulative contribution to hazard curves');
xlabel('Number of earthquake scenarios');

% Plot cumulative contribution vs number of eqs for only those selected eqs
numeqs=1:maxj;      % Vector counting number of eqs
figure;
plot(numeqs,Msort(1:maxj,7));
ylim([0 1]);
title('Cumulative contribution vs number of earthquakes');
ylabel('Cumulative contribution to hazard curves');
xlabel('Number of earthquake scenarios');
end

