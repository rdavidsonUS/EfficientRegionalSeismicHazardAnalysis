function cc=spatialcorre(Event)

% This function calculates the spatial correlation ?_(a,b) coefficient of 
% ground motion values at every pair of sites (a,b) in the set of ground 
% motion maps in Event

sitenum=length(Event(1,1).GM(:,1));     % Number of sites
count=length(Event);                    % Number of ground motion maps
for i=1:sitenum
    for j=1:count
        site(i,j)=Event(1,j).GM(i,1);
    end
end
site=site';
for i=1:sitenum
    for j=i:sitenum
        for k=1:count
            cc1(k)=Event(k).AP*site(k,i)*site(k,j);
            cc2(k)=Event(k).AP*site(k,i);
            cc3(k)=Event(k).AP*site(k,j);
            cc4(k)=Event(k).AP*site(k,i)*site(k,i);
            cc5(k)=Event(k).AP*site(k,j)*site(k,j);
        end
        cc(i,j)=(sum(cc1)-sum(cc2)*sum(cc3))/(sqrt(sum(cc4)-sum(cc2)*sum(cc2))*sqrt(sum(cc5)-sum(cc3)*sum(cc3)));
    end
end