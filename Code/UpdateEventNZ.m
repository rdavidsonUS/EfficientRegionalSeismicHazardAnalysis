function [RedEvent]=UpdateEventNZ(Event)

% This function takes the results of optimization applied to ground motion
% maps (2nd implementation of optimization) and produces the final set of 
% maps (Event) (only those with Pj>0

% Read in Pj's output by optimization
% First row has number of candidates in optimization run
% Second row has allowable number of eqs in reduced set 
% Remaining rows have Pj values for j=1 to Numcand

Numcand=dlmread('Pjs2.txt','',[0 0 0 0]);     
Jred=dlmread('Pjs2.txt','',[1 0 1 0]);
Pjs=dlmread('Pjs2.txt','',2,0);

% Update annual occurrence probabilities for retained ground motion maps to 
% be the Pj values

count=length(Event);

k=1;
for i=1:count
    if Pjs(i)>0
        RedEvent(k)=Event(i);
        RedEvent(k).AP=Pjs(i);
        k=k+1;
    end
end

end
