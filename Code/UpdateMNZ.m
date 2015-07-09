function [RedM2,RedP2]=UpdateMNZ(RedM1,RedP1)

% This function updates the M and P structures to keep only those eq 
% scenarios and exceedence probabilities with Pj>0 (where Pj's are assigned 
% eq scenario annual probabilities output from optimization)

% Read in Pj's output by optimization
% First row has number of candidates in optimization run
% Second row has allowable number of eqs in reduced set 
% Remaining rows have Pj values for j=1 to Numcand

Numcand=dlmread('Pjs.txt','',[0 0 0 0]);     
Jred=dlmread('Pjs.txt','',[1 0 1 0]);
Pjs=dlmread('Pjs.txt','',2,0);

% Update annual occurrence probabilities for retained eq scenarios to be
%  the Pj values

RedM2=RedM1(Pjs>0,:);
RedM2(:,5)=Pjs(Pjs>0);
RedP2=RedP1(Pjs>0);

end
