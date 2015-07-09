function PlotHazCurve(site,HCD,trueRGM,returnperiods)

% This function plots a hazard curve for a specified site
% site is the site for which the hazard curve will be plotted
% HCD is teh hazard curve information for the run of interest
% trueRGM contains the "true" hazard curve points being matched
% returnperiods is vector of return periods trueRGM is available for

gm=0.001:0.001:3;

figure;
semilogy(gm, HCD(site,:));
hold all;
AEP=1./returnperiods;
plot(trueRGM(site,:),AEP,'o');
hold off;
title('Hazard curve for site with worst HCE vs true hazard curve points');
xlabel('Peak ground acceleration (g)');
ylabel('Annual probability of exceedence');

end
