clear all
close all
clc;
a = 'Frank'; b = 'V3A'; c = '29May2019'; d = 'TT';
[monkeyName, area, date, tt] = whichTT(a,b,c,d);
[filenames] = fileSeeker(monkeyName, area, date, tt);
[units, unitAmt] = unitStats_allTT(filenames);
linearCorr = compareUnits(units, unitAmt);

% figure;
% for i = 1:unitAmt
%     inData = units(i).ave;
%     inBounds = NaN(2, 4);
%     inBounds(:, 3) = [0; 18];
%     inBounds(:, 4) = [0; 360];
%     Tilts = linspace(0, 315, 8);
%     [FitParameters, CIs] = runVonMises4Tuning(inData,inBounds,Tilts);
%     vmFit(i).FitParameters = FitParameters;
%     vmFit(i).CIs = CIs;
%     
%     vonMisesFun = @(Params,xdata) Params(1)+Params(2)*exp(Params(3)*(cosd(xdata-Params(4))-1));
%     Params = FitParameters;
%     theta = linspace(0, 2 * pi, 200);
%     rho = vonMisesFun(Params,theta * 180 / pi);
%     polarplot(theta, rho); hold on
% %     polarplot([0, atan2(imag(units(i).polarNorm), real(units(i).polarNorm))], [0, units(i).sumLength],...
% %         'LineWidth', 2);hold on
% end
% title('von Mises fits');