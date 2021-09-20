function [units, unitAmt] = unitStats_allTT(filenames)
%Calculate average RF for each condition and saccade discrmination index
%for each unit on an electrode (all the tetrodes).
%(This function is modified based on function 'unitStats'.) 
% 
% Inputs:
%     filename - The filenames of desired 3D pose files. 
%     
% Outputs:
%     units.ttID - The tetrode number. 
%     units.unitID - The unit number.
%     units.ave - Average RFs of all the conditions.
%     units.sdi - Saccade discrimination index of the unit.
%     units.polarVec - Vectors of average RFs of all conditions.
%     units.polarSum - Sum of units.polarVec. Indicate the preferred 
%         direction.
%     units.sumLength - Length of units.polarSum.
%     units.sumAngle - Angle of units.polarSum.
%     units.polarNorm - Normalized units.polarSum. Comparable across units.
% 
%     unitAmt - The amount fo units.
% 
% See also:
% 
% Author: Rosenberg Lab
% email: ari.rosenberg@wisc.edu
% Website: https://neuro.wisc.edu/staff/rosenberg-ari/
% Created: Sept 06 2021, ZKZ
% Editting history: 
% 06-Sep-2021, ZKZ: Created the function;
% 09-Sep-2021, ZKZ: Added the von Mises fit;
% 14-Sep-2021, ZKZ: Convert 'units.sumAngle' from radian to degree;
%                   added code to tell if 'isSig_ANOVA' is significant;
 
%------------- BEGIN CODE --------------
%% 
% 'unitAmt' is the amount of usable units
unitAmt = size(filenames, 1);

tilts = linspace(0, 7 / 4 * pi, 8);

% figure;
for i = 1:unitAmt
    % Load datafile.
    load(filenames(i));
    fprintf('filename: %s \n',(filenames(i)));
    
    % Extract the tetrode and unit IDs from its filename.
    sp = strsplit(filenames(i), {'TT','Unit'});
    sp_t = convertStringsToChars(sp(2));
    sp_u = convertStringsToChars(sp(3));
    % Save the tetrode and unit IDs.
    units(i).ttID = sp_t(1);
    units(i).unitID = sp_u(1);
    if saccadeData(1).isSig_ANOVA <= 0.05
        units(i).sig = 1;
    else
        units(i).sig = 0;
    end
    % Sort all the data by conditions, and save them into cell 'trials'.
    trials = cell(1, 8);
    for j = 1:size(saccadeData, 2)
        col_id = saccadeData(j).tilt / 45 + 1;
        trials{1, col_id} = [trials{1, col_id}, saccadeData(j).SaccadeFR];
    end
    
    % Get the average FR and sum of squared errors of each condition.
    errors = [];
    sumSqErrors = [];
    units(i).ave = NaN(1, 8);
    for j = 1:8
%         disp(trials{1, j});
        units(i).ave(1,j) = mean(trials{1, j}, 'all');

        errors = trials{1, j} - units(i).ave(1, j);
        sumSqErrors(1, j) = sum(errors.^2);
    end
    
    % Calculate the saccade discrimination index.
    rmax = max(units(i).ave);
    rmin = min(units(i).ave);
    sse = sum(sumSqErrors);
    units(i).sdi = (rmax - rmin) / (rmax - rmin + ...
        2 * sqrt(sse / (size(saccadeData, 2) - 8)));
    
    % Calculate the average vectors of each conditions in polar
    % coordinates with Euler's formula.
    units(i).polarVec = units(i).ave .* exp(1i * tilts);
    % Sum the vectors up. 
    units(i).polarSum = sum(units(i).polarVec);
    % Get the rho and theta of the sum vector. 
    units(i).sumLength = norm(units(i).polarSum);
    % BE CAREFUL WITH THE INPUT ORDER!!! In atan2, if we want the angle of 
    % (x, y), the input order should be (y, x). 
    units(i).sumAngle = atan2(imag(units(i).polarSum), ...
        real(units(i).polarSum));
    % Make all the preferred direction positive values. If that angles 
    % greater than pi were set as negative values, the angle coordinate 
    % would be incontinuous. This will create mistakes when we compute the
    % difference btw preferred directions later. 
    if units(i).sumAngle < 0
        units(i).sumAngle = units(i).sumAngle + 2 * pi;
    end
    % Convert angles from radian to degree.
    units(i).sumAngle = units(i).sumAngle * 180 / pi;
    
%     % Plot the tuning curve and the sum vector in polar coordinates.
%     subplot(2,round(unitAmt / 2), i);
%     
%     
% %     Plot the von Mise fit. 
% %      Set up inputs for function 'runVonMises4Tuning'. 
%     plotTitle = strcat('Tuning curve and von Mises fit of TT', units(i).ttID, ...
%         ' unit', units(i).unitID);
%     title(plotTitle);
    inData = units(i).ave;
    inBounds = NaN(2, 4);
    inBounds(:, 3) = [0; 18];
    inBounds(:, 4) = [0; 360];
    Tilts = linspace(0, 315, 8);
    
    % Calculate the von Mises fit. 
    [FitParameters, CIs] = runVonMises4Tuning(inData,inBounds,Tilts);
    vmFit(i).FitParameters = FitParameters;
    vmFit(i).CIs = CIs;
    
    units(i).vmMiu = FitParameters(4);
    units(i).vmKappa = FitParameters(3);
    % Make all the preferred direction positive values. If that angles 
    % greater than 360 were set as negative values, the angle coordinate 
    % would be incontinuous. This will create mistakes when we compute the
    % difference btw preferred directions later. 
    if units(i).vmMiu < 0
        units(i).vmMiu = units(i).vmMiu + 360;
    end
    
    % Normalize each sum vector.
    units(i).polarNorm = (units(i).polarSum / sum(units(i).ave - min(units(i).ave)))';

%     if units(i).vmKappa > 4 && norm(units(i).polarNorm) < 0.25
%         figure();
%         % Plot the tuning curve.
%         thetas = linspace(0, 2 * pi, 9);
%         rhos = [units(i).ave, units(i).ave(1)];
%         polarplot(thetas, rhos, 'x-', 'LineWidth', 1.5); hold on
%         
%         % Plot the sum vector.
%         polarplot([0, units(i).sumAngle * pi / 180], [0, units(i).sumLength], 'r', ...
%             'LineWidth', 2);
%         
%         % Plot the von Mises fit.
%         vonMisesFun = @(Params,xdata) Params(1)+Params(2)*exp(Params(3)*(cosd(xdata-Params(4))-1));
%         Params = FitParameters;
%         theta = linspace(0, 2 * pi, 200);
%         rho = vonMisesFun(Params,theta * 180 / pi);
%         peak = vonMisesFun(FitParameters, FitParameters(4));
%         polarplot(theta, rho, 'LineWidth', 2); hold on
%         
%         % Plot the peak of von Mises. 
%         polarplot([0, FitParameters(4) * pi / 180], [0, peak], 'k','LineWidth', 2);
%         legend('Tuning curve', 'Sum vector', 'Von Mises fit', 'Von Mises peak');
%         Kappa = num2str(units(i).vmKappa);
%         Mag = num2str(units(i).sumLength);
%         Norm_mag = num2str(norm(units(i).polarNorm));
%         SDI = num2str(norm(units(i).sdi));
%         plotTitle = ['Kappa = ', Kappa, ', Mag = ', Mag, ', Norm mag = ', Norm_mag,  ', SDI = ', SDI];
%         title(plotTitle);
%     end

end

% % Plot the normalized vector in the same coordinates. 
% subplot(2,fix(unitAmt / 2) + 1, unitAmt + 1); 
% for i = 1:unitAmt
%     polarplot([0, atan2(imag(units(i).polarNorm), ...
%         real(units(i).polarNorm))], [0, norm(units(i).polarNorm)],...
%         'LineWidth', 2);hold on
% end
% title('Nomalized neural vectors');
% hold off

end

%------------- END OF CODE --------------