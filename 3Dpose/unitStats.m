function [units, unitAmt] = unitStats(filenames)
%Calculate average RF for each condition and saccade discrmination index
%for each unit.
% 
% Inputs:
%     filename - The filenames of desired 3D pose files. 
%     
% Outputs:
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
% Created: Sept 04 2021, ZKZ
% Editting history: 
% 04-Sep-2021, ZKZ: Created the function;
% 05-Sep-2021, ZKZ: Corrected the mistake that discarded the first unit
%     data; Added an output 'unitAmt'; Commented the normalized vector 
%     plot; Made all the preferred directions to be positve angles. 
 
%------------- BEGIN CODE --------------
%% 
% 'unitAmt' is the amount of usable units
unitAmt = size(filenames, 1);

tilts = linspace(0, 7 / 4 * pi, 8);

for i = 1:unitAmt
    % Load datafile.
    load(filenames(i));
    fprintf('filename: %s \n',(filenames(i)));
    
    % Extract the unit number from its filename.
    % Save the unit ID.
    sp = strsplit(filenames(i), 'Unit');
    sp_c = convertStringsToChars(sp(2));
    units(i).unitID = sp_c(1);
    
    sectNum = size(saccadeData, 2) / 8;
    
    % Sort all the data by conditions, and save them into matrix 'trials'.
    trials = zeros(sectNum, 8);
    for j = 1:size(saccadeData, 2)
        col_id = saccadeData(j).tilt / 45 + 1;
        row_id = saccadeData(j).BlockNum;
        trials(row_id, col_id) = saccadeData(j).SaccadeFR;
    end
    
    % Get the average FR of each condition.
    units(i).ave = mean(trials);
    
    % Calculate the saccade discrimination index.
    rmax = max(units(i).ave);
    rmin = min(units(i).ave);
    errors = trials - repmat(units(i).ave, 7, 1);
    sqErrors = errors.^2;
    sse = sum(sqErrors, 'all');
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
    % Make all the preferred direction positive values. If angles greater
    % than pi were set as negative values, the angle coordinate would be
    % incontinuous. This will create mistakes when we compute the
    % difference btw preferred directions later. 
    if units(i).sumAngle < 0
        units(i).sumAngle = units(i).sumAngle + 2 * pi;
    end
    
    % Plot the tuning curve and the sum vector in polar coordinates. 
    subplot(2,round(unitAmt / 2), i);
    
    % Plot the tuning curve. 
    thetas = linspace(0, 2 * pi, 9);
    rhos = [units(i).ave, units(i).ave(1)];
    polarplot(thetas, rhos, 'x-', 'LineWidth', 1.5); hold on
    plotTitle = strcat('Tuning curve of unit', units(i).unitID);
    title(plotTitle)
    % Plot the sum vector.
    polarplot([0, units(i).sumAngle], [0, units(i).sumLength], ...
        'LineWidth', 2); hold off
    
    % Normalize each sum vector.
    units(i).polarNorm = units(i).polarSum / sum(units(i).ave)';
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