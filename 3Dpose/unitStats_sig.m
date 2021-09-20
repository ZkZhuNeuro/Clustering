function [units, unitAmt, sigAmt] = unitStats_sig(filenames)


%% 
% 'unitAmt' is the amount of usable units
unitAmt = size(filenames, 1);
sigAmt = unitAmt;

tilts = linspace(0, 7 / 4 * pi, 8);

for i = 1:unitAmt
    
    % Load datafile.
    load(filenames(i), 'saccadeData');
    fprintf('filename: %s \n',(filenames(i)));
    
    % Extract the tetrode and unit IDs from its filename.
    sp = strsplit(filenames(i), {'TT','Unit'});
    sp_t = convertStringsToChars(sp(2));
    sp_u = convertStringsToChars(sp(3));
    % Save the tetrode and unit IDs.
    
    if saccadeData(1).isSig_ANOVA <= 0.05
        units(i).ttID = sp_t(1);
        units(i).unitID = sp_u(1);
        units(i).sig = 1;
        
        % Sort all the data by conditions, and save them into cell
        % 'trials'.
        trials = cell(1, 8);
        for j = 1:size(saccadeData, 2)
            col_id = saccadeData(j).tilt / 45 + 1;
            trials{1, col_id} = [trials{1, col_id}, ...
                saccadeData(j).SaccadeFR];
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
        % BE CAREFUL WITH THE INPUT ORDER!!! In atan2, if we want the angle
        % of (x, y), the input order should be (y, x).
        units(i).sumAngle = atan2(imag(units(i).polarSum), ...
            real(units(i).polarSum));
        % Make all the preferred direction positive values. If that angles
        % greater than pi were set as negative values, the angle coordinate
        % would be incontinuous. This will create mistakes when we compute 
        % the difference btw preferred directions later.
        if units(i).sumAngle < 0
            units(i).sumAngle = units(i).sumAngle + 2 * pi;
        end
        % Convert angles from radian to degree.
        units(i).sumAngle = units(i).sumAngle * 180 / pi;
        
        % Normalize each sum vector. The 'new' nm. 
        units(i).polarNorm = (units(i).polarSum / sum(units(i).ave -...
            min(units(i).ave)))';
        
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
    else 
        sigAmt = sigAmt - 1;
        units(i).unitID = 0;
        fprintf('Tuning is not significant! \n')
    end
    
end

end