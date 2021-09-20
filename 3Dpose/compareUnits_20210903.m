function [unit_vec] = compareUnits_20210903(filenames)
% Measure the correlation among units.
% Firstly, collect all the SaccadeFR data from the 3Dpose files.
% Secondly, compute the averages of SaccadeFRs for each condition.
% Thirdly, measure the linear correlation among the average signals.

%% Get SaccadeFR data from 3Dpose files.

% 'unit_str' is for naming.
unit_str = []; 
% 'units' is a matrix which contains all the SaccadeFR data in the trial
% order (which is pseudo-random).
units = [];
% 'units_many' is the amount of usable units
units_many = size(filenames, 1) - 1;

for i = 2:size(filenames, 1)
    
    % unit_sum is a temporary array that saves the sum of the FR in
    % each condition
    unit_sum = zeros(1, 8);
    % unit_sum is a temporary array that saves all the FR in
    % each condition
    unit_i = [];
    
    sp = strsplit(filenames(i), 'Unit');
    sp_c = convertStringsToChars(sp(2));
    % Extract the unit number from its filename.
    % Note that the unit numbers 'Unit_str(i)' are still str here.
    load(filenames(i));
    fprintf('filename: %s \n',(filenames(i)));
    unit_str = [unit_str, sp_c(1)];
    sectNum = size(saccadeData, 2) / 8;
    for j = 1:size(saccadeData, 2)
        sect = fix((j - 1)/8) + 1;
        % Append a SaccadeFR to 'unit_i'.
        unit_i = [unit_i, saccadeData(j).SaccadeFR];
        
        % In order to get the average FR, I go through all the data first,
        % and put all the data into struct 'unit_mat' in the order of 
        % conditions. 
        if saccadeData(j).tilt == 0
            unit_mat(sect).tilt1 = saccadeData(j).SaccadeFR;
        elseif saccadeData(j).tilt == 45
            unit_mat(sect).tilt2 = saccadeData(j).SaccadeFR;
        elseif saccadeData(j).tilt == 90
            unit_mat(sect).tilt3 = saccadeData(j).SaccadeFR;
        elseif saccadeData(j).tilt == 135
            unit_mat(sect).tilt4 = saccadeData(j).SaccadeFR;
        elseif saccadeData(j).tilt == 180
            unit_mat(sect).tilt5 = saccadeData(j).SaccadeFR;
        elseif saccadeData(j).tilt == 225
            unit_mat(sect).tilt6 = saccadeData(j).SaccadeFR;
        elseif saccadeData(j).tilt == 270
            unit_mat(sect).tilt7 = saccadeData(j).SaccadeFR;
        elseif saccadeData(j).tilt == 315
            unit_mat(sect).tilt8 = saccadeData(j).SaccadeFR;
        end
    end
    
    % Get the average of each tilt. Save them in struct 'units_ave'.
    for j = 1:8
        units_ave(i-1).tilt1 = sum([unit_mat(:).tilt1]) / sectNum;
        units_ave(i-1).tilt2 = sum([unit_mat(:).tilt2]) / sectNum;
        units_ave(i-1).tilt3 = sum([unit_mat(:).tilt3]) / sectNum;
        units_ave(i-1).tilt4 = sum([unit_mat(:).tilt4]) / sectNum;
        units_ave(i-1).tilt5 = sum([unit_mat(:).tilt5]) / sectNum;
        units_ave(i-1).tilt6 = sum([unit_mat(:).tilt6]) / sectNum;
        units_ave(i-1).tilt7 = sum([unit_mat(:).tilt7]) / sectNum;
        units_ave(i-1).tilt8 = sum([unit_mat(:).tilt8]) / sectNum;
    end
    % Pile 'unit_i's up to get the 'units'. 
    if ~isempty(units)
        units = cat(1, units, unit_i);
    else
        units = [units, unit_i];
    end
    
    % Compute the tilt discrimination index (TDI).
    % Find the maximum of averages.
    rmax = max([units_ave(i-1).tilt1, units_ave(i-1).tilt2, ...
        units_ave(i-1).tilt3, units_ave(i-1).tilt4, units_ave(i-1).tilt5, ...
        units_ave(i-1).tilt6, units_ave(i-1).tilt7, units_ave(i-1).tilt8]);
    
    % Find the minimum of averages.
    rmin = min([units_ave(i-1).tilt1, units_ave(i-1).tilt2, ...
        units_ave(i-1).tilt3, units_ave(i-1).tilt4, units_ave(i-1).tilt5, ...
        units_ave(i-1).tilt6, units_ave(i-1).tilt7, units_ave(i-1).tilt8]);
    
    % Calculate SSE.
    sse = 0;
    for j = 1:sectNum
        square_error(j).tilt1 = (unit_mat(j).tilt1 - units_ave(i-1).tilt1)^2;
        square_error(j).tilt2 = (unit_mat(j).tilt2 - units_ave(i-1).tilt2)^2;
        square_error(j).tilt3 = (unit_mat(j).tilt3 - units_ave(i-1).tilt3)^2;
        square_error(j).tilt4 = (unit_mat(j).tilt4 - units_ave(i-1).tilt4)^2;
        square_error(j).tilt5 = (unit_mat(j).tilt5 - units_ave(i-1).tilt5)^2;
        square_error(j).tilt6 = (unit_mat(j).tilt6 - units_ave(i-1).tilt6)^2;
        square_error(j).tilt7 = (unit_mat(j).tilt7 - units_ave(i-1).tilt7)^2;
        square_error(j).tilt8 = (unit_mat(j).tilt8 - units_ave(i-1).tilt8)^2;
        sse = sse + square_error(j).tilt1 + square_error(j).tilt2 ...
            + square_error(j).tilt3 + square_error(j).tilt4 ...
            + square_error(j).tilt5 + square_error(j).tilt6 ...
            + square_error(j).tilt7 + square_error(j).tilt8;
    end
    
    % Calculate TDI, and save it to struct 'unit_vec'.
    unit_vec(i-1).TDI = (rmax - rmin) / (rmax - rmin + ...
        2 * sqrt(sse / (size(saccadeData, 2) - 8)));
end
% 
% disp(units_ave);
% plot(units_ave');

%% Generate tuning curves and neural vectors for the units. 

% Create a struct 'unit_vec' which contains all the information of the
% tuning curves
for i = 1:units_many
    % 'unit_length_total' is a temporary variable that saves the sum of
    % vectors' length. 
    unit_length_total = sum([units_ave(i).tilt1, units_ave(i).tilt2, ...
        units_ave(i).tilt3, units_ave(i).tilt4, units_ave(i).tilt5, ...
        units_ave(i).tilt6, units_ave(i).tilt7, units_ave(i).tilt8]);
    
    % 'tilt_x/y' is the vector's x/y value in Cartesian coordinates. 
    unit_vec(i).tilt1_x = units_ave(i).tilt1 * cosd(0);
    unit_vec(i).tilt1_y = units_ave(i).tilt1 * sind(0);
    unit_vec(i).tilt2_x = units_ave(i).tilt2 * cosd(45);
    unit_vec(i).tilt2_y = units_ave(i).tilt2 * sind(45);
    unit_vec(i).tilt3_x = units_ave(i).tilt3 * cosd(90);
    unit_vec(i).tilt3_y = units_ave(i).tilt3 * sind(90);
    unit_vec(i).tilt4_x = units_ave(i).tilt4 * cosd(135);
    unit_vec(i).tilt4_y = units_ave(i).tilt4 * sind(135);
    unit_vec(i).tilt5_x = units_ave(i).tilt5 * cosd(180);
    unit_vec(i).tilt5_y = units_ave(i).tilt5 * sind(180);
    unit_vec(i).tilt6_x = units_ave(i).tilt6 * cosd(225);
    unit_vec(i).tilt6_y = units_ave(i).tilt6 * sind(225);
    unit_vec(i).tilt7_x = units_ave(i).tilt7 * cosd(270);
    unit_vec(i).tilt7_y = units_ave(i).tilt7 * sind(270);
    unit_vec(i).tilt8_x = units_ave(i).tilt8 * cosd(315);
    unit_vec(i).tilt8_y = units_ave(i).tilt8 * sind(315);
    
    % 'tilt_sum_x/y' is the sum value of x/y in Cartesian coordinates. 
    unit_vec(i).tilt_sum_x = unit_vec(i).tilt1_x ...
        + unit_vec(i).tilt2_x + unit_vec(i).tilt3_x ...
        + unit_vec(i).tilt4_x + unit_vec(i).tilt5_x ...
        + unit_vec(i).tilt6_x + unit_vec(i).tilt7_x ...
        + unit_vec(i).tilt8_x;
    unit_vec(i).tilt_sum_y = unit_vec(i).tilt1_y ...
        + unit_vec(i).tilt2_y + unit_vec(i).tilt3_y ...
        + unit_vec(i).tilt4_y + unit_vec(i).tilt5_y ...
        + unit_vec(i).tilt6_y + unit_vec(i).tilt7_y ...
        + unit_vec(i).tilt8_y;
    
    % Plot the tuning curve and the sum vector in polar coordinates. 
    subplot(2,3,i);
    
    % Plot the tuning curve.
    thetas = linspace(0, 2 * pi, 9);
    rhos = [units_ave(i).tilt1, units_ave(i).tilt2, ...
        units_ave(i).tilt3, units_ave(i).tilt4, units_ave(i).tilt5, ...
        units_ave(i).tilt6, units_ave(i).tilt7, units_ave(i).tilt8, ... 
        units_ave(i).tilt1];
    polarplot(thetas, rhos, 'x-');hold on
    plotTitle = strcat('Sum of neural vectors in unit', unit_str(i));
    title(plotTitle);
    
    % Plot the sum vector. 
    [t, r] = cart2pol(unit_vec(i).tilt_sum_x, unit_vec(i).tilt_sum_y);
    polarplot([0, t], [0, r], 'LineWidth', 2)
    hold off
    
   
%     unit_vec(i).unit_length = sqrt((unit_vec(i).tilt_sum_x)^2 ...
%         + (unit_vec(i).tilt_sum_y)^2);
%     
%     unit_vec(i).unit_length_normalized = unit_vec(i).unit_length ...
%         / unit_length_total;

    % The next step is to plot the normalized vectors. 
    % 'scaled_x/y' is the normalized x/y value in Cartesian coordinates. 
    % (Here, scaled_x is equal to sum of vector x divided by sum of total
    % lengths, instead of sum of vector x divided by sum of x lengths,
    % because we want to keep the original directions, and only adjust the
    % lengths of the vectors. If we do later, x and y directions will be
    % scaled by different factors, which would change the direction of the
    % vector.)
    unit_vec(i).scaled_x = unit_vec(i).tilt_sum_x ...
        / unit_length_total;
    unit_vec(i).scaled_y = unit_vec(i).tilt_sum_y ...
        / unit_length_total;
end

%% Plot the normalized neural vectors. 
subplot(2,3,units_many + 1)
for i = 1:units_many
    [t, r] = cart2pol(unit_vec(i).scaled_x, unit_vec(i).scaled_y);
    polarplot([0, t], [0, r], 'LineWidth', 2);hold on
end

end