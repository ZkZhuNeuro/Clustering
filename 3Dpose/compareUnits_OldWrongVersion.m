function [sigComp] = compareUnits_OldWrongVersion(filenames)
% Measure the correlation among units. 
    % Firstly, collect all the SaccadeFR data from the 3Dpose files. 
    % Secondly, compute the averages of SaccadeFRs for each condition. 
    % Thirdly, measure the linear correlation among the average signals. 

    %% Get SaccadeFR data from 3Dpose files. 
    unit_str = []; % unit_str is for naming. 
    units = []; % units is a matrix which saves all the SaccadeFR data. 
    for i = 1:size(filenames, 1)
        sp = strsplit(filenames(i), 'Unit');
        sp_c = convertStringsToChars(sp(2));
        
        % Discard the first unit, which is noise. 
        if str2num(sp_c(1)) > 1 
            
            % Extract the unit number from its filename. 
            % Note that the unit numbers 'Unit_str(i)' are still str here. 
            load(filenames(i));
            fprintf('filename: %s \n',(filenames(i)));
            unit_i = [];
            unit_str = [unit_str, sp_c(1), ' '];
            
            % Save each unit's saccadeFR to unit_i. 
            % After this for loop, unit_i contains all the saccadeFR data
            % of one data file. 
            for j = 1:size(saccadeData, 2)
                unit_i = [unit_i, saccadeData(j).SaccadeFR];
            end
            
        end
%             % Append the SaccadeFR to marix 'units' so that in the end of 
%             % the big for loop, 'units' is the matrix of all SaccadeFRs. 
%             if ~isempty(units)
%                 units = cat(1, units, unit_i);
%             else
%                 units = [units, unit_i];
%             end
%         end
%     end
% %     plot(units')
% 
%     %% Compare each pair of the units. 
%     
%     %Get the average signal matrix 'units_ave'.
%     units_ave = zeros(size(units, 1), 8);
%     for i = 1:size(units, 1)
%         allData = units(i, :);
%         sect = size(allData, 2) / 8;
%         matData = reshape(allData, 8, sect);
%         matData = matData';
%         aveData = mean(matData);
%         units_ave(i, :) = aveData;
%     end
%     fprintf('Average signals of all the units: \n');
%     disp(units_ave);
%     plot(units_ave');
%     
%     % Generate a pairwise list ('compL') for naming after the comparing. 
%     unit_n = str2num(unit_str);
%     compL = nchoosek(unit_n, 2);
%     
%     % Linear comparison.
%     % stats is a matrix that saves the comparison results. 
%     stats = zeros(size(compL, 1), 4);
%     k = 0;
%     for i = 1:(size(units_ave, 1) - 1)
%         for j = (i + 1):size(units_ave, 1)
%             k = k + 1;
%             unit1 = units_ave(i, :);
%             unit2 = units_ave(j, :);
%             
%             %Choose 'Pearson' / ''Kendall' / 'Spearman':
%             [r, p] = corr(unit1', unit2', 'Type', 'Pearson');
%             stats(k, :) = [compL(k, 1),compL(k, 2), r, p];
%         end
%     end
%     
%     % sigComp is a struct that saves the comparison results. 
%     for i = 1:size(stats, 1)
%         sigComp(i).unit1 = stats(i, 1);
%         sigComp(i).unit2 = stats(i, 2);
%         sigComp(i).r_value = stats(i, 3);
%         sigComp(i).p_value = stats(i, 4);
%     end
%     
end