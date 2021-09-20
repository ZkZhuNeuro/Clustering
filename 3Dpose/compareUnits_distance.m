function [linearCorr] = compareUnits_distance(units, unitAmt)
%

figure;
%% Plot the normalized vector in the same coordinates.
subplot(2, 2, 1);
for i = 1:unitAmt
    polarplot([0, atan2(imag(units(i).polarNorm), ...
        real(units(i).polarNorm))], [0, norm(units(i).polarNorm)],...
        'LineWidth', 2); hold on
end
title('Nomalized neural vectors');
hold off

%% Plot the differences among units in a rose plot (polarhistogram). Then, calculated their linear correlations. 
theta = NaN(1, (unitAmt - 1) * unitAmt / 2);
linearCorr = NaN((unitAmt - 1) * unitAmt / 2, 2);
d_sdi = NaN(1, (unitAmt - 1) * unitAmt / 2);
d_mag = NaN(1, (unitAmt - 1) * unitAmt / 2);
k = 0;
% A double for-loop is used here for comparion. 

d_sdi = cell(1, 8);
d_mag = cell(1, 8);
d_dir = cell(1, 8);
for i = 1:(unitAmt - 1)
    for j = (i + 1):unitAmt
        k = k + 1;
        % Put all the angles of vectors into matrix theta. 
        theta(1, k) = abs(units(i).sumAngle - units(j).sumAngle);
        [r, p] = corr(units(i).ave',units(j).ave');
        linearCorr(k, :) = [r, p];
        fprintf('Correlation coeffient of TT %s unit %s and TT %s unit %s is %.3f, p-value is %.3f; \n', ...
            units(i).ttID, units(i).unitID, units(j).ttID, ...
            units(j).unitID, r, p);
        
        dis = abs(str2num(units(i).ttID) - str2num(units(j).ttID));
        if dis == 0
            d_sdi{1, 1} = [d_sdi{1, 1}, abs(units(i).sdi - units(j).sdi)];
            d_mag{1, 1} = [d_mag{1, 1}, abs(norm(units(i).polarNorm) ...
                - norm(units(j).polarNorm))];
            d_dir{1, 1} = [d_dir{1, 1}, abs(units(i).sumAngle - ...
                units(j).sumAngle)];
        elseif dis == 1
            d_sdi{1, 2} = [d_sdi{1, 2}, abs(units(i).sdi - units(j).sdi)];
            d_mag{1, 2} = [d_mag{1, 2}, abs(norm(units(i).polarNorm) ...
                - norm(units(j).polarNorm))];
            d_dir{1, 2} = [d_dir{1, 2}, abs(units(i).sumAngle - ...
                units(j).sumAngle)];
        elseif dis == 2
            d_sdi{1, 3} = [d_sdi{1, 3}, abs(units(i).sdi - units(j).sdi)];
            d_mag{1, 3} = [d_mag{1, 3}, abs(norm(units(i).polarNorm) ...
                - norm(units(j).polarNorm))];
            d_dir{1, 3} = [d_dir{1, 3}, abs(units(i).sumAngle - ...
                units(j).sumAngle)];
        elseif dis == 3
            d_sdi{1, 4} = [d_sdi{1, 4}, abs(units(i).sdi - units(j).sdi)];
            d_mag{1, 4} = [d_mag{1, 4}, abs(norm(units(i).polarNorm) ...
                - norm(units(j).polarNorm))];
            d_dir{1, 4} = [d_dir{1, 4}, abs(units(i).sumAngle - ...
                units(j).sumAngle)];
        elseif dis == 4
            d_sdi{1, 5} = [d_sdi{1, 5}, abs(units(i).sdi - units(j).sdi)];
            d_mag{1, 5} = [d_mag{1, 5}, abs(norm(units(i).polarNorm) ...
                - norm(units(j).polarNorm))];
            d_dir{1, 5} = [d_dir{1, 5}, abs(units(i).sumAngle - ...
                units(j).sumAngle)];
        elseif dis == 5
            d_sdi{1, 6} = [d_sdi{1, 6}, abs(units(i).sdi - units(j).sdi)];
            d_mag{1, 6} = [d_mag{1, 6}, abs(norm(units(i).polarNorm) ...
                - norm(units(j).polarNorm))];
            d_dir{1, 6} = [d_dir{1, 6}, abs(units(i).sumAngle - ...
                units(j).sumAngle)];
        elseif dis == 6
            d_sdi{1, 7} = [d_sdi{1, 7}, abs(units(i).sdi - units(j).sdi)];
            d_mag{1, 7} = [d_mag{1, 7}, abs(norm(units(i).polarNorm) ...
                - norm(units(j).polarNorm))];
            d_dir{1, 7} = [d_dir{1, 7}, abs(units(i).sumAngle - ...
                units(j).sumAngle)];
        elseif dis == 7
            d_sdi{1, 8} = [d_sdi{1, 8}, abs(units(i).sdi - units(j).sdi)];
            d_mag{1, 8} = [d_mag{1, 8}, abs(norm(units(i).polarNorm) ...
                - norm(units(j).polarNorm))];
            d_dir{1, 8} = [d_dir{1, 8}, abs(units(i).sumAngle - ...
                units(j).sumAngle)];
        end
        
%         if dis == 0
%             d_mag{1, 1} = [d_mag{1, 1}, abs(norm(units(i).polarNorm) ...
%                 - norm(units(j).polarNorm))];
%         elseif dis == 1
%             d_mag{1, 2} = [d_mag{1, 2}, abs(norm(units(i).polarNorm) ...
%                 - norm(units(j).polarNorm))];
%         elseif dis == 2
%             d_mag{1, 3} = [d_mag{1, 3}, abs(norm(units(i).polarNorm) ...
%                 - norm(units(j).polarNorm))];
%         elseif dis == 3
%             d_mag{1, 4} = [d_mag{1, 4}, abs(norm(units(i).polarNorm) ...
%                 - norm(units(j).polarNorm))];
%         elseif dis == 4
%             d_mag{1, 5} = [d_mag{1, 5}, abs(norm(units(i).polarNorm) ...
%                 - norm(units(j).polarNorm))];
%         elseif dis == 5
%             d_mag{1, 6} = [d_mag{1, 6}, abs(norm(units(i).polarNorm) ...
%                 - norm(units(j).polarNorm))];
%         elseif dis == 6
%             d_mag{1, 7} = [d_mag{1, 7}, abs(norm(units(i).polarNorm) ...
%                 - norm(units(j).polarNorm))];
%         elseif dis == 7
%             d_mag{1, 8} = [d_mag{1, 8}, abs(norm(units(i).polarNorm) ...
%                 - norm(units(j).polarNorm))];
%         end
        
    end
end
subplot(2, 2, 2);
% edges_rose = linspace(0, 2 * pi, 13);
% polarhistogram(theta,'BinEdges', edges_rose + 15 / 180 * pi);
for i = 1:8
%     figure;
    scatter(300 * (i - 1) * ones(1, size(d_dir{1, i}, 2)), ...
        d_dir{1, i} * 180 / pi, 'LineWidth', 1); hold on
    scatter(300 * (i - 1), mean(d_dir{1, i} * 180 / pi), '_k', 'LineWidth', 2);
end
xticks(linspace(0, 2100, 8));
title('Distribution of differences of preferred directions');
xlabel('Distance');
ylabel('ΔDirection');


%% Plot the distributions of differences of SDIs and normalized magnitudes. 

subplot(2, 2, 3);
for i = 1:8
%     figure;
    scatter(300 * (i - 1) * ones(1, size(d_sdi{1, i}, 2)), d_sdi{1, i}, ...
        'LineWidth', 1); hold on
%     boxplot(d_sdi{1, i}, 300 * (i - 1)); 
    hold on
    scatter(300 * (i - 1), mean(d_sdi{1, i}), '_k', 'LineWidth', 2);
end
title('Distribution of differences of saccade discrimination indices');
xticks(linspace(0, 2100, 8));
% ylim([0 1]);
xlabel('Distance');
ylabel('ΔSDI');

subplot(2, 2, 4);
for i = 1:8
%     figure;
    scatter(300 * (i - 1) * ones(1, size(d_mag{1, i}, 2)), d_mag{1, i}, ...
        'LineWidth', 1); hold on
    scatter(300 * (i - 1), mean(d_mag{1, i}), '_k', 'LineWidth', 2);
end
title('Distribution of differences of normalized magnitudes');
xticks(linspace(0, 2100, 8));
% ylim([0 1]);
xlabel('Distance');
ylabel('Δmagnitude');

end