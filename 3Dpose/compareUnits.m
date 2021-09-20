function [utc] = compareUnits(units, unitAmt, sigAmt)
%Prepare the data for pairwise comparison.
% Inputs:
%     units - This is the output of function 'unitStats'.
%     unitAmt - The number of units
%
% Outputs:
%     utc - A strcut that saves unit data for comparison.
%
% See also:
%
% Author: Rosenberg Lab
% email: ari.rosenberg@wisc.edu
% Website: https://neuro.wisc.edu/staff/rosenberg-ari/
% Created: Sept 06 2021, ZKZ
% Editting history:
% 06-Sep-2021, ZKZ: Created the function;
% 07-Sep-2021, ZKZ: Changed distributions of sdi and magnitudes to
%     distribution of differences of sdi and magnitudes.
% 14-Sep-2021, ZKZ: Combine all the output to a struct 'utc'.
% 18-Sep-2021, ZKZ: Add condition that: 
%     units(i).unitID > 0 && units(j).unitID > 0

%------------- BEGIN CODE --------------
% figure;
%% Plot the normalized vector in the same coordinates.
% subplot(2, 2, 1);
% for i = 1:unitAmt
%     polarplot([0, atan2(imag(units(i).polarNorm), ...
%         real(units(i).polarNorm))], [0, norm(units(i).polarNorm)],...
%         'LineWidth', 2); hold on
% end
% title('Nomalized neural vectors');
% % hold off

%% Plot the differences among units in a rose plot (polarhistogram). Then, calculated their linear correlations.
utc.d_dir = NaN(1, (sigAmt - 1) * sigAmt / 2);
utc.d_dir_vm = NaN(1, (sigAmt - 1) * sigAmt / 2);
% linearCorr = NaN((sigAmt - 1) * sigAmt / 2, 2);
utc.d_sdi = NaN(1, (sigAmt - 1) * sigAmt / 2);
utc.d_mag = NaN(1, (sigAmt - 1) * sigAmt / 2);
k = 0;
% A double for-loop is used here for comparion.
for i = 1:(unitAmt - 1)
    for j = (i + 1):unitAmt
        if units(i).unitID > 0 && units(j).unitID > 0
            k = k + 1;
            % Put all the angles of vectors into matrix theta.
            if units(i).sumAngle < units(j).sumAngle
                temp = units(i).sumAngle + 360;
                utc.d_dir(1, k) = temp - units(j).sumAngle;
            else
                utc.d_dir(1, k) = units(i).sumAngle - units(j).sumAngle;
            end
            
            if units(i).vmMiu < units(j).vmMiu
                temp_vm = units(i).vmMiu + 360;
                utc.d_dir_vm(1, k) = temp_vm - units(j).vmMiu;
            else
                utc.d_dir_vm(1, k) = units(i).vmMiu - units(j).vmMiu;
            end
            
            %         utc.d_dir_raw(1, k) = units(i).sumAngle - units(j).sumAngle;
            
            %         disp(units(i).vmMiu);
            %         disp(units(j).vmMiu);
            %         utc.d_dir_vm_raw(1, k) = units(i).vmMiu - units(j).vmMiu;
            %         [r, p] = corr(units(i).ave',units(j).ave');
            %         linearCorr(k, :) = [r, p];
            %         fprintf('Correlation coeffient of TT %s unit %s and TT %s unit %s is %.3f, p-value is %.3f; \n', ...
            %             units(i).ttID, units(i).unitID, units(j).ttID, ...
            %             units(j).unitID, r, p);
            
            utc.d_sdi(1, k) = abs(units(i).sdi - units(j).sdi);
            utc.d_mag(1, k) = abs(norm(units(i).polarNorm) - norm(units(j).polarNorm));
        end
    end
end
% subplot(2, 2, 2);
% edges_rose = linspace(0, 2 * pi, 13);
% polarhistogram(theta,'BinEdges', edges_rose + 15 / 180 * pi); hold on
% title('Differences between preferred directions');

%% Plot the distributions of differences of SDIs and normalized magnitudes. 
% edges_hist = linspace(0, 1, 11);
% 
% subplot(2, 2, 3);
% histogram(d_sdi, edges_hist); hold on
% title('Distribution of differences of saccade discrimination indices');
% xlabel('SDI');
% ylabel('Count');
% 
% subplot(2, 2, 4);
% histogram(d_mag, edges_hist); hold on
% title('Distribution of differences of normalized magnitudes');
% xlabel('Normalized magnitude');
% ylabel('Count');

end

%------------- END OF CODE --------------
