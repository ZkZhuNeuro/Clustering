function [utc] = compareUnits_refCell(units, unitAmt)
%This function collect the pairwise data for comparison. 
% 
% Inputs:
%     units - The data I calculated by function 'uniStats_allTT.m'. 
%     unitAmt - The amount of units within a tetrode. 
%     
% Outputs:
%     units.dir_i_x/y - sum vector's direction. 
%     utc.sdi_i_x/y - sdi.
%     utc.mag_i_x/y - normalized magnitudes.
%     utc.dir_vm_i_x/y - von Mises' miu.
%     utc.bw_i_x/y - von Mises' kappa.
%
% See also:
%
% Author: Rosenberg Lab
% email: ari.rosenberg@wisc.edu
% Website: https://neuro.wisc.edu/staff/rosenberg-ari/
% Created: Sept 14 2021, ZKZ
% Editting history:
% 14-Sep-2021, ZKZ: Created the function;
% 18-Sep-2021, ZKZ: Add condition that: 
%     units(i).unitID > 0 && units(j).unitID > 0


%------------- BEGIN CODE --------------
%
% figure;
%% Plot the normalized vector in the same coordinates.
% subplot(2, 2, 1);
% for i = 1:unitAmt
%     polarplot([0, atan2(imag(units(i).polarNorm), ...
%         real(units(i).polarNorm))], [0, norm(units(i).polarNorm)],...
%         'LineWidth', 2); hold on
% end
% title('Nomalized neural vectors');
% hold off

%% Assign data to x/y for comparison.

k = 0;
for i = 1:(unitAmt - 1)
    for j = (i + 1):unitAmt
        if units(i).unitID > 0 && units(j).unitID > 0
            k = k + 1;
            % Put all the angles of vectors into matrix theta.
            
            %         [r, p] = corr(units(i).ave',units(j).ave');
            %         linearCorr(k, :) = [r, p];
            %         fprintf('Correlation coeffient of TT %s unit %s and TT %s unit %s is %.3f, p-value is %.3f; \n', ...
            %             units(i).ttID, units(i).unitID, units(j).ttID, ...
            %             units(j).unitID, r, p);
            utc.dir_i_x(1, k) = units(i).sumAngle;
            utc.dir_i_y(1, k) = units(j).sumAngle;
            
            % Shift the y data for 360 and see which is closer to the
            % identity line.
            d1 =  abs(utc.dir_i_y(1, k) - 360 - utc.dir_i_x(1, k));
            d2 =  abs(utc.dir_i_y(1, k) - utc.dir_i_x(1, k));
            d3 =  abs(utc.dir_i_y(1, k) + 360 - utc.dir_i_x(1, k));
            if min([d1, d2, d3]) == d1
                utc.dir_i_y(1, k) = utc.dir_i_y(1, k) - 360;
            elseif min([d1, d2, d3]) == d3
                utc.dir_i_y(1, k) = utc.dir_i_y(1, k) + 360;
            end
            
            utc.sdi_i_x(1, k) = units(i).sdi;
            utc.sdi_i_y(1, k) = units(j).sdi;
            utc.mag_i_x(1, k) = norm(units(i).polarNorm);
            utc.mag_i_y(1, k) = norm(units(j).polarNorm);
            utc.dir_vm_i_x(1, k) = units(i).vmMiu;
            utc.dir_vm_i_y(1, k) = units(j).vmMiu;
            
            % Shift the y data for 360 degrees and see which is closer to the
            % identity line.
            d4 = abs(utc.dir_vm_i_y(1, k) - 360- utc.dir_vm_i_x(1, k));
            d5 = abs(utc.dir_vm_i_y(1, k) - utc.dir_vm_i_x(1, k));
            d6 = abs(utc.dir_vm_i_y(1, k) + 360 - utc.dir_vm_i_x(1, k));
            if min([d4, d5, d6]) == d4
                utc.dir_vm_i_y(1, k) = utc.dir_vm_i_y(1, k) - 360;
            elseif min([d4, d5, d6]) == d6
                utc.dir_vm_i_y(1, k) = utc.dir_vm_i_y(1, k) + 360;
            end
            
            utc.bw_i_x(1, k) = units(i).vmKappa;
            utc.bw_i_y(1, k) = units(j).vmKappa;
        end
    end
end
% subplot(2, 2, 2);
% 
% scatter(theta_x * 180 / pi, theta_y * 180 / pi); hold on
% plot([0, 360], [0, 360], '-k');
% xticks(linspace(0, 360, 13));
% yticks(linspace(0, 360, 13));
% xlim([0 360]);
% ylim([0 360]);
% title('Differences between preferred directions');

%% Plot the distributions of differences of SDIs and normalized magnitudes. 
% 
% subplot(2, 2, 3);
% scatter(d_sdi_x, d_sdi_y); hold on
% plot([0, 1], [0, 1], '-k');
% title('Distribution of differences of saccade discrimination indices');
% xlim([0 1]);
% ylim([0 1]);
% 
% subplot(2, 2, 4);
% scatter(d_mag_x, d_mag_y); hold on
% plot([0, 1], [0, 1], '-k');
% title('Distribution of differences of normalized magnitudes');
% xlim([0 1]);
% ylim([0 1]);

end
