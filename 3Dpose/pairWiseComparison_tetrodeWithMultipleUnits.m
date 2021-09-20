function [pwc, mag, kappa, magNs, kappaNs] = pairWiseComparison_tetrodeWithMultipleUnits(tetrodes)
%Get the properties I want to compare from my units.  
% 
% Inputs:
%     tetrodes - Cell data which has filenames of units from the same tetrode.  
%     
% Outputs:
%     pwc - Data ready for pairwise comparison. 
%     mag - Normalized magnitudes of all units. 
%     kappa - von Mises' kappa for all units. 
% 
% See also:
% 
% Author: Rosenberg Lab
% email: ari.rosenberg@wisc.edu
% Website: https://neuro.wisc.edu/staff/rosenberg-ari/
% Created: Sept 14 2021, ZKZ
% Editting history: 
% 14-Sep-2021, ZKZ: Created the function;

 
%------------- BEGIN CODE --------------
%% 'pwc.xxx' is the property I want to compare. 
pwc.d_dir = [];
pwc.d_dir_vm = [];
mag = [];
kappa = [];
magNs = [];
kappaNs = [];
pwc.dir_vm = [];
pwc.dir_sumVec = [];

% Following x/y files are for scatter plots. 
pwc.dir_x = [];
pwc.sdi_x = [];
pwc.mag_x = [];
pwc.dir_vm_x = [];
pwc.bw_x = [];

pwc.dir_y = [];
pwc.sdi_y = [];
pwc.mag_y = [];
pwc.dir_vm_y = [];
pwc.bw_y = [];

%% Collect data from all units from a single tetrode. 
for i = 1:size(tetrodes, 2)
    filenames = tetrodes{i};
    [units, unitAmt] = unitStats_allTT(filenames);
    for j = 1: unitAmt
        pwc.dir_vm = [pwc.dir_vm, units(j).vmMiu];
        pwc.dir_sumVec = [pwc.dir_sumVec, units(j).sumAngle];
    end
    [utcr] = compareUnits_refCell(units, unitAmt);
    
    % Stack my data to the final matrix for comparison. 
    pwc.dir_x = [pwc.dir_x, utcr.dir_i_x];
    pwc.dir_y = [pwc.dir_y, utcr.dir_i_y];
    pwc.sdi_x = [pwc.sdi_x, utcr.sdi_i_x];
    pwc.sdi_y = [pwc.sdi_y, utcr.sdi_i_y];
    pwc.mag_x = [pwc.mag_x, utcr.mag_i_x];
    pwc.mag_y = [pwc.mag_y, utcr.mag_i_y];
    pwc.dir_vm_x  = [pwc.dir_vm_x, utcr.dir_vm_i_x];
    pwc.dir_vm_y  = [pwc.dir_vm_y, utcr.dir_vm_i_y];
    pwc.bw_x = [pwc.bw_x, utcr.bw_i_x];
    pwc.bw_y = [pwc.bw_y, utcr.bw_i_y];
    
    [utc] = compareUnits(units, unitAmt);
    pwc.d_dir = [pwc.d_dir, utc.d_dir];
%     pwc.d_dir_raw = [pwc.d_dir_raw, utc.d_dir_raw];
    pwc.d_dir_vm = [pwc.d_dir_vm, utc.d_dir_vm];
%     pwc.d_dir_vm_raw = [pwc.d_dir_vm_raw, utc.d_dir_vm_raw];
    
    for j = 1:unitAmt
        if units(j).sig == 1
            mag = [mag, norm(units(j).polarNorm)];
            kappa = [kappa, units(j).vmKappa];
        else
            magNs = [magNs, norm(units(j).polarNorm)];
            kappaNs = [kappaNs, units(j).vmKappa];
        end
    end
end

end