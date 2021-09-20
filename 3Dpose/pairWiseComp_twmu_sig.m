function [pwc, mag, kappa] = pairWiseComp_twmu_sig(tetrodes)


%% 
mag = [];
kappa = [];

% 'pwc.xxx' is the property I want to compare. 
pwc.d_dir = [];
pwc.d_dir_vm = [];

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

%%
for i = 1:size(tetrodes, 2)
    
    filenames = tetrodes{i};
    [units, unitAmt, sigAmt] = unitStats_sig(filenames);
    
    % 'sigAmt' is used to tell if this tetrode has more than one
    % significantly tuned units. 
    if sigAmt > 1
        [utcr] = compareUnits_refCell(units, unitAmt);
        
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
        
        [utc] = compareUnits(units, unitAmt, sigAmt);
        pwc.d_dir = [pwc.d_dir, utc.d_dir];
        %     pwc.d_dir_raw = [pwc.d_dir_raw, utc.d_dir_raw];
        pwc.d_dir_vm = [pwc.d_dir_vm, utc.d_dir_vm];
        %     pwc.d_dir_vm_raw = [pwc.d_dir_vm_raw, utc.d_dir_vm_raw];
        
        for j = 1:unitAmt
            if units(j).unitID > 0
                mag = [mag, norm(units(j).polarNorm)];
                kappa = [kappa, units(j).vmKappa];
            end
        end
    end
    
end

end