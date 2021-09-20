function [tetrodes] = fileFilter_tetrodeWithMultipleUnits(names)
%Find the tetrodes with multiple units, and save the filenames in a cell. 
% Inputs:
%     names - all the filenames from the desired animals' target areas. 
%     (outout of 'fileseeker').
% Outputs:
%     tetrodes - a cell that saves unit filenames within the same tetrode. 
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
tetrodes = {};
tt = 1;
while true 
    % 'twmu' is a matrix that temporarily saves the filenames of the files
    % from the same tetrode. 
    twmu = names(1, 1);
    strspl_1 = strsplit(names(1, 1), '_');
    for j = 2:size(names, 1)
        strspl_j = strsplit(names(j, 1), '_');
        if strspl_1(1, 1) == strspl_j(1, 1) && ...
                strspl_1(1, 3) == strspl_j(1, 3) && ...
                strspl_1(1, 5) == strspl_j(1, 5) && ...
                strspl_1(1, 4) == strspl_j(1, 4)
            twmu = [twmu; strsplit(names(j, 1))];
        end
    end 
    
    % save the 'twmu' with more than one filenames. 
    if size(twmu, 1) > 1
        tetrodes{tt} = twmu;
        tt = tt + 1;
    end
    
    % Remove the scanned file names. 
    for j = 1:size(twmu, 1)
        names(names == twmu(j, 1)) = [];
        
    end
    twmu = [];
    
    % If there is less than one file left, break. 
    if size(names, 2) < 1
        break;
    end
end

end