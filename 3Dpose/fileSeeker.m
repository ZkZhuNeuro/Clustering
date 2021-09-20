function [names] = fileSeeker(Subj, RecordingArea)
%Find the desired files
% Inputs:
%     Subj - monkey names. 
%     RecordingArea - target areas. 
% 
% Outputs:
%     names - filenames of the desired files. 
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
FileName = {};
PathName = {};
PosePathName = {};
for iSubj = 1:length(Subj)
    % Set path
    switch RecordingArea
        case 'V3A'
            initPathName = ['P:\', Subj{iSubj}, '\NeuroRecording\', RecordingArea ,'\FinalPostFit\backup_forJohn\'];
        case 'CIP'
            initPathName = ['P:\', Subj{iSubj}, '\NeuroRecording\', RecordingArea ,'\FinalPostFit\3Dpose_ChoiceSO\'];
        otherwise
            error('ERROR! The recording area is not recognized. Must be V3A or CIP');
    end
    
    [tmpFileName, tmpPathName] = uigetfile({ '*Rate.mat' ; '*.*' }, 'Select Files',initPathName, 'MultiSelect', 'on');
    
    if ischar(tmpFileName) == 1
        tmpFileName = cellstr(tmpFileName);
    end
    
    if ischar(tmpPathName) == 1
        tmpPathName = cellstr(tmpPathName);
    end
    
    PathName = [PathName, repmat(tmpPathName, 1, length(tmpFileName))];
    FileName = [FileName, tmpFileName];
end
% 
% k = 0; 
names = [];
for i = 1:size(FileName, 2)
    names = [names; convertCharsToStrings(FileName{1, i})];
end

end