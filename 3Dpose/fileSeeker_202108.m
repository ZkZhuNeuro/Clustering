function [filenames] = fileSeeker_202108(monkeyName, area, date, tt)
% Find the files that define by the function whichTT.
    % It seems the function dir cannot be directly used for this task. 
    allFiles = ls(pwd);
    Flength = size(allFiles, 1);
    filenames = [];
    for i = 3:Flength
%         fprintf('round# %i \n', i);
        name = convertCharsToStrings(allFiles(i,:));
%         disp(name);
        if contains(name, monkeyName)
            if contains(name, area)
                if contains(name, date)
                    if contains(name, tt)
                        filenames = [filenames; name];
%                         disp(filenames);
                    end
                end
            end
        end
    end
% 	fprintf('desired files: %s \n', filenames);
end