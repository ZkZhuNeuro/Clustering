function [] = testFunc(x)
% This is a function for testing. 
    p = inputParser; 
    addRequired(p, 'x', @isnumeric);
    bool = parse(p, x);
%     disp(p);
    
end