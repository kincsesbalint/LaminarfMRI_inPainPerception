function BK_checkfile(filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if exist(filename, 'file') == 2
%         disp(['File ', filename, ' exists.']);
    else
        disp(['File ', filename, ' does not exist.']);
    end
end