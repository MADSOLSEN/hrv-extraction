% name: HRV.Classes.Features.Dynamic.RR_pre
% input: objCan containing all heart beats segments. 
% Output: Feature output 
% Description: RR_pre is the preceding RR interval

classdef RR_pre < HRVLibrary.Features.Dynamic
    %ECTO_RRPREV Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = 'previous RR interval'
    end
    
    methods (Static)
        function [feature] = extract_feature(objCan)
            
            RR = [1 diff(objCan.qrs_loc_up)]; % RR locations
            feature = (RR(objCan.getFeatNonDisIDX)/objCan.FS_up);
            
        end
    end
    
end

