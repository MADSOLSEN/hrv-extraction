% name: HRV.Classes.Features.Dynamic.RR_post
% input: objCan containing all heart beats segments. 
% Output: Feature output 
% Description: RR_post is the posterior RR interval 

classdef RR_post < HRVLibrary.Features.Dynamic
    %ECTO_RRPOST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = 'post RR interval'
    end
    
    methods (Static)
        function [feature] = extract_feature(objCan)
            
            RR = [diff(objCan.qrs_loc_up) 1];
            feature = (RR(objCan.getFeatNonDisIDX)/objCan.FS_up);
            
        end
    end
    
end

