% name: HRV.Classes.Features.Dynamic.RR_bigLocal
% input: objCan containing all heart beats segments. 
% Output: Feature output 
% Description: RR_bigLoc is the median of 300 consequtive RR intervals.

classdef RR_bigLocal < HRVLibrary.Features.Dynamic
    %ECTO_RRLOC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = 'Big Local average';
    end
    
    methods (Static)
        function [feature] = extract_feature(objCan)
            
            % parameter
            bigLocal = 5*60;
            
            % RR tachogram
            RR = diff(objCan.qrs_loc_up);
            
            % process
            RRenlong = HRVLibrary.Features.Dynamic.mirrorEdges(RR,bigLocal,'start');
            RRmedian = medfilt1(RRenlong,bigLocal);
            RRtrim = RRmedian(bigLocal:end);
            
            % extract feature
            feature = (RRtrim(objCan.getFeatNonDisIDX)/objCan.FS_up);
            
        end
    end
    
end

