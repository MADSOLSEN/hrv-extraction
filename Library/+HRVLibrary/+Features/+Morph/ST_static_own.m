% name: HRV.Classes.Features.Morph.ST_static_own
% Description: Object containing properties for ST_static_own feature

classdef ST_static_own < HRVLibrary.Features.Morphologic
    %MORPH_DTW_QRS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = 'DTW ST static own';
        P_MEDLEN = 0.6;
        prefixed_Marking = -0.040;
        postfixed_Marking = 0.38;
    end
    
end

