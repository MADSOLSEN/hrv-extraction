% name: HRV.Classes.Features.Morph.QRS_static_own
% Description: Object containing properties for QRS_static_own feature

classdef QRS_static_own < HRVLibrary.Features.Morphologic
    %MORPH_DTW_QRS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = 'DTW QRS static own';
        P_MEDLEN = 1;
        prefixed_Marking = 0.250;
        postfixed_Marking = 0.380;
    end
    
end

