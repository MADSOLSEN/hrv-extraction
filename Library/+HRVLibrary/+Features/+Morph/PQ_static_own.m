% name: HRV.Classes.Features.Morph.PQ_static_own
% Description: Object containing properties for PQ_static_own feature

classdef PQ_static_own < HRVLibrary.Features.Morphologic
    %MORPH_DTW_QRS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = 'DTW PQ static own';
        P_MEDLEN = 0.4;
        prefixed_Marking = 0.250;
        postfixed_Marking = -0.040; 
    end
end

