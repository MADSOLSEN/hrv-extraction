% name: HRV.Classes.Features.Morph.Energy
% input: objCan containing all heart beats segments. 
% Output: Feature output 
% Description: The energy between a candidate heart beat and a template

classdef Energy < HRVLibrary.Features.Morphologic
    %MORPH_COVAR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = 'Covariance';
        P_MEDLEN = 1;
        prefixed_Marking = 0.250;
        postfixed_Marking = 0.380;
        
    end
    
    methods
        function [Energy] = calcFeature(~,objCan)
            
            % Preallocations
            Energy = zeros(sum(objCan.getFeatNonDisIDX),1);
            CanIDX = find(objCan.getNonDisIDX);
            E1 = objCan.TemplateEnergy;
            
            for n = 1:length(Energy)    
                E2 = objCan.CandidateMatrix(:,CanIDX(n))' * objCan.CandidateMatrix(:,CanIDX(n));
                Energy(n) = E2/E1;
            end

        end
    end
    
end




