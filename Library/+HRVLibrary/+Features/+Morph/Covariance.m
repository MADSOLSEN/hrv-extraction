% name: HRV.Classes.Features.Morph.Covariance
% input: objCan containing all heart beats segments. 
% Output: Feature output 
% Description: The covariance between a candidate heart beat and a template

classdef Covariance < HRVLibrary.Features.Morphologic
    %MORPH_COVAR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = 'Covariance';
        P_MEDLEN = 1;
        prefixed_Marking = 0.250;
        postfixed_Marking = 0.380;
        
    end

    methods
        function [Morph_covar] = calcFeature(~,objCan)
            
            % Preallocations
            Morph_covar = zeros(sum(objCan.getFeatNonDisIDX),1);
            CanIDX = find(objCan.getNonDisIDX);
            
            for n = 1:length(Morph_covar)
                temp = cov(zscore([objCan.Template',objCan.CandidateMatrix(:,CanIDX(n))]));
                Morph_covar(n) = temp(2);
                %plot(objCan.Template); hold on;
                %plot(zscore(objCan.CandidateMatrix(:,CanIDX(n)))); hold off;
            end
        end
    end
    
end

