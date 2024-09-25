% name: HRV.Classes.CollectFeatures
% Description: class collecting all featuers

classdef CollectFeatures < handle
    %COLLECTFEATURES Summary of this class goes here
    %   Detailed explanation goes here
    
    % This function is to hold information about each feature. The features
    % themselves are described in a separate class. So this function should
    % be able to:
    % 1) Extract the feature values, based on the data.
    % 2) Extract labels for each feature
    % 3) Transform data
    % 4) Normalize data
    
    properties
        NumFeatures = 0;
        features = {};
    end
    
    methods
        
        function addFeature(obj,feature)
            obj.NumFeatures = obj.NumFeatures + 1;
            obj.features{obj.NumFeatures} = feature;
        end
        
        function featureNames = getFeatureNames(obj)
            featureNames = cell(1,obj.NumFeatures);
            for n = 1:obj.NumFeatures
                featureNames{n} = obj.features{n}.name;
            end
        end
        
        function features = extract_feature(obj,objCan)
            
            % Determin Discard if any:
            FeatSize = sum(objCan.getFeatNonDisIDX);
            
            % Preallocate feature matrix
            features = single(zeros(FeatSize,obj.NumFeatures));
            
            for n = 1:obj.NumFeatures
                thisFeature = obj.features{n}.extract_feature(objCan);
                features(:,n) = thisFeature(:);
            end
            
        end    
        
        function features = Transform(obj,features)
           % Tranformation log modulus
            
           T_IDX = ismember(obj.getFeatureNames,{'DTW PQ static own','DTW QRS static own','DTW ST static own'});
           FeatNum = find(T_IDX); 
           
           for n = FeatNum
               features(:,n) = log(features(:,n));
           end
           
        end
        
        function N_feature = PhysNormalize(obj,features)
            
            try
                load('NormalizationNumbersWSC_AllPatients')
            catch
                error('No NormalizationNumbersWSC_AllPatients')
            end
            
            N_feature = zeros(size(features));
            
            for n = 1:obj.NumFeatures
                FeatIDX = ismember(NormalizationNumbersWSC.FeatureNames,obj.features{n}.name);
                
                N_feature(:,FeatIDX) = (features(:,FeatIDX) - NormalizationNumbersWSC.Quan20min(FeatIDX)) ...
                    /(NormalizationNumbersWSC.Quan80max(FeatIDX) - NormalizationNumbersWSC.Quan20min(FeatIDX))*2 - 1 ;
            end
            
        end
    
    end
    
end
