% name: HRV.Classes.Features.Dynamic
% Input: objCan containing all heart beats segments. 
% Output: Feature output 
% Description: superclass assigning properties and methods for dynamic features

classdef Dynamic < handle
    %ECTOPIC_FEATURES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract)
        name;
    end
    
    methods (Abstract)
        [feature] = extract_feature(ObjCan)
    end
    
    % Normalization methods
    methods (Static)
        
        % zscore normalization
        function [N_feature] = zNormalize(feature)
            N_feature = zscore(feature);
        end
        
        % Soft normalization
        function [N_feature] = softNormalize(feature)
            Y = quantile(feature,[0.15 0.85]);
            feature_ = bsxfun(@minus,feature,Y(1,:));
            N_feature = bsxfun(@rdivide,feature_,Y(2,:)-Y(1,:));
        end
        
        function Xmirrored = mirrorEdges(X,N,type)
            %mirrorEdges mirrors the edges function X with N samples.
            %
            % Call: Xmirrored = mirrorEdges(X,N,type)
            %
            % Input:  X     The matric or vector to be mirrored
            %         N     Number of samples mirrored
            %         type  Option of mirror type. can be 'start'
            %               'end' or 'middle'.
            
            X = X(:)';
            
            switch type
                case 'start'
                    Xmirrored = [fliplr(X(1:N)), X];
                case 'end'
                    Xmirrored = [ X, fliplr(X(end-N+1:end))];
                case 'middle'
                    Xmirrored = [fliplr(X(1:N)), X, fliplr(X(end-N+1:end))];
            end
            
        end
    end
end

