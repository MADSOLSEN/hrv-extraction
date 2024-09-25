% name: HRV.Classes.Features.Morpholigic
% Input: objCan containing all heart beats segments.
% Output: Feature output
% Description: superclass assigning properties and methods for morpholocial features

classdef Morphologic < handle
    %ECTOPIC_MORPH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract)
        name;
        P_MEDLEN;
        prefixed_Marking
        postfixed_Marking
    end
    
    methods
        
        function feature = extract_feature(obj,objCan)
            feature = obj.calcFeature(objCan);
            
        end
        function DTW = calcFeature(obj,objCan)
            
            % Parameters
            N = sum(objCan.getFeatNonDisIDX); % Length of candidates
            w = 20*obj.P_MEDLEN; % window size
            
            % Preallocations
            DTW(1:N) = 0;
            
            % updateMarkers:
            middle_idx = round(objCan.MEDIAN_LEN-(objCan.postfixed_Marking*objCan.FS_down));
            start_idx =  max([1,middle_idx-ceil(obj.prefixed_Marking*objCan.FS_down)]);
            end_idx = min([middle_idx+floor(obj.postfixed_Marking*objCan.FS_down),size(objCan.CandidateMatrix,1)]);
            
            % Template stuff
            t = objCan.Template;
            NormPara = quantile(t,[0.1 1]);
            t = (t - NormPara(1))/(diff(NormPara));
            t = t(start_idx:end_idx);
            
            %
            CanIDX = find(objCan.getFeatNonDisIDX);
            
            % Calc DTW
            for n = 1:N
                
                
                % Candidate
                CanTemp = objCan.CandidateMatrix(:,CanIDX(n));
                
                % normalize candidate
                NormPara = quantile(CanTemp,0.1);
                NormPara(2) = CanTemp(middle_idx);
                CanTemp = (CanTemp - NormPara(1))/(diff(NormPara));
                
                % Extract candidate
                s = CanTemp(start_idx:end_idx)';
                
                dist_diff = size(s,2)/size(t,2);
                
                % DTW
                temp = HRVLibrary.Features.Morphologic.DTW_own(s,t,w)/dist_diff; % Normalize by difference in length.
                temp = temp * (1/obj.P_MEDLEN); % Normalize by length
                DTW(n) = temp; % Normalize by relative energy
                
            end
            
        end
        
    end
    
    methods (Static)
        
        function [D] = DTW_own(s,t,w)
            
            DM = HRVLibrary.Features.Morphologic.gp_dist(s,t);
            
            ns=size(s,2);
            nt=size(t,2);
            
            w=max(w, abs(ns-nt)); % adapt window size
            
            test = ones(size(DM));
            test2 = triu(test,(w+2))>0;
            test3 = tril(test,-(w+2))>0;
            DM((test2+test3)>0) = inf;
            
            if ns >= nt
                [d2,~] = min(DM,[],2);
            else
                [d2,~] = min(DM,[],1);
            end
            
            % output
            D = sum(d2(d2>0));
            
        end
        
        function dist12=gp_dist(x1,x2)
            % Euclid dist between x1 and x2
            % x1  d*N1
            % x2  d*N2
            x12=sum(x1.*x1,1);
            x22=sum(x2.*x2,1);
            %
            dist12=x12'*ones(1,size(x2,2))+ones(size(x1,2),1)*x22 - 2*x1'*x2;
            
            dist12(dist12<0) = 0;
            dist12 = sqrt(dist12);
        end
        
        
    end
    
end

