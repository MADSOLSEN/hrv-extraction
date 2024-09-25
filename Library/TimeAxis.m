% name: ar.classes.Features.All
% Description: time axis object 

classdef TimeAxis
    %TIMEAXIS2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name;
        FS;
        tm;
        discard;
        ectopic;
        target;
        hypnogram;
        Apnea;
        LM;
        Arousal;
        Transition;
        type;
        x;
    end
    
    methods 
        
        % Constructor
        function obj = TimeAxis(name)
            obj.name = name;
        end
        
        function [obj] = addTimeAxis(obj,x,FS,tm,startMargin)
            
            if isempty(tm)
                obj.tm = ((0:length(x)-1)/FS + startMargin);
                obj.FS = (FS);
                obj.discard = zeros(1,length(obj.tm))>0;
                obj.type = 'series';
            else
                obj.tm = (tm(:)');
                obj.FS = (FS);
                obj.discard = zeros(1,length(obj.tm))>0;
                obj.type = 'sample';
            end
            
        end
        
    end
    
end

