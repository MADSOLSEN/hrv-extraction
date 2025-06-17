% name: ar.classes.Features.All
% Destription: collects all features in a matrix for a file.
% operations for the time axis object

classdef CollectTimeAxis < handle
    %TIMEAXIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N_TimeAxis = 0;
        timeAxis ={};
        tm_discard = [];
    end
    
    methods
        function obj = CollectTimeAxis(x,FS)
            
            obj.addTimeAxis(x,FS,[],0,'EKG');
            obj.setX(x,'EKG');
            
        end
        
        function addTimeAxis(obj,varargin)
            % ADDTIMEAXIS adds another time axis.
            %   obj.ADDTIMEAXIS(x,FS,[],startMargin,name) will create a
            %   timeAxis with length(x) indexes with sample frequency FS
            %   that starts in startMargin and has the name name.
            %
            %   obj.ADDTIMEAXIS([],FS,tm,[],name) will create a timeAxis
            %   with tm with the name name. FS is the sample frequency from
            %   the where the time series was extracted.
            %
            %   Copyright May 2016, All Rights Reserved.
            %   Mads Olsen <mads_olsen123@hotmail.com>
            
            if ~ismember(varargin{5},obj.getTimeAxisNames)
                obj.N_TimeAxis = obj.N_TimeAxis + 1;
                obj.timeAxis{obj.N_TimeAxis} = TimeAxis(varargin{5}); % creating the object
                obj.timeAxis{obj.N_TimeAxis} = obj.timeAxis{obj.N_TimeAxis}.addTimeAxis(varargin{1:4}); % create timeAxis
                
                % Update discard
                %obj.oneDiscard;
            end
        end
        
        function rmTimeAxis(obj,name)
            % RMTIMEAXIS will remove the most recent timeAxis with 'name'.
            %
            %   Copyright May 2016, All Rights Reserved.
            %   Mads Olsen <mads_olsen123@hotmail.com>
            
            tm_IDX = find(ismember(obj.getTimeAxisNames,name));
            obj.timeAxis(tm_IDX(end)) = [];
            obj.N_TimeAxis = obj.N_TimeAxis - 1;
        end
        
        function timeAxisNames = getTimeAxisNames(obj)
            % GETTIMEAXISNAMES will get the names for each timeAxis
            %
            %   Copyright May 2016, All Rights Reserved.
            %   Mads Olsen <mads_olsen123@hotmail.com>
            
            timeAxisNames = cell(obj.N_TimeAxis,1);
            for n = 1:obj.N_TimeAxis
                timeAxisNames{n} = obj.timeAxis{n}.name;
            end
            
        end
        
        function HRVextract(obj,rFS,DCadjust)
            
            if sum(ismember(obj.getTimeAxisNames,'qrs up'))
                
                obj.setRR_tachogram;
                obj.setHRV(rFS,DCadjust);
                
            else
                disp('requires qrs loc signal');
            end
            
        end
        
        
        function spreadDis(obj,from,to)
            % Only the large areas should be removed.
            
            % Extract tm IDX
            if sum(ismember({from,to},obj.getTimeAxisNames)) == 2;
                
                TMFrom = obj.getTM(from);
                TMTo = obj.getTM(to);
                disFrom = obj.getDis(from);
                
                % Morphologic Operations
                se1 = strel('line',13,1); % Over 1.5 second on each side
                se2 = strel('line',41,1); % Over five seconds on each side
                im = imclose(disFrom,se1);
                im2 = imopen(im,se2);
                labels = bwlabel(im2);
                
                toDis = zeros(1,length(obj.getTM(to)));
                
                for n = 1:max(labels)
                    
                    segment = TMFrom(labels == n);
                    [~, locstart] = min(abs(segment(1) - TMTo));
                    [~, locstop] = min(abs(segment(end) - TMTo));
                    toDis(locstart:locstop) = 1;
                    
                end
                obj.setDis(toDis,to);
            else
                disp(['Missing one of the timeAxis: ',from,' OR ',to]);
            end
        end
        %
        
        function alignHypnogram(obj,hypObj,name)
            
            % Hypnogram file
            hyp_file = hypObj.getHypnogramSplit(30/(1/obj.getFS(name)));
            hyp_file_tm = (1:length(hyp_file))/obj.getFS(name);
            
            % timeAxis
            hyp_tm = obj.getTM(name);
            
            [~,startLoc] = min(abs(hyp_tm(1) - hyp_file_tm));
            [~,stopLoc] = min(abs(hyp_tm(end) - hyp_file_tm));
            
            hyp = hyp_file(startLoc:stopLoc);
            
            % Set hypnogram
            obj.setHyp(hyp,name);
            
        end
        
        function [ MatrixIDX ] = getTargetMatrixIDX( obj , windowSize , name )
            %AROUSALMATRIX will output the indexes of the Events to extraction for a
            % signal
            %
            % Call: [ MatrixIDX ] = getTargetMatrixIDX(windowSize, name )
            
            % getTarget
            targets = obj.getTarget(name);
            labels = bwlabel(targets);
            
            N_preSamples = round(windowSize(1) * obj.getFS(name));
            N_postSamples = round(windowSize(2) * obj.getFS(name));
            
            % Parameters
            MatrixIDX = zeros(N_preSamples + N_postSamples + 1 , max(labels));
            dummy = 1:length(targets);
            
            for n = 1:max(labels)
                segment = dummy(labels == n);
                startIDX = segment(1) - N_preSamples;
                stopIDX = segment(1) + N_postSamples;
                if stopIDX <= max(dummy)
                    MatrixIDX(:,n) = startIDX : stopIDX;
                end
            end
            MatrixIDX(:,(sum(MatrixIDX,1) == 0)) = [];
            
        end
        
        function setEvent2Target(obj,Events,name,widening,refPeriod)
            % This function will set the target of the time series
            
            target = zeros(size(obj.getTM(name)));
            
            for n = 1:length(Events)
                [Amp1,startIDX]= min(abs((Events(n).time_second - widening(1)) - obj.getTM(name)));
                [Amp2,stopIDX]= min(abs((Events(n).time_second_stop + widening(2)) - obj.getTM(name)));
                if (Amp1 < 0.5) && (Amp2 < 0.5)
                    target(startIDX:stopIDX) = 1;
                    target(max([1, stopIDX+1]):min([length(target), stopIDX + round(refPeriod*obj.getFS(name))])) = nan;
                end
            end
            
            obj.setTarget(target,name);
            
        end
        
        function EventTM = getEventTM(obj,Events,widening,refPeriod)
            % This function will set the target of the time series
            
            name = 'Feature';
            EventTM = zeros(size(obj.getTM(name)));
            
            % If arousals
            if ismember(Events(1).type,{'LM','Respiratory Event','Spontaneous'}) % arousals
                for n = 1:length(Events)
                    EventFactor = find(ismember({'LM','Respiratory Event','Spontaneous'},Events(n).type));
                    [Amp1,startIDX]= min(abs((Events(n).time_second - widening(1)) - obj.getTM(name)));
                    [Amp2,stopIDX]= min(abs((Events(n).time_second_stop + widening(2)) - obj.getTM(name)));
                    if (Amp1 < 0.5) && (Amp2 < 0.5) % makes sure only arousals within the valid HRV signal is considered.
                        EventTM(startIDX:stopIDX) = 1 * EventFactor;
                        EventTM(max([1, stopIDX+1]):min([length(EventTM), stopIDX + round(refPeriod*obj.getFS(name))])) = nan;
                    end
                end
            else
                % in other than arousal
                for n = 1:length(Events)
                    [Amp1,startIDX]= min(abs((Events(n).time_second - widening(1)) - obj.getTM(name)));
                    [Amp2,stopIDX]= min(abs((Events(n).time_second_stop + widening(2)) - obj.getTM(name)));
                    if (Amp1 < 0.5) && (Amp2 < 0.5) % makes sure only arousals within the valid HRV signal is considered.
                        EventTM(startIDX:stopIDX) = 1;
                        EventTM(max([1, stopIDX+1]):min([length(EventTM), stopIDX + round(refPeriod*obj.getFS(name))])) = nan;
                    end
                end
            end
            
        end
        
        function setEditedRRintervals(obj)
            % Segments with artefact or ectopic beats are replaced by
            % interpolated RR intervals based on the surrounding RR
            % intervals. This is done nomatter the amount of missing RR
            % intervals, since too long segments will be removed later
            % anyway.
            
            RRdis = obj.getDis('RR');
            RR = obj.getX('RR');
            RRtm = obj.getTM('RR');
            
            labels = bwlabel(RRdis);
            dummy = 1:length(RR);
            
            for n = fliplr(1:max(labels))
                seg = dummy(n == labels);
                
                if seg(1) > 1 && seg(end)<length(dummy) % Check not start/end
                    
                    % beat before and after the segment of interest
                    loc1 = seg(1)-1;
                    loc2 = seg(end)+1;
                    
                    % Estimate local median
                    DisLN = RRtm(loc2) - RRtm(loc1);
                    medIDX = (max([1, loc1-10]):min([loc1+10, length(RR)]));
                    RRmed = median(RR(medIDX(RRdis(medIDX)==0)));
                    
                    % Estimate number of missing beats
                    NBeats = round(DisLN/RRmed) - 1; % 1 beat is median...
                    
                    % Linear interpolation of missing beats
                    RRNew = linspace(RR(loc1),RR(loc2),NBeats+2);
                    RRtmNew = RRtm(loc1) + cumsum(RRNew(2:end-1));
                    
                    % Insert new beats in time series vectors
                    RRtm = [RRtm(1:loc1),RRtmNew,RRtm(loc2:end)];
                    RR = [RR(1:loc1),RRNew(2:end-1),RR(loc2:end)];
                    RRdis = [RRdis(1:loc1),zeros(1,NBeats),RRdis(loc2:end)];
                    
                end
            end
            
            % Save RR tachogram
            obj.addTimeAxis([],1,RRtm,[],'RRnew');
            obj.setDis(RRdis>0,'RRnew');
            obj.setX(RR,'RRnew');
            
        end
        
        
        % get functions
        function tm = getTM(obj,name,varargin)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            tm = obj.timeAxis{tm_IDX}.tm;
            if nargin > 2
                tm = obj.timeAxis{dis_IDX}.tm;
            end
        end
        
        function hyp = getHyp(obj,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            hyp = obj.timeAxis{tm_IDX}.hypnogram;
        end
        
        function dis = getDis(obj,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            dis = obj.timeAxis{tm_IDX}.discard;
        end
        
        function FS = getFS(obj,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            FS = obj.timeAxis{tm_IDX}.FS;
        end
        
        function x = getX(obj,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            x = obj.timeAxis{tm_IDX}.x;
        end
        
        function ect = getEct(obj,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            ect = obj.timeAxis{tm_IDX}.ectopic;
        end
        
        function target = getTarget(obj,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            target = obj.timeAxis{tm_IDX}.target;
        end
        
        function target = getTransition(obj,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            target = obj.timeAxis{tm_IDX}.Transition;
        end
        
        function target = getApnea(obj,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            target = obj.timeAxis{tm_IDX}.Apnea;
        end
        
        function target = getLM(obj,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            target = obj.timeAxis{tm_IDX}.LM;
        end
        
        function target = getArousal(obj,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            target = obj.timeAxis{tm_IDX}.Arousal;
        end
        
        
        % set functions
        function setTarget(obj,x,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            obj.timeAxis{tm_IDX}.target = (x(:)');
        end
        
        function setApnea(obj,x,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            obj.timeAxis{tm_IDX}.Apnea = (x(:)');
        end
        
        function setLM(obj,x,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            obj.timeAxis{tm_IDX}.LM = (x(:)');
        end
        
        function setTransition(obj,x,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            obj.timeAxis{tm_IDX}.Transition = (x(:)');
        end
        
        function setArousal(obj,x,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            obj.timeAxis{tm_IDX}.Arousal = (x(:)');
        end
        
        function setTM(obj,x,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            obj.timeAxis{tm_IDX}.tm = (x(:)');
        end
        
        function setHyp(obj,x,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            obj.timeAxis{tm_IDX}.hypnogram = (x(:)');
        end
        
        function setDis(obj,x,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            obj.timeAxis{tm_IDX}.discard = (x(:)');
        end
        
        function setEct(obj,x,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            obj.timeAxis{tm_IDX}.ectopic = (x(:)');
        end
        
        function setFS(obj,x,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            obj.timeAxis{tm_IDX}.FS = (x(:)');
        end
        
        function setX(obj,x,name)
            tm_IDX = ismember(obj.getTimeAxisNames,name);
            obj.timeAxis{tm_IDX}.x = (x(:)');
        end
    end
    
    methods (Access = protected)
        
        function [RR] = setRR_tachogram(obj)
            % This function will
            
            % RR tachogram
            Rloc = obj.getTM('qrs up')/obj.getFS('qrs up');
            RR = diff(Rloc);
            RRtm = cumsum(RR) + Rloc(1) - RR(1)/2;
            
            % Discard IDX
            RDiscard = obj.getDis('qrs up');
            RRDiscard = (RDiscard(1:end-1) + RDiscard(2:end))>0;
            
            % Additional RR editing
            RRclass = obj.RRectopicDetection(RR);
            RRDiscard = (RRDiscard + RRclass)> 0;
            
            % Save RR tachogram
            obj.addTimeAxis([],1,RRtm,[],'RR');
            obj.setDis(RRDiscard>0,'RR');
            obj.setX(RR,'RR');
            
        end
        
        function setHRV(obj,rFS,adjust)
            % This will perform HRV extraction
            
            % Extract HRV
            RR = obj.getX('RR');
            RR_tm = round(obj.getTM('RR')*obj.getFS('qrs up'));
            RR_spline = RR_tm(1):RR_tm(end);
            startMargin = RR_tm(1)/obj.getFS('qrs up');
            
            RRdis = obj.getDis('RR');
            RR_tm([0 RRdis(2:end-1) 0]>0) = [];
            RR([0 RRdis(2:end-1) 0]>0) = [];
            
            RRex_mask = obj.getDis('RR');
            RRex_tm = round(obj.getTM('RR')*obj.getFS('qrs up'));
            [value, number] = hist(RRex_tm,unique(RRex_tm));
            IDX = ismember(RRex_tm,number(value==2));
            RRex_tm(IDX) = [];
            RRex_mask(IDX) = [];
            
            % Cubic spline interpolation
            %sRR = spline(RR_tm,RR,RR_spline);
            sRR = pchip(RR_tm, RR, RR_spline);
            sRR_ex = spline(double(RRex_tm),double(RRex_mask),double(RR_spline));
            sRR_ex = sRR_ex>0.5;
            
            % Resample
            HRV = sRR(1:obj.getFS('qrs up')/rFS:end);
            HRVex = sRR_ex(1:obj.getFS('qrs up')/rFS:end);
            
            % add time series
            obj.addTimeAxis(HRV,rFS,[],startMargin,'HRV');
            obj.setDis(HRVex,'HRV');
            
            % Edit HRV with DC in discarded areas:
            if adjust
                HRV = obj.insertDC(HRV);
            end
            
            % add signal
            obj.setX(HRV,'HRV');
            
        end

        function setEDR(obj,rFS,adjust)
            
            % Extract Rpeak from ECG
            Rloc = obj.getTM('qrs up');
            ECG_up = obj.getX('EKG up');
            R_peak = ECG_up(Rloc);

            % extract Phase Space Reconstruction
            interval = round(obj.getFS('qrs up') * 0.04);
            ECG_up = obj.getX('EKG up');
            R_area = zeros(size(Rloc));
            for n = 1:length(Rloc)
                QRS_interval_temp = ECG_up(Rloc(n) - interval:Rloc(n) + interval);
                psr = phasespace(QRS_interval_temp);
                R_area(n) = polyarea(psr(:,1), psr(:,2));
            end
            
            % time axis
            R_tm = obj.getTM('qrs up');
            R_spline = R_tm(1):R_tm(end);
            
            % Discard axis
            R_peak(obj.getDis('qrs up')>0) = [];
            R_area(obj.getDis('qrs up')>0) = [];
            R_tm(obj.getDis('qrs up')>0) = [];
            
            % Interpolation
            sRR_EDR_peak = pchip(R_tm, R_peak, R_spline); % EDR_peak
            sRR_EDR_psr = pchip(R_tm, R_area, R_spline); % EDR_psr
            
            % Resample
            EDR_peak = sRR_EDR_peak(1:obj.getFS('qrs up')/rFS:end);
            EDR_psr = sRR_EDR_psr(1:obj.getFS('qrs up')/rFS:end);

            % add time series
            obj.addTimeAxis(EDR_peak,rFS,[],startMargin,'EDR_peak');
            obj.setDis(HRVex,'EDR_peak');
            obj.addTimeAxis(EDR_psr,rFS,[],startMargin,'EDR_psr');
            obj.setDis(HRVex,'EDR_psr');

            % Edit HRV with DC in discarded areas:
            if adjust
                EDR_peak = obj.insertDC(EDR_peak);
                EDR_psr = obj.insertDC(EDR_psr);
            end

            % add signal
            obj.setX(EDR_peak,'EDR_peak');
            obj.setX(EDR_psr,'EDR_psr');
        end
        
        function HRVa = insertDC(obj,HRV)
            % This function will perform interpolation points
            
            IDX = obj.getDis('HRV');
            %se1 = strel('line',13,1); % Over 1.5 second on each side
            se2 = strel('line',41,1); % Over five seconds on each side
            %im = imclose(IDX,se1);
            im2 = imopen(IDX,se2);
            labels = bwlabel(im2);
            HRVa = HRV;
            
            for n = 1:max(labels)
                segment = HRV(labels == n);
                HRVa(labels == n) = linspace(segment(1),segment(end),length(segment));
            end
        end
    end
    
    methods (Static)
        
        function RRclass = RRectopicDetection(RR)
            % This function will carry out RRtm < 0.5 or RR > 1.5 validation
            
            RRclass = zeros(size(RR));
            
            % estimate mean of 21 windows with 0.5 overlap
            windows = CollectTimeAxis.window_create(RR,21,0.8);
            RRmu = median(RR(windows),1);
            
            LI = (bsxfun(@minus, RR(windows),RRmu/2) < 0);
            HI = (bsxfun(@minus,RR(windows),RRmu*1.5) > 0);
            
            LI_cand = windows(LI);
            HI_cand = windows(HI);
            
            [LI_amount,LI_cand] = hist(LI_cand,unique(LI_cand));
            [HI_amount,HI_cand] = hist(HI_cand,unique(HI_cand));
            
            LI_dis = LI_cand(LI_amount > 2); % More than one occurence
            HI_dis = HI_cand(HI_amount > 2);
            
            RRclass(LI_dis) = 1;
            RRclass(HI_dis) = 1;
            
        end
        
        function [windows, W_STEP] = window_create(signal, window_length, overlap)
            % This function will create windows for indexing.
            %
            % Call: [windows, signal, sig_trunc, w_step] = window_create(signal, window_length, overlap)
            %
            % Input:    signal          The signal to be windowed
            %           window_length   The length of the window segments. Must be
            %                           shorter than the length of the signal.
            %           overlap         The amount of overlap. Must be between 0 and
            %                           0.99
            %
            % Output:   windows         A matrix of the window indexes. Row length will
            %                           be equal to the window_length and coloumns
            %                           continue until signal has ended.
            %           signal          A truncated version of the signal, this is to
            %                           match the window_length in the last coloumn.
            %           w_step          The step size for between each consecutive
            %                           windows.
            %
            
            % Errorcheck on overlap
            while overlap > 0.999 || overlap < 0 % Overlap check
                overlap = input('\nInvalid overlap, must be between 0 and 0.99.\nOverlap: ');
            end
            
            % Parameters
            SIG_LEN = length(signal);
            W_STEP = round(window_length * (1-overlap)); % Step for each window
            N_WINDOWS = floor((SIG_LEN-window_length)/(W_STEP)) + 1; % number of windows
            
            % Preallocation
            windows = zeros(window_length,N_WINDOWS); % preallocate windows
            
            % Generate the window matrix to index the signal
            for n = 1:N_WINDOWS
                windows(:,n) = (1+(n-1)*W_STEP: window_length+(n-1)*W_STEP);
            end
        end
        
    end
    
end

