% name: HRV.Classes.ecgCandidate
% Input: objCan containing all heart beats segments.
% Output: Feature output
% Description: class handling for the ecgCandidate object that hold all ecg
% segments.

classdef ecgCandidate < handle
    %ECG_CANDIDATE Summary of this class goes here
    %   Detailed explanation goes here
    
    % This function is to hold information about heart beat segments (HBS).
    % Each HBS is stored in a matrix, and a template is calculated.
    % be able to:
    % 1) setCandidateMatrix
    % 2) setTemplate.
    % 3) downSample ecg.
    % 4) Local HR ecg.
    
    properties (Constant)
        prefixed_Marking = 0.250;
        postfixed_Marking = 0.380;
    end
    
    properties
        FS_down;
        FS_up;
        qrs_loc_up;
        Markers;
        MEDIAN_LEN;
        CandidateMatrix;
        CandidateR;
        Template;
        TemplateStatus;
        TemplateEnergy;
        DisIDX;
    end
    
    methods
        
        function [obj] = ecgCandidate(tm)
            obj.qrs_loc_up = tm.getTM('qrs up');
            obj.FS_up = tm.getFS('qrs up');
            obj.setCandidateMatrix(tm);
        end
        
        function [obj] = setCandidateMatrix(obj,tm)
            % Parameters
            obj.FS_down = tm.getFS('EKG');
            
            % Find markers
            obj.setMarkers(obj.qrs_loc_up(:),obj.FS_up);
            
            % Set median
            obj.MEDIAN_LEN = round(obj.prefixed_Marking*obj.FS_down) + round(obj.postfixed_Marking*obj.FS_down) + 1;
            
            % Preallocate
            Rcandidates = single(zeros(obj.MEDIAN_LEN,length(obj.qrs_loc_up)));
            
            % Calculate
            for n = 1:length(obj.qrs_loc_up)
                
                candidate = obj.getCandidate(tm.getX('EKG up'),n);
                LEN_CAN = length(candidate);
                FS_SPLINE = obj.MEDIAN_LEN/LEN_CAN*obj.FS_up;
                
                ds_candidate = single(obj.downsample_ecg(obj.FS_up,FS_SPLINE,candidate));
                Rcandidates(:,n) = ds_candidate(1:obj.MEDIAN_LEN);
                
            end
            
            % set candidateMatrix
            obj.CandidateMatrix = Rcandidates;
            
            % set template
            obj.setTemplate(Rcandidates,obj.qrs_loc_up,obj.FS_up);
            
            % Update DisIDX
            obj.setDisIDX(tm);
        end
        
        function setDisIDX(obj,tm)
            obj.DisIDX = tm.getDis('qrs up');
        end
        
        function IDX = getFeatDisIDX(obj)
            IDX = obj.DisIDX;
            %se = strel('line',3,1);
            %im = imdilate(FDisIDX(2:end-1),se);
            %IDX = [1, im, 1]>0;
        end
        
        function IDX = getFeatNonDisIDX(obj)
            IDX = obj.DisIDX == 0;
            %se = strel('line',3,1);
            %im = imdilate(FDisIDX(2:end-1),se);
            %IDX = [0, im==0, 0]>0;
        end
        
        function IDX = getDisIDX(obj)
            IDX = obj.DisIDX;
        end
        
        function IDX = getNonDisIDX(obj)
            IDX = obj.DisIDX == 0;
        end
        
    end
    
    methods (Access = private)
        
        
        function [obj] = setMarkers(obj, qrs_loc,FS)
            % This function will set markers
            
            start_Markers = qrs_loc - ceil(obj.prefixed_Marking*FS);
            end_Markers = qrs_loc + floor(obj.postfixed_Marking*FS);
            mid_Markers = qrs_loc - start_Markers;
            
            obj.Markers = [start_Markers, mid_Markers, end_Markers];
            
        end
        
        
        function candidate = getCandidate(obj,data,n)
            
            % This function will get the Candidate
            candidate = data(obj.Markers(n,1): obj.Markers(n,3));
            
            % Figure out what to do with the candidates
            if isempty(candidate)
                candidate = rand(obj.MEDIAN_LEN,1);
            end
        end
        
        function [obj] = setTemplate(obj,Rcandidates,qrs_loc,FS)
            % This section will find the template
            
            % Local HR
            local_HR = obj.Local_HR_ecg(qrs_loc,FS);
            
            % RR tachogram
            RR = diff(qrs_loc);
            
            % Generate windows for standard deviation extraction
            WINDOW_LEN = 10; % Window length of 10 consecutive beats
            OVERLAP = (1-1/WINDOW_LEN); % We move 1 beat each time
            [windows,~] = HRVLibrary.ecgCandidate.window_create(RR,WINDOW_LEN,OVERLAP);
            
            % Sort by std
            X_w = RR(windows);
            sigma_w = std(X_w);
            sigma_idx = sigma_w < 50;
            
            % index physiological HR:
            HR_idx1 = 0.7< local_HR(6:end-5);
            HR_idx2 = 1.5 >local_HR(6:end-5);
            HR_idx = (HR_idx1+HR_idx2)>1;
            
            Templ_win_idx = (sigma_idx + HR_idx) > 1;
            Templ_idx = unique(windows(:,Templ_win_idx));
            
            % Set template
            if sum(Templ_win_idx) == 0
                disp('Unable to find template');
                obj.TemplateStatus = 2;
                obj.TemplateEnergy = median(Rcandidates,2)' * median(Rcandidates,2);
                obj.Template = zscore(median(Rcandidates,2))';
            else
                obj.TemplateEnergy = median(Rcandidates(:,Templ_idx(1:end-1)),2)' * median(Rcandidates(:,Templ_idx(1:end-1)),2);
                obj.Template = zscore(median(Rcandidates(:,Templ_idx(1:end-1)),2))';
                obj.TemplateStatus = 1;
            end
            
        end
    end
    
    
    methods (Access = private, Static = true)
        
        function [varargout] = downsample_ecg(FS_old,FS_new,varargin)
            % This function will resample the FS_new.
            %
            try
                % Errorcheck
                if nargin > 4
                    error('only two signals allowed')
                elseif nargin < 3
                    error('no signal inputted')
                else
                    data = varargin{1};
                    
                    % Parameters
                    LEN = length(data);
                    LEN_UP = round(FS_new/FS_old*LEN);
                    STEP = 1/(LEN_UP-1)*(LEN-1); % stepsize
                    
                    % Function
                    data_up = spline(1:LEN,data,1:STEP:LEN);
                    varargout{1} = data_up;
                    
                    % Conditional
                    if nargin > 3
                        qrs_loc = varargin{2};
                        qrs_loc_up = qrs_loc*LEN_UP/LEN;
                        varargout{2} = round(qrs_loc_up);
                        
                        annotation = varargin{3};
                        annotation_up = annotation*LEN_UP/LEN;
                        varargout{3} = round(annotation_up);
                    end
                end
            catch
                varargout = varargin;
            end
            
        end
        
        function [local_HR] = Local_HR_ecg(qrs_loc,FS)
            % This function will find the local heart rate of the the ecg
            %
            
            % Parameters
            WINDOW_LEN = 10; % Window length of 10 consecutive beats
            OVERLAP = 0.9; % We move 1 beat each time
            
            % Initialization
            RR = diff(qrs_loc);
            HR = FS./RR;
            [windows,~] = HRVLibrary.ecgCandidate.window_create(RR,WINDOW_LEN,OVERLAP);
            
            % Local HR
            local_HR = median(HR(windows)); % Median heart rate local
            local_HR = [repmat(local_HR(1),1,5), local_HR, repmat(local_HR(end),1,5)];
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

