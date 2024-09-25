classdef RRTachogram
    %RRTACHOGRAM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Static = true)
        
        function tm = wrapper(x,FS)
            %wrapper: runs all HRV related functions.
            % Will carry out the
            % Inputs: X datavector
            %         FS sampling frequency
            %
            % Output: time axis object containing:
            % type: "tm.getTimeAxisNames" to see what it contains.
            %
            
            % paths
            %addpath(genpath('C:\Users\mads_\OneDrive - Danmarks Tekniske Universitet\Dokumenter\MATLAB\physionet_ecg'));
            
            
            %% Initialize time axis object
            % Organises all time axis that are created in the module
            tm = CollectTimeAxis(x,FS);
            
            
            %% Bandpass filter
            %--------------------------------------------------------------
            RRTachogram.Bandpass_filter(tm);
            
            %% R peak detection
            %--------------------------------------------------------------
            
            % Feature extraction and classification
            RRTachogram.R_peak_detection(tm);
            
            % PostProcessing
            FS_up = 256*3;
            RRTachogram.Template_matching(tm,FS_up);
            
            
            % Storage
            %--------------------------------------------------------------
            objCan = HRVLibrary.ecgCandidate(tm);
            
            
            %% Artefact detection
            %--------------------------------------------------------------
            
            % Features
            objArti = HRVLibrary.CollectFeatures;
            
            objArti.addFeature(HRVLibrary.Features.Morph.Covariance)
            objArti.addFeature(HRVLibrary.Features.Morph.Energy);
            objArti.addFeature(HRVLibrary.Features.Dynamic.RR_pre);
            objArti.addFeature(HRVLibrary.Features.Dynamic.RR_post);
            
            [Artefeatures] = objArti.extract_feature(objCan);
            
            % Classification
            [Arte_IDX] = RRTachogram.ClassifyArtefacts(Artefeatures,objCan);
            
            
            % Update
            %--------------------------------------------------------------
            tm.setDis(Arte_IDX,'qrs up');
            objCan.setDisIDX(tm);
            
            
            %% Ectopic beat detection
            %--------------------------------------------------------------
            
            % Features
            objEcto = HRVLibrary.CollectFeatures;
            
            %Dynamic
            objEcto.addFeature(HRVLibrary.Features.Dynamic.RR_pre);
            objEcto.addFeature(HRVLibrary.Features.Dynamic.RR_post);
            objEcto.addFeature(HRVLibrary.Features.Dynamic.RR_smallLocal);
            objEcto.addFeature(HRVLibrary.Features.Dynamic.RR_bigLocal);
            
            % Morphologic
            objEcto.addFeature(HRVLibrary.Features.Morph.QRS_static_own)
            objEcto.addFeature(HRVLibrary.Features.Morph.PQ_static_own)
            objEcto.addFeature(HRVLibrary.Features.Morph.ST_static_own);
            
            [Ectofeatures] = objEcto.extract_feature(objCan);
            [Ectofeatures] = objEcto.Transform(Ectofeatures);
            [Ectofeatures] = objEcto.PhysNormalize(Ectofeatures);
            
            % Classification
            [Ecto_IDX]  = RRTachogram.ClassifyEctopic(Ectofeatures,objCan,tm);
            
            
            % Update
            %--------------------------------------------------------------
            tm.setDis(Ecto_IDX,'qrs up');
            objCan.setDisIDX(tm);
            
            
            %% Cubic spline interpolation and resampling
            %--------------------------------------------------------------
            
            % Extrac interpolated and resampled RR tachogram
            DCadjust = 1;
            rFS = 4;
            tm.HRVextract(rFS,DCadjust);
            
            % edited RR intervals calculated
            tm.setEditedRRintervals
            
            % trimming
            %tm.rmTimeAxis('EKG up');
            %tm.rmTimeAxis('qrs up');
        end
        
        
    end
    
    methods (Static = true, Access = protected)
        
        %% Bandpass filter
        %------------------------------------------------------------------
        function Bandpass_filter(tm)
            %Bandpass_filter bandpass filters the datavector
            % Equiripple bandpass filter designed in fdatool for signals
            % with sampling frequency 200 Hz. Parameters:
            %
            % Fpass1: 0.8 Hz      Fpass2: 35.0 Hz
            % Fstop1: 0.3 Hz      Fstop2: 35.5 Hz
            % Astop1: 80 dB       Astop2: 80 dB
            %
            % Call: Bandpass_filter(tm)
            %
            % Input: tm         Time axis object. This holds the EKG
            %                   datavetor that is processed. The time axis
            %                   object is updated directly in the function
            
            % Load Bandpass filter
            load('./Resources/DataMatrix/Bandpass.mat')
            a = [1, zeros(1,length(Bandpass)-1)];
            b = Bandpass;
            
            % filtering
            y = filter(b,a,tm.getX('EKG'));
            
            % delay correction
            delay = mean(grpdelay(Bandpass));
            y = [ y(delay+1 : end), zeros(1,delay)];
            
            % Set signal to filtered
            tm.setX(y,'EKG');
            
        end
        
        
        %% R peak detection
        %------------------------------------------------------------------
        
        % Features and classification
        function [qrs_loc, delay, T_high,T_low,feature] = R_peak_detection(tm)
            %R_peak_detection identifies the R peaks in the EKG signal.
            % This function conducts the qrs detection algorithm described
            % by Saadi et al (2015).
            %
            % Call: [qrs_loc, delay, T_high,T_low,feature] = R_peak_detection(tm)
            %
            % Input:  tm        Time axis object. This holds the EKG
            %                   datavetor that is processed. The time axis
            %                   object is updated directly in the function
            %
            % Output: qrs_loc   Locations of the R peaks.
            %         delay     The group delay introduced by the filters.
            %         T_high    High adaptive threshold
            %         T_low     Low adaptive threshold
            %         feature   Feature analysed for R peak detection.
            %                   processed with filters and other processing
            %                   steps
            
            
            % Filtering of ecg data
            [feature, delay] = RRTachogram.ecgfilter(tm);
            
            % Window initialization
            Threshold_window = 2 * tm.getFS('EKG');
            window_exp = Threshold_window - mod(length(feature),Threshold_window);
            feature = [feature, zeros(1,window_exp)];
            windows = reshape(1:length(feature), ...
                [Threshold_window,length(feature)/Threshold_window]); % window indexing
            
            % Parameters
            RRscale = 1.5; % scaling parameter
            Ttheta = 35; % Threshold. # of samples from low to high variability
            alpha = 0.7; % scaling factor omkring 0.7-0.8 estimated
            beta = 0.6; % 70
            Tref = ceil(0.25 * tm.getFS('EKG')); % refractory period
            theta = 0;
            
            % Preallocation doubles
            m = 0; % window number
            case1 = 0; % special case 1
            case2 = 0;
            recentQRS = 0;
            candidateDetected = 0;
            candidatePos = 0;
            candidate = 0;
            refEnd = 0;
            T = nan;
            RRlong = 0;
            Thigh = 0; % higher treshold
            Tlow = 0; % lower treshold
            RRmax = 0;
            QRScounter = 0;
            
            % Preallocation [vectors]
            QRS = zeros(1,length(feature))>0; % QRS complexes
            RRsb = nan(1,8); % Search back RR intervals last 8
            T_high = zeros(1,length(feature));
            T_low = zeros(1,length(feature));
            QRS_loc = zeros(1,length(feature));
            
            for n = 1:(length(feature))
                
                T_high(n) = Thigh;
                T_low(n) = Tlow;
                
                % Threshold updates - will be done no matter what
                
                if rem(n,Threshold_window) == 0  % New window
                    m = m + 1; % Window number
                    
                    % High threshold
                    if m > 8
                        Thigh = median(max(feature(windows(:,m-8:m-1)))) * alpha;
                    end
                    T = Thigh;
                    
                    % Low threshold
                    if m > 2
                        s1 = sum(sum(QRS(windows(:,m-2:m-1)))); % Check it will work
                        if s1 < 1
                            s1 = 1;
                        elseif s1 > 8
                            s1 = 8;
                        end
                        
                        s2 = 10;
                        if theta > Ttheta % High variability
                            s2 = 12;
                        end
                        
                        Tlow_temp = median(mean(feature(windows(:,m-2:m-1)))) * s2/s1;
                        
                        if Tlow_temp > Thigh * beta
                            Tlowprev = Tlow;
                            Tlow = Thigh * beta;
                        else
                            Tlowprev = Tlow;
                            Tlow =  Tlow_temp;
                        end
                    end
                    
                    % RRmax
                    if RRlong ~= 0
                        if theta < Ttheta || theta == Ttheta % Low variability
                            RRtemp = median(RRlong);
                        else % High variability
                            
                            if sum(RRsb>0) > 8
                                RRtemp = min(median(RRsb), median(RRshort));
                            else
                                RRtemp = median(RRshort);
                            end
                            
                        end
                        RRmaxprev = RRmax;
                        RRmax = RRtemp * RRscale;
                        if round(RRmax) < RRn && RRmaxprev > RRn
                            case1 = 1;
                        end
                    end
                    
                end
                
                RRn = n - recentQRS; % Current RR interval
                
                if m > 0
                    if ~isnan(RRmax) && candidateDetected == 0
                        
                        if round(RRmax) > RRn
                            T = Thigh;
                        elseif round(RRmax) == RRn || case1 == 1 % searchback
                            
                            if min(Tlowprev,Tlow) < max(feature(recentQRS+Tref:n))
                                
                                [~,candidatePosR] = max(feature(recentQRS+Tref:n));
                                candidatePos = candidatePosR + recentQRS+Tref - 1;
                                
                                if candidatePos < n - Tref
                                    QRS(candidatePos) = 1;
                                    QRScounter = QRScounter +1;
                                    QRS_loc(QRScounter) = candidatePos;
                                    recentQRSsb = candidatePos;
                                    candidatePos = 0;
                                    case2 = 1;
                                else % If near searchbackEND
                                    refEnd = candidatePos + Tref;
                                    candidate = feature(candidatePos); % new candidate
                                    candidateDetected = 1;
                                end
                            end
                            
                            if case2 > 0
                                RRsb = [RRsb(2:end), recentQRSsb - recentQRS]; % Searchback RR
                                recentQRS = recentQRSsb;
                            end
                            case1 = 0;
                            case2 = 0;
                            
                        else % Low threshold
                            T = Tlow;
                        end
                    end
                end
                
                if (n - recentQRS > Tref) && m > 8 % Check refractory period
                    
                    if candidateDetected == 1
                        
                        if refEnd == n
                            
                            QRS(candidatePos) = 1;
                            QRScounter = QRScounter +1;
                            QRS_loc(QRScounter) = candidatePos;
                            recentQRS = candidatePos;
                            candidateDetected = 0;
                            candidatePos = 0;
                            candidate = 0;
                            
                            if QRScounter > 34
                                
                                RRlong = diff(QRS_loc(QRScounter-34:QRScounter)); %Latest 34 RR intervals
                                epsilon = abs(RRlong - median(RRlong)); % Error
                                epsilonSort = sort(epsilon);
                                
                                % Only when variability is high.
                                theta = mean(epsilonSort(1:end-2)); % Current variable measure
                                RRshort = RRlong(end-7:end);
                            end
                            
                        elseif feature(n) > candidate
                            candidate = feature(n);
                            candidatePos = n;
                            
                        end
                    elseif feature(n) > T
                        candidate = feature(n); % new candidate
                        candidatePos = n; % new candidate position;
                        candidateDetected = 1;
                        refEnd = n + Tref;
                    end
                    
                end
            end
            
            % Trim QRS with uncompleted window
            QRS = QRS(1:end-window_exp);
            T_high = T_high(1:end-window_exp);
            T_low = T_low(1:end-window_exp);
            feature = feature(1:end-window_exp);
            
            % Finds accurate qrs peaks
            qrs_loc = find(QRS);
            qrs_loc = qrs_loc - delay;
            
            while qrs_loc(end) > length(feature)-tm.getFS('EKG')
                qrs_loc(end) = [];
            end
            
            tm.addTimeAxis([],tm.getFS('EKG'),qrs_loc,[],'qrs');
            
        end
        
        function [feature, delay]= ecgfilter(tm)
            %ecgfilter extract the feature used for R peak detection.
            % This function conducts the feature extraction described
            % by Saadi et al (2015).
            %
            % Call: [feature, group_delay]= ecgfilter(tm)
            %
            % Input:  tm        Time axis object. This holds the EKG
            %                   datavetor that is processed. The time axis
            %                   object is updated directly in the function
            %
            % Output: feature   Feature analysed for R peak detection.
            %                   processed with filters and other processing
            %                   steps
            %         delay     The group delay introduced by the filters.
            
            % First bandpass filter
            L = round(10*tm.getFS('EKG')/512);
            M = round(2*L/10);
            block = [1*ones(1,M) zeros(1,L-2*M) (-1)*ones(1,M)];
            hBP1 = [fliplr(block) block];
            
            % Second bandpass filter
            L = round(14*tm.getFS('EKG')/512);
            M = round(2*L/14);
            block = [1*ones(1,M) zeros(1,L-2*M) (-1)*ones(1,M)];
            hBP2 = [fliplr(block) block];
            
            % First lowpass filter
            L = round(16*tm.getFS('EKG')/512);
            hLP1 = ones(1,L)/L;
            
            % Second lowpass filter
            L = round(8*tm.getFS('EKG')/512);
            hLP2 = ones(1,L)/L;
            
            % dffir configuration
            HBP1 = dfilt.dffir(hBP1);
            HBP2 = dfilt.dffir(hBP2);
            HLP = dfilt.dffir(hLP1);
            Hcas = dfilt.cascade(HBP1,HBP2,HLP); % Cascade of FIR filters.
            
            % 1) BP filtering
            yBP = filter(Hcas,tm.getX('EKG')*1000);
            delay1 = grpdelay(Hcas); % First group delay
            
            % 2) Absolute values
            yabs = abs(yBP);
            
            % 3) Smoothening
            HLP2 = dfilt.dffir(hLP2);
            feature = filter(HLP2,yabs);
            
            delay2 = grpdelay(HLP2); % Second group delay
            delay = floor(delay2(1)+delay1(1)); % combined group delay
            
            % plot
            %flag = 0;
            %HRV.plot.Rfilter(flag);
        end
        
        % post processing
        function Template_matching(tm,FS_up)
            %Template_matching corrects the R peak location by template matching
            % Template based cross-corrrelation correction of the R peak
            % location.
            %
            % Call: Template_matching(tm,FS_up)
            %
            % Input:  tm        Time axis object. This holds the EKG
            %                   datavetor that is processed. The time axis
            %                   object is updated directly in the function
            %         FS_up     Sampling frequency of upsampled signal
            
            % upsample signal
            RRTachogram.upsample_ecg(tm,FS_up);
            
            % Find template
            [template ] = RRTachogram.template_ecg(tm);
            
            % Find xcorr
            RRTachogram.fiduacial_adjustment(tm,template);
            
        end
        
        function upsample_ecg(tm,FS_up)
            %upsamle_ecg upsamples the ecg signal
            % Upsamples the ecg datavector from FS to FS_up by spline
            % interpolation.
            %
            % Call: upsample_ecg(tm,FS_up)
            %
            % Input:  tm        Time axis object. This holds the EKG
            %                   datavetor that is processed. The time axis
            %                   object is updated directly in the function
            %         FS_up     Sampling frequency of upsampled signal
            
            % Upsampled parameters
            LEN = length(tm.getX('EKG'));
            LEN_UP = round(FS_up/tm.getFS('EKG')*LEN);
            
            % Upsample
            STEP = 1/(LEN_UP-1)*(LEN-1); % stepsize
            data_up = spline(1:LEN,double(tm.getX('EKG')),1:STEP:LEN);
            
            % qrs location
            qrs_loc_up = round(double(tm.getTM('qrs'))*FS_up/tm.getFS('EKG'));
            
            % add to tm
            tm.addTimeAxis(data_up,FS_up,[],0,'EKG up');
            tm.setX(data_up,'EKG up');
            tm.addTimeAxis([],FS_up,qrs_loc_up,[],'qrs up');
            
        end
        
        function [template] = template_ecg(tm)
            %template_ecg extracts a template from the ecg signal
            % a template of 10 consecutive R peak candidates is found. The
            % 10 R peaks are identified as the 10 consecutive RR peaks with
            % the least variation.
            %
            % Call: [template] = template_ecg(tm)
            %
            % Input:  tm        Time axis object. This holds the EKG
            %                   datavetor that is processed. The time axis
            %                   object is updated directly in the function
            
            % extract vectors from the time axis object
            data_up = tm.getX('EKG up');
            qrs_loc_up = tm.getTM('qrs up');
            RR = diff(qrs_loc_up);
            
            % Generate windows for standard deviation extraction
            WINDOW_LEN = 10; % Window length of 10 consecutive beats
            OVERLAP = 0.90; % We move 1 beat each time
            [windows,~] = RRTachogram.window_create(RR,WINDOW_LEN,OVERLAP);
            
            % Sort by std
            X_w = RR(windows);
            sigma_w = std(X_w);
            [~,stdloc] = sort(sigma_w);
            
            % Consider template...
            k = stdloc(1); % Lowest std
            
            % Find template
            template_loc = round(qrs_loc_up(windows(:,k))); % R peak of each 10 qrs complex.
            localHalfWidth = round(tm.getFS('EKG up') * 0.08);
            segments = zeros(localHalfWidth*2+1,WINDOW_LEN); % preallocation
            for m = 1:WINDOW_LEN; segments(:,m) = data_up(template_loc(m)-localHalfWidth:template_loc(m)+localHalfWidth); end
            
            % extract template from the mean
            template = zscore(mean(segments,2));
            
        end
        
        function fiduacial_adjustment(tm,template)
            %fiducial_adjustment finds the corrected R peak location
            % Template based cross-corrrelation correction of the R peak
            % location.
            %
            % Call: fiduacial_adjustment(tm,template)
            %
            % Input:  tm        Time axis object. This holds the EKG
            %                   datavetor that is processed. The time axis
            %                   object is updated directly in the function
            %         template  template extracted from template_ecg
            %                   function
            
            % RR tachogram
            HR = tm.getFS('qrs up')./gradient(tm.getTM('qrs up')); % beats pr sek
            
            % physiological check
            HR = min(HR, ones(size(HR))*4); % max 4 beats pr sekond
            HR = max(HR, ones(size(HR))*0.25); % min 0.25 beats pr sekond
            
            local_HR = medfilt1(HR,10);
            
            % extract stuff
            data_up = tm.getX('EKG up');
            qrs_loc_up = tm.getTM('qrs up');
            
            %
            local_width_half = round(1./local_HR*0.23*tm.getFS('EKG up'));
            
            % Preallocation
            adj_qrs_loc_up = zeros(size(qrs_loc_up));
            
            for n = 1:length(qrs_loc_up)
                
                % initialize loop
                qrs_loc_old = qrs_loc_up(n);
                qrs_loc_new = qrs_loc_old;
                
                first = 1;
                count = 0;
                
                while qrs_loc_old ~= qrs_loc_new && count < 20 || first
                    
                    qrs_loc_old = qrs_loc_new;
                    
                    % Extract candidate
                    candidate = zscore(data_up(qrs_loc_old-local_width_half(n) : ...
                        qrs_loc_old+local_width_half(n)));
                    rel_qrs_loc = ceil(length(candidate)/2);
                    [~,qrs_loc_templ] = max(template);
                    
                    % Correlation
                    [amp,lag] = xcorr(candidate,template);
                    [~,R_loc] = max(amp);
                    LAG_DIFF = lag(R_loc);
                    rel_diff = ( rel_qrs_loc - qrs_loc_templ - LAG_DIFF);
                    
                    % adjusted qrs location
                    qrs_loc_new = qrs_loc_old -rel_diff;
                    count = count +1;
                    first = 0;
                    
                end
                
                adj_qrs_loc_up(n) = qrs_loc_new;
                
            end
            
            adj_qrs_loc_up = unique(adj_qrs_loc_up);
            
            tm.setTM(adj_qrs_loc_up,'qrs up');
            tm.setDis(zeros(1, length(adj_qrs_loc_up)), 'qrs up');
            
        end
        
        
        %% Artefact detection
        %------------------------------------------------------------------
        
        % Classification
        function [Discard_IDX] = ClassifyArtefacts(features,objCan)
            %ClassifyArtefacts classification of the artefacts
            % Classifies in a rule based manner from the features inputted
            %
            % Call: [Discard_IDX] = ClassifyArtefacts(Morphfeature,objCan)
            %
            % Input:  features  Feature matrix
            %         objCan    Object with stored Template, Heart
            %                   beat candidates and R locations.
            
            % Find areas of artefacts
            Discard_IDX1 = (features(:,1) < 0.35); % Correlation
            Discard_IDX2 = (features(:,2) < 0.2); % Energy
            Discard_IDX3 = (features(:,2) > 15); % Energy
            Discard_IDX4 = (features(:,3) < 0.25); % Pre RR
            Discard_IDX5 = (features(:,3) > 2.5); % Pre RR
            Discard_IDX6 = (features(:,4) < 0.25); % Post RR
            Discard_IDX7 = (features(:,4) > 2.5); % Post RR
            
            Discard_IDX = (Discard_IDX1 + Discard_IDX2 + Discard_IDX3 + ...
                Discard_IDX4 + Discard_IDX5 + Discard_IDX6 + Discard_IDX7)>0;
            
            % Add to the right place
            ObjCanDis = objCan.getDisIDX;
            ObjCanDis(objCan.getFeatNonDisIDX) = Discard_IDX;
            
            Discard_IDX = ObjCanDis;
        end
        
        
        %% Ectopic beat detection
        %------------------------------------------------------------------
        
        % Classification
        function [Discard_IDX] = ClassifyEctopic(features,objCan,tm)
            %ClassifyEctopic classification of the ectopic beats
            % Classifies ectopic beats with a feed forward neural network
            % that has been trained, validated and tested on the MIT
            % Arrhythmia database. Performance is specified in the article.
            % 3 categories are detected. normal, Atrial premature, and
            % ventricular premature beats.
            %
            % Call: [Discard_IDX] = ClassifyEctopic(features,objCan,tm)
            %
            % Input:  features  Feature matrix
            %         objCan    Object with stored Template, Heart
            %                   beat candidates and R locations.
            %         tm        Time axis object. This holds the
            %                   datavetor that are processed. Ectopic
            %                   vector is updated in the function
            
            % Classification
            load('./Resources/nets/EctopicBeats');
            y = net(features');
            classesfull = vec2ind(y);
            ectopic_beats = classesfull>1;
            
            % add to ectopic with type
            % Below does not have any artefacts removed
            EctopicIDX = zeros(size(objCan.getDisIDX));
            EctopicIDX(objCan.getFeatNonDisIDX) = classesfull;
            tm.setEct(EctopicIDX,'qrs up');
            
            % Add to the right place
            % Below does have artefacts removed
            ObjCanDis = objCan.getDisIDX;
            ObjCanDis(objCan.getFeatNonDisIDX) = ectopic_beats;
            
            Discard_IDX = ObjCanDis;
            
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

