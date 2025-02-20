function TrialData = simple_check_single_trial_responses(TrialData, StimInfo, Ops, recompute)
%Compute all single trial peak responses using function simple_check_if_responsive
%If peak data already found in TrialData, default is to skip unless recompute == 1

if nargin < 4
    recompute = 0;
end

compute = 1; %Default is to compute

%If peak data already exists and recompute == 0, don't compute
if isfield(TrialData, 'Trial_IsResponsive')
    if recompute == 0
        compute = 0;
        disp('Peak responses already found for TrialData. Skipping...')
    end
end

if compute
    for b = 1:size(TrialData,2) %Blocks
        block_trials = TrialData(b).Trials;
        
        %Preallocate for all variables in PeakData
        %Right now this is hardcoded, maybe one day we can make it automatic
        [IsResponsive, HasPeak, Peak, Peak_3fr, Peak_win, Peak_avg, Peak_Onset, Peak_Latency, Peak_Width, Peak_AUC,...
            HasTrough, Trough, Trough_3fr, Trough_win, Trough_avg, Trough_Onset, Trough_Latency, Trough_Width, Trough_AUC] = deal(nan(size(block_trials,1), size(block_trials,2)));
        Response_Type = strings(size(IsResponsive));

        for c = 1:size(block_trials,1) %Cells
            trials = squeeze(block_trials(c,:,:));
            
            %Perform peak detection on single trials
            [PeakData, ~] = simple_check_if_responsive(trials, StimInfo.nBaselineFrames, StimInfo.fs, Ops.Z_level, Ops.AUC_level, Ops.smTime, 0);
            
            %Store response parameters
            IsResponsive(c,:)   = PeakData.IsResponsive;
            Response_Type(c,:)  = PeakData.ResponseType;
            HasPeak(c,:)        = PeakData.HasPeak;
            Peak(c,:)           = PeakData.Peak;
            Peak_3fr(c,:)       = PeakData.Peak_3fr;
            Peak_win(c,:)       = PeakData.Peak_win;
            Peak_avg(c,:)       = PeakData.Peak_avg;
            Peak_Onset(c,:)     = PeakData.Peak_Onset;
            Peak_Latency(c,:)   = PeakData.Peak_Latency;
            Peak_Width(c,:)     = PeakData.Peak_Width;
            Peak_AUC(c,:)       = PeakData.Peak_AUC;
            HasTrough(c,:)      = PeakData.HasTrough;
            Trough(c,:)         = PeakData.Trough;
            Trough_3fr(c,:)     = PeakData.Trough_3fr;
            Trough_win(c,:)     = PeakData.Trough_win;
            Trough_avg(c,:)     = PeakData.Trough_avg;
            Trough_Onset(c,:)   = PeakData.Trough_Onset;
            Trough_Latency(c,:) = PeakData.Trough_Latency;
            Trough_Width(c,:)   = PeakData.Trough_Width;
            Trough_AUC(c,:)     = PeakData.Trough_AUC;
        end
        
        TrialData(b).Trial_IsResponsive = logical(IsResponsive);
        TrialData(b).Trial_ResponseType = Response_Type;
        TrialData(b).HasPeak            = logical(HasPeak);
        TrialData(b).Peak               = single(Peak);
        TrialData(b).Peak_3fr           = single(Peak_3fr);
        TrialData(b).Peak_win           = single(Peak_win);
        TrialData(b).Peak_avg           = single(Peak_avg);
        TrialData(b).Peak_Onset         = single(Peak_Onset);
        TrialData(b).Peak_Latency       = single(Peak_Latency);
        TrialData(b).Peak_Width         = single(Peak_Width);
        TrialData(b).Peak_AUC           = single(Peak_AUC);
        TrialData(b).HasTrough          = logical(HasTrough);
        TrialData(b).Trough             = single(Trough);
        TrialData(b).Trough_3fr         = single(Trough_3fr);
        TrialData(b).Trough_win         = single(Trough_win);
        TrialData(b).Trough_avg         = single(Trough_avg);
        TrialData(b).Trough_Onset       = single(Trough_Onset);
        TrialData(b).Trough_Latency     = single(Trough_Latency);
        TrialData(b).Trough_Width       = single(Trough_Width);
        TrialData(b).Trough_AUC         = single(Trough_AUC);
    end
end