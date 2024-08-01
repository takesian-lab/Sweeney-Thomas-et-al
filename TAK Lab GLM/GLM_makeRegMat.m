% make regressor mat :
function [Xmat, mask, maskLabel,params] = GLM_makeRegMat(params,Ops,ntfilt,sz,t)
Xmat = [];

    for p =  1:size(params.eventTime,2); % unique regressors (does not include regressors that are neuron traces - those get added later)
        if convertCharsToStrings(params.eventLabel{1,p}) == convertCharsToStrings('Motor Activity')
            if Ops.ShiftLoco
                Stim =  params.eventTime{1,p}';
                XmatTemp = params.eventTime{1,p}';
                  XmatTemp = hankel([zeros(ntfilt(2)-1,1);Stim(1:end-(ntfilt(2)-1))],...
                Stim(end-(ntfilt(2)-1):end)); % temporally shift data by 60 frames
               mask(p) = abs(ntfilt(2));
            else
            %             Stim = params.eventTime{1,p}';
            XmatTemp = params.eventTime{1,p}';
            mask(p) = 1;
            end
            
            maskLabel{p} = 'Motor Activity';
        elseif convertCharsToStrings(params.eventLabel{1,p}) == convertCharsToStrings('Licks')
            if Ops.ShiftLicks
                Stim =  params.eventTime{1,p}';
                XmatTemp = params.eventTime{1,p}';
                XmatTemp = hankel([zeros(ntfilt(2)-1,1);Stim(1:end-(ntfilt(2)-1))],...
                    Stim(end-(ntfilt(2)-1):end)); % temporally shift data by 60 frames
                mask(p) = abs(ntfilt(2));
            else
                XmatTemp = params.eventTime{1,p}';
                mask(p) = 1;
            end
        else
            Stim = zeros(sz(2),1);
            mask(p) = abs(ntfilt(2));
            maskLabel{p} = params.eventLabel{1,p};
            for time=1:length(params.eventTime{p}) % timestamps of events (i.e. sound time)
                sound = params.eventTime{p}(time);
                
                if Ops.ByTrials
                    closest_frame_sound = sound;
                else
                    [~, closest_frame_sound] = min(abs(t(:)-sound));
                end
                
                
                if Ops.useMag && convertCharsToStrings(params.eventLabel{1,p}) == convertCharsToStrings('All V1; Analog')
                    
                    Stim(closest_frame_sound) = Data.([s]).TrialData(bl).Stim_V1(time); % putting stimuli into t
                elseif Ops.useMag && convertCharsToStrings(params.eventLabel{1,p}) == convertCharsToStrings('All V2; Analog')
                    Stim(closest_frame_sound) = Data.([s]).TrialData(bl).Stim_V2(time); % putting stimuli into t
                else
                    Stim(closest_frame_sound) = 1; % putting stimuli into t
                end
            end
            XmatTemp = hankel([zeros(ntfilt(2)-1,1);Stim(1:end-(ntfilt(2)-1))],...
                Stim(end-(ntfilt(2)-1):end)); % temporally shift data by 60 frames
        end
        %-----------------------------------------------------------------------
        Xmat = [Xmat, XmatTemp]; %concat regressors
    end


%-----------------------------------------------------------------------
    if Ops.useAllCells 
        fillAll = find(params.eventCode == 4);
        mask(fillAll) = 1;
        for f = 1:length(fillAll)
        maskLabel{1,fillAll(f)} = 'All Neurons in FOV';
        end
    end
%-----------------------------------------------------------------------
if Ops.Residuals | Ops.CorrResid
           fillResid = find(params.eventCode == 5);
            mask(fillResid) = 1;
            for f = 1:length(fillResid)
            maskLabel{1,fillResid(f)} = 'Residual traces';
            end
end
%-----------------------------------------------------------------------
if Ops.useCorrCell
    fillCorr = find(params.eventCode == 7);
    if Ops.CorrMean
    for i = 1:2
        if i ==1
            maskLabel{fillCorr(i)} = 'Positively correlated cells';
            mask(fillCorr(i)) = 1;
        else
           maskLabel{fillCorr(i)} = 'Negatively correlated cells';
             mask(fillCorr(i)) = 1;
        end
    end
    else
        % todo!!!! Fill in what to do if there are many neurons from the
        % correlation analysis, not just averaged ones!
    end
end



Xmat = [ones(sz(2),1),Xmat]; % constant parameter

end
