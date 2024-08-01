%%helper_match_loco_and_sound
% How many Loco-responsive cells are also Sound-responsive?
% Compute and plot pie charts

%% Save summary sheets only from all ED files to make them easier to load

% cd('\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice\New stim analysis\September 19th corrected BW and reliability')
% 
% load('ExtractedData_FM_20230412-010119_newStimAnalysis.mat');
% Summary = ExtractedData.Summary.All;
% save('ExtractedData_Summary_FM.mat', 'Summary');
% 
% load('ExtractedData_RF_20230412-031712_newStimAnalysis.mat');
% Summary = ExtractedData.Summary.All;
% save('ExtractedData_Summary_RF.mat', 'Summary');
% 
% load('ExtractedData_SAM_20230412-031314_newStimAnalysis.mat')
% Summary = ExtractedData.Summary.All;
% save('ExtractedData_Summary_SAM.mat', 'Summary');
% 
% load('ExtractedData_SAMfreq_20230412-041442_newStimAnalysis.mat');
% Summary = ExtractedData.Summary.All;
% save('ExtractedData_Summary_SAMfreq.mat', 'Summary');
% 
% cd('\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice')
% 
% load('ExtractedData_NoiseITI_20230412-001517.mat')
% Summary = ExtractedData.Summary.All;
% save('ExtractedData_Summary_NoiseITI.mat', 'Summary');
% 
% load('ExtractedData_Loco_20230918-173839.mat')
% Summary = ExtractedData.Summary.All;
% save('ExtractedData_Summary_Loco.mat', 'Summary');

%% Load only the summary sheets

cd('\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice\New stim analysis\September 19th corrected BW and reliability\Summary');

load('ExtractedData_Summary_Loco.mat');
Loco = Summary;

load('ExtractedData_Summary_NoiseITI.mat');
Noise = Summary;

load('ExtractedData_Summary_FM.mat');
FM = Summary;

load('ExtractedData_Summary_SAM.mat');
SAM = Summary;

load('ExtractedData_Summary_SAMfreq.mat');
SAMfreq = Summary;

load('ExtractedData_Summary_RF.mat');
RF = Summary;

%% Compute percentage overlap between responses

%Find matching cell in all other sound files
for ii = 1:5
    switch ii
        case 1
            temp = FM;
        case 2
            temp = RF;
        case 3
            temp = SAM;
        case 4
            temp = SAMfreq;
        case 5
            temp = Noise;
    end

    isResp = nan(size(Loco,1),5);

    for i = 1:size(Loco,1)
            
        %Go through every row of Loco and try to find that cell in sound stim
        AP = Loco.AnalysisPath(i);
        cell = Loco.CellNumber(i);
        
        ind = find((strcmp(temp.AnalysisPath, AP) + (temp.CellNumber == cell)) == 2);
        
        if isempty(ind)
            %No matching cell found
            continue;
        elseif length(ind) > 1
            error('maryse')
        elseif length(ind) == 1
            isResp(i, ii) = temp.RF(ind);
        end

    end
    
    switch ii
        case 1
            Loco.FM = isResp(:,ii);
        case 2
            Loco.RecField = isResp(:,ii);
        case 3
            Loco.SAM = isResp(:,ii);
        case 4
            Loco.SAMfreq = isResp(:,ii);
        case 5
            Loco.Noise = isResp(:,ii);
    end
end

%% Sum up

soundmat = [Loco.FM, Loco.RecField, Loco.SAM, Loco.SAMfreq, Loco.Noise];
allnanind = sum(isnan(soundmat),2) == 5; %all rows have sound equivalents

soundmat_nonan = soundmat;
soundmat_nonan(isnan(soundmat)) = 0;

Loco.AnyIsResp = (sum(soundmat_nonan,2)) > 0;

%% Compare

recode = double(Loco.AnyIsResp);
recode(recode == 1) = 2;

merge = Loco.RF + recode;

tabulate(merge);

%% Divide by group

NDNF = merge(strcmp(Loco.Group,'NDNF'));
VIP = merge(strcmp(Loco.Group,'VIP'));
PYR = merge(strcmp(Loco.Group,'PYR'));

N = tabulate(NDNF)
V = tabulate(VIP)
P = tabulate(PYR)

figure;
subplot(1,3,1)
pie(P(:,3), num2str(P(:,1)))
title('PYR')

subplot(1,3,2)
pie(V(:,3), num2str(V(:,1)))
title('VIP')

subplot(1,3,3)
pie(N(:,3), num2str(N(:,1)))
title('NDNF')

sgtitle('0 = none, 1 = loco only, 2 = sound only, 3 = both')