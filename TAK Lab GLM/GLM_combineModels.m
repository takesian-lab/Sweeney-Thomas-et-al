% wrapper script for GLM_ExDat
% output of GML_ExDat is RidgeData: load the RD from each stim protocol,
% identify matching cells, and average the CC between holdout and predicted
% data.

StimTypes = {'FM', 'RF','SAMnoise','SAMfreq'};
RidgePath = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\GLM\Example Features Nov2023\AllResidualsInFOV_allStim';

cd(RidgePath)
for s = 1:length(StimTypes)
    ridgename = strcat('RidgeRegression_', string(StimTypes{s}));
    load(ridgename)
    Data.([StimTypes{s}]) = RidgeData;
    clear RidgeData
end

%% create a list of all cells
CellList = table;
fn = fieldnames(Data);
for r = 1:length(fn)
    sn = fn{r};
    [celldata] = cellList_from_ridge(Data.([sn]),sn);
    CellList = [CellList; celldata];
end
clear celldata

%% work in progress... plot and get cc using mask codes. 
% working on creating a function that will get a cc from hold out and
% predicted data/ will also generate a graph. 
groups = unique(CellList.Group);
for G = 1: length(groups)
    gIX = find(strcmp(groups{G},CellList.Group));
Allmatches = unique(CellList.cellMatch);
matches = Allmatches(gIX);
Mean_CC_All = [];
GroupList = [];
for m = 1:length(matches)
    cellOI = matches(m); % 
    cellRows = find(CellList.cellMatch == cellOI); % rows that correspond to the same cell
    figure;
    drawnow
    t = tiledlayout('flow')
    for c = 1:length(cellRows) % get all the data specific for this cell for this stim type
        tblock = CellList.BlockName{cellRows(c)}; % blockname for CellOI for the matching row
        tstim = CellList.StimType(cellRows(c)); % stim type
        tmouse = CellList.MouseID(cellRows(c));
        tmatch = string(CellList.cellMatch(cellRows(c)));
        tname = strcat(tmouse,tstim,tmatch); % for figure
        tdata = RD.([tstim]); % Ridge data for a given stim type (i.e. FM)
        brow = find(strcmp(tblock,tdata.BlockName)); % row for the block that matches tblock
        T = tdata.TestTrace{brow}(CellList.cellIDX(cellRows(c)),:); % Ca data to test
        Ws =  RD.([tstim]).Filter{brow}; % filter from validation
        const = Ws(1,CellList.cellIDX(cellRows(c)));
        maskCode = [2];
        maskidx = find(ismember(RD.([tstim]).params(brow).eventCode,maskCode));
        MaskTraces = squeeze(RD.([tstim]).MaskTraces{brow}(CellList.cellIDX(cellRows(c)),:,:));
        [p_sound,tCC] = cc_by_mask('hasMaskTraces',1,maskidx,const,MaskTraces,T);
        sound_ccByMasktemp(c) = tCC;
        [p_all,aCC] = cc_by_mask(1:length(RD.([tstim]).params(brow).eventCode),const,MaskTraces,T);
    
        
%        ccByMask = [];
       nexttile
       plot(T); hold on
       plot(p_all)
       plot(p_sound)
       title(tname);
       title(t,tmouse);
    end
%            mean_cc = mean(ccByMasktemp,1);
    end
    
 
    clear ccByMasktemp
    %-----

 % remove nan data:
%  if ~isnan(mean_cc)
% 
%   GroupList = [GroupList; CellList.Group(cellRows(1))];
%   Mean_CC_All = [Mean_CC_All;mean_cc];
% 
%  end



end



% from this point on, it's old code that I havent touched. use as a
% reference:
%% Loop through matching cells and find CC for each leave one out data
matches = unique(CellList.cellMatch);
Mean_CC_All = [];
GroupList = [];
for m = 1:length(matches)
    cellOI = matches(m); % 
    cellRows = find(CellList.cellMatch == cellOI); % rows that correspond to the same cell
    for c = 1:length(cellRows)
        tblock = CellList.BlockName{cellRows(c)}; % blockname for CellOI for the matching row
        tstim = CellList.StimType(cellRows(c)); % stim type
        tdata = RD.([tstim]); % Ridge data for a given stim type (i.e. FM)
        brow = find(strcmp(tblock,tdata.BlockName)); % row for the block that matches tblock
        uniqueMaskCode = unique(tdata.params(brow).eventCode); % mask codes
        T = tdata.TestTrace{brow}(CellList.cellIDX(cellRows(c)),:); % Ca data to test
        Ws =  RD.([tstim]).Filter{brow}; % filter from validation
        mask =  RD.([tstim]).mask{brow};
        maskLabel = RD.([tstim]).maskLabel{brow};
       %-----
        ccByMask = [];
        for j = 1:length(uniqueMaskCode)+1%(mask,2)
            maskset = 1:size(mask,2); 
            if j<length(uniqueMaskCode)+1
            codeCheck = uniqueMaskCode(j);%RidgeData.params(bl).eventCode(j); % which masks should be grouped together. If codes match
            
           matchCodes = find(codeCheck == RD.([tstim]).params(brow).eventCode);
           maskset(matchCodes) = []; % leave out this code (aka regressor)
            end
            % loop through cells in FOV :

%                 T = testTrace(CellList.cellIDX(cellRows(c)),:); %
                const = Ws(1,CellList.cellIDX(cellRows(c)));
                leaveone = 0;
                for k = 1:length(maskset) % masks with one mask code removed
                    beta = squeeze(RD.([tstim]).MaskTraces{brow}(CellList.cellIDX(cellRows(c)),maskset(k),:)); 
                    leaveone = leaveone +beta;
                end

                leaveone = leaveone +const;
                R = corrcoef(leaveone,T');
                ccByMasktemp(c,j) = R(1,2);
                
          
        end
    end
    
    mean_cc = mean(ccByMasktemp,1);
    clear ccByMasktemp
    %-----

 % remove nan data:
 if ~isnan(mean_cc)

  GroupList = [GroupList; CellList.Group(cellRows(1))];
    Mean_CC_All = [Mean_CC_All;mean_cc];

 end
end

    



%% Look at the CC across regressors/groups
Groups = unique( GroupList);
 leg = {'Motor','Sound', 'Residuals', 'All'};
for G = 1:length(Groups)
    figure;
    gIX = strcmp(GroupList,Groups(G));
    MM = Mean_CC_All(gIX,:);
    for i = 1:size(MM,2)
         cdfplot(MM(:,i)); hold on 
    end
    title(Groups(G))
    legend(leg)
end

%% each regressor new figure; look at groups 

Groups = unique( GroupList);
 leg = {'Motor','Sound', 'Residuals', 'All'};
for i = 1:size(MM,2)
    figure;
for G = 1:length(Groups)

    gIX = strcmp(GroupList,Groups(G));
    MM = Mean_CC_All(gIX,i);
cdfplot(MM); hold on
   
end
 title(leg{i})
    legend(Groups)
end
%% add the betas... different way (order?) of looking at the data
% instead of "leave one out" we are going to keep putting one back in.
% Start with sound to see how well we can model the sound. Then we will add
% in additional Betas until we the "all" cdf

% calculate CC for just sound. Then calculate CC for sound + motor. Then
% Sound + motor + residuals...

% this analysis is a little hard coded in that it runs in this specific
% order

matches = unique(CellList.cellMatch);
Mean_CC_All = [];
GroupList = [];
for m = 1:length(matches)
    cellOI = matches(m);
    cellRows = find(CellList.cellMatch == cellOI);
    for c = 1:length(cellRows)
        tblock = CellList.BlockName{cellRows(c)};
        tstim = CellList.StimType(cellRows(c));
        tdata = RD.([tstim]);
        brow = find(strcmp(tblock,tdata.BlockName));
%         uniqueMaskCode = unique(tdata.params(brow).eventCode);
        uniqueMaskCode = [2,1,5]; % put in the order that I want ... sound+motor+residuals
        T = tdata.TestTrace{brow}(CellList.cellIDX(cellRows(c)),:); % Ca data to test
        Ws =  RD.([tstim]).Filter{brow}; % filter from validation
        mask =  RD.([tstim]).mask{brow};
        maskLabel = RD.([tstim]).maskLabel{brow};
        %-----
        ccByMask = [];
        for j = 1:length(uniqueMaskCode)%(mask,2)
            maskset = 1:size(mask,2); % all the masks in dataset; remove last one because that is the whole model.
            if j == 1
                codeCheck = uniqueMaskCode(j);%RidgeData.params(bl).eventCode(j); % which masks should be grouped together. If codes match
                unmatchCodes = find(codeCheck ~= RD.([tstim]).params(brow).eventCode);
                maskset(unmatchCodes) = []; % leave otu regressors taht arent sounds
            elseif j ==2
                 codeCheck = uniqueMaskCode(j+1);
                 matchCodes = find(codeCheck ~= RD.([tstim]).params(brow).eventCode);
                 maskset(matchCodes) = []; % remove residuals (or whichever matchCode we put third) from the analysis. This essentially just gives you sound +motor
            end
            
            % loop through cells in FOV :
            
            %                 T = testTrace(CellList.cellIDX(cellRows(c)),:); %
            const = Ws(1,CellList.cellIDX(cellRows(c)));
            leaveone = 0;
            for k = 1:length(maskset) % masks with one mask code removed
                beta = squeeze(RD.([tstim]).MaskTraces{brow}(CellList.cellIDX(cellRows(c)),maskset(k),:));
                leaveone = leaveone +beta;
            end
            
            leaveone = leaveone +const;
            R = corrcoef(leaveone,T');
            ccByMasktemp(c,j) = R(1,2);
            
            
        end
    end
    
    mean_cc = mean(ccByMasktemp,1);
    clear ccByMasktemp
    %-----
    
    % remove nan data:
    if ~isnan(mean_cc)
        
        GroupList = [GroupList; CellList.Group(cellRows(1))];
        Mean_CC_All = [Mean_CC_All;mean_cc];
        
    end
end

%% cdf with groups for method above
% plot Group(sound) and Group(whole model) to observe shifts

Groups = unique( GroupList);
 leg = {'NDNF-Beta-Sound','NDNF - Beta-All','PYR-Beta-Sound','PYR - Beta-All','VIP-Beta-Sound','VIP - Beta-All'};
   figure;
for G = 1:length(Groups)
    gIX = strcmp(GroupList,Groups(G));
    MM = Mean_CC_All(gIX,[1 3]);
    for i = 1:size(MM,2)
         cdfplot(MM(:,i)); hold on 
    end
    title(Groups(G))
    legend(leg)
end

%%

FolderName = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\GLM\For Paper Aug2023\All residuals as regressor\Figures all Stim';
% FolderName = save_path;
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = num2str(get(FigHandle, 'Number'));
    set(0, 'CurrentFigure', FigHandle);
    %   saveas(FigHandle, strcat(FigName, '.png'));
    saveas(FigHandle, fullfile(FolderName,strcat(FigName, '.fig'))); % specify the full path
    saveas(FigHandle, fullfile(FolderName,strcat(FigName, '.png'))); % specify the full path
end

close all
