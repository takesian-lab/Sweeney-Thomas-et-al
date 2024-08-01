function  [params] = GLM_simple_setRegs(ExtractedData,D,bl,ztraces,Ops,bsfr, Compiled_Table)

% params.eventcode is used to group together similar type stim for masks. %
% codes are:
% motor = 1
% sound = 2
% outcome = 3
% cell activity = 4
% Residuals  = 5; 
% Average activity from positive/negative correlated neurons = 7; To
% seperate out pos from neg activity, set negative to 8

%------------------------------------------------------------------
% determine if these data are from FibPhot block
%------------------------------------------------------------------
if strcmp(ExtractedData.BlockInfo.AnalysisPath(bl),'FibPhot')
    FP = true;
    FPch = Ops.FibPhotChannel;
else
    FP = false;
end

%------------------------------------------------------------------
% set trial times for splitting data into testing and training sets -
%------------------------------------------------------------------
ts = ExtractedData.TrialData(bl).Timestamp; % timestamps for the experiment
sizebout = Ops.dataChunk*ExtractedData.BlockInfo.Framerate(bl); % size of chunk of data
divTime = floor(length(ts)/sizebout); % break up timestamp into groups for analysis.

if Ops.ByTrials
    divTime = numel(find(ExtractedData.TrialData(bl).IsBlank ==0));
    sizebout= size(ExtractedData.TrialData(bl).Trials,3);
    start = 0;
    for i = 1:divTime
        ttime(i,1) = start+1;
        ttime(i,2) = start+sizebout;
        start =  ttime(i,2);
    end
else
    start = 1;
    sizebout = Ops.dataChunk*ExtractedData.BlockInfo.Framerate(bl); % size of chunk of data
    divTime = floor(length(ts)/sizebout); % break up timestamp into groups for analysis.
    for i = 1:divTime
        ttime(i,1) = start;
        ttime(i,2) = start+(sizebout-5); %todo: magic number - 5 frames works such that events dont overlap, but that we dont loose too much data
        start = start+sizebout;
    end
end
params.trigTime = ttime;
%--------------------------------------------
% determine which rows (i.e. cell traces) to use from extracted data.
% for active cells, they are chosen as those that are
% activated/prolonged/suppressed in extracted data.
%--------------------------------------------
if FP
    params.cellnum = ExtractedData.CellList.CellNumber(D.rows(FPch));
    params.cellIDX = ExtractedData.CellList.CellOrder(D.rows(FPch));
else
    if Ops.activeCells
        if Ops.MultipleActiveStim
            active_cells = [];
            for d = 1:length(D.rows)
                RFstim = sum(ExtractedData.CellDataByStim.All(D.rows(d)).RF);
                if RFstim>=Ops.stimnumthesh;
                    active_cells = [active_cells;d];
                end
            end
            
            
        else
            active_cells = find(ExtractedData.CellData.All.ResponseType(D.rows) ~= 'none');
        end
        
        params.cellnum = ExtractedData.CellList.CellNumber(D.rows(active_cells));
        params.cellIDX = ExtractedData.CellList.CellOrder(D.rows(active_cells));
        params.cellMatch = ExtractedData.CellList.MatchedRow(D.rows(active_cells));
    else
        params.cellnum = ExtractedData.CellList.CellNumber(D.rows);
        params.cellOrder = ExtractedData.CellList.CellOrder(D.rows);
        params.cellIDX = ExtractedData.CellList.CellOrder(D.rows);
        params.cellMatch = ExtractedData.CellList.MatchedRow(D.rows);
    end
end

%--------------------------------------------
FF = 1; % FF will keep track of number of events that go into the design matrix
% 
if Ops.useLocomotor
    if Ops.ByTrials
        CatLoc = [];
        for i = 1:size(ExtractedData.TrialData(bl).Trials,2)
            if ExtractedData.TrialData(bl).IsBlank(i) == 0
                locdat = ExtractedData.TrialData(bl).Loco(i,:);
                if Ops.NormLoco
                    locbase = locdat(bsfr);
                    bm = mean(locbase,1,'omitnan');
                    bst = std(locbase);
                    bst(bst==0) = 1; % correct if std ==0
                    newlocdata = (locdat - bm)./bst;
                elseif Ops.Acceleration
                    zeroloc = zeros(size(locdat));
                    zeroloc(2:end) = diff(locdat);
                    newlocdata = zeroloc;
                elseif Ops.LocoBin
                   locdat(locdat>0) = 1;
                   newlocdata = locdat;
                else
                    newlocdata = locdat;
                end
                CatLoc = [CatLoc,newlocdata];
                
            end
        end
        params.eventTime{1,FF} = CatLoc;
       
    else
      if Ops.NormLoco          
          params.eventTime{1,FF} = ExtractedData.TrialData(bl).Full_Loco./max(ExtractedData.TrialData(bl).Full_Loco);
      else
        params.eventTime{1,FF} = ExtractedData.TrialData(bl).Full_Loco;
      end
    end
    params.eventLabel{1,FF} = 'Motor Activity';
    params.eventCode(1,FF) = 1;
    
    
    if Ops.LocoBin & ~ Ops.PredictMotor
        params.eventTime{1,FF}(params.eventTime{1,FF}>0) = 1;
    end
    FF = FF+1;
end
%--------------------------------------------
% if Ops.usePupil
%     if Ops.ByTrials
%         CatLoc = [];
%         for i = 1:size(ExtractedData.TrialData(bl).Trials,2)
%             if ExtractedData.TrialData(bl).IsBlank(i) == 0
%                 locdat = ExtractedData.TrialData(bl).Loco(i,:);
%                 if Ops.NormLoco
%                     locbase = locdat(bsfr);
%                     bm = mean(locbase,1,'omitnan');
%                     bst = std(locbase);
%                     bst(bst==0) = 1; % correct if std ==0
%                     newlocdata = (locdat - bm)./bst;
%                 elseif Ops.Acceleration
%                     zeroloc = zeros(size(locdat));
%                     zeroloc(2:end) = diff(locdat);
%                     newlocdata = zeroloc;
%                 elseif Ops.LocoBin
%                    locdat(locdat>0) = 1;
%                    newlocdata = locdat;
%                 else
%                     newlocdata = locdat;
%                 end
%                 CatLoc = [CatLoc,newlocdata];
%                 
%             end
%         end
%         params.eventTime{1,FF} = CatLoc;
%        
%     else
%       if Ops.NormLoco          
%           params.eventTime{1,FF} = ExtractedData.TrialData(bl).Full_Loco./max(ExtractedData.TrialData(bl).Full_Loco);
%       else
%         params.eventTime{1,FF} = ExtractedData.TrialData(bl).Full_Loco;
%       end
%     end
%     params.eventLabel{1,FF} = 'Motor Activity';
%     params.eventCode(1,FF) = 1;
%     
%     

%--------------------------------------------
if Ops.useLicks
    if Ops.ByTrials
        CatLick = [];
        for i = 1:size(ExtractedData.TrialData(bl).Trials,2)
            if ExtractedData.TrialData(bl).IsBlank(i) == 0
                lickdat = full(ExtractedData.TrialData(bl).Licks(i,:));
                CatLick = [CatLick,lickdat];
            end
        end
        params.eventTime{1,FF} = CatLick;
        
    else
      params.eventTime{1,FF} = ExtractedData.TrialData(bl).Full_licks;
    end
end
params.eventLabel{1,FF} = 'Licks';
params.eventCode(1,FF) = 9;
FF = FF+1;


%--------------------------------------------
if Ops.useSounds
    if Ops.ByTrials
        times = params.trigTime(:,1) + bsfr-1;
    else
        times = ExtractedData.TrialData(bl).Sound_Time(D.stimIDX);
    end
    params.eventTime{1,FF} = times;
    
    params.eventLabel{1,FF} = 'All Sounds';
    params.eventCode(1,FF) = 2;
    FF = FF+1;
end
%--------------------------------------------
if Ops.useV1
    V1list = D.uniqueV1;
    if Ops.useMag %use value of V1 as regressor - makes one regressor
        if Ops.ByTrials
            times = params.trigTime(:,1) + bsfr-1;
        else
            times = ExtractedData.TrialData(bl).Sound_Time(D.stimIDX);
        end
        params.eventTime{1,FF} = times;
        
        params.eventLabel{1,FF} = 'All V1; Analog';
        params.eventCode(1,FF) = 2;
        FF = FF+1;
    else % use binary indicator for V1, each unique V1 will be a seperate regressor
        for i = 1:length(V1list)
            if Ops.ByTrials
                params.eventTime{1,FF} = params.trigTime(D.V1_IDX{i},1) + bsfr-1;
                params.eventLabel{1,FF} = num2str(V1list(i));
            else
                params.eventTime{1,FF} = ExtractedData.TrialData(bl).Sound_Time(D.V1_IDX{i});
                params.eventLabel{1,FF} = num2str(V1list(i));
            end
            params.eventCode(1,FF) = 2;
            FF = FF+1;
        end
    end
    
end
%--------------------------------------------
if Ops.useV2
    V2list = D.uniqueV2;
    if Ops.useMag %use value of V2 as regressor - makes one regressor
        if Ops.ByTrials
            times = params.trigTime(:,1) + bsfr-1;
        else
            times = ExtractedData.TrialData(bl).Sound_Time(D.stimIDX);
        end
        
        params.eventTime{1,FF} = times;
        params.eventLabel{1,FF} = 'All V2; Analog';
        params.eventCode(1,FF) = 2;
        FF = FF+1;
    else % use binary indicator for V2, each unique V1 will be a seperate regressor
        for i = 1:length(V2list)
            if Ops.ByTrials
                params.eventTime{1,FF} = params.trigTime(D.V2_IDX{i},1) + bsfr -1;
                params.eventLabel{1,FF} = num2str(V2list(i));
                if Ops.useOutcome & ExtractedData.Ops.stim_protocol ==7
                     params.eventCode(1,FF) = 3;
                else
                    params.eventCode(1,FF) = 2;
                end
            else
                params.eventTime{1,FF} = ExtractedData.TrialData(bl).Sound_Time(D.V2_IDX{i});
                params.eventLabel{1,FF} = num2str(V2list(i));
                
                 if Ops.useOutcome & ExtractedData.Ops.stim_protocol ==7
                      params.eventCode(1,FF) = 3;
                else
                    params.eventCode(1,FF) = 2;
                end
            end
            FF = FF+1;
        end
    end
    
end

%------------------------------------------------
if Ops.useOutcome % for freqdisc data, Outcome = V2. currently accounted for in that section. Future updates my have a seperate part here. 
    if ExtractedData.Ops.stim_protocol ~=7
    error('Outcome programmed for FreqDisc analysis only. Talk to Carolyn. March 2024')
    params.eventCode(1,FF) = 3;
FF = FF+1
end
%------------------------------------------------
% if Ops.useCorrCell && Ops.CorrMean
%    for f = 1:2
%        if f ==1
%            params.eventCode(1,FF) = 7; % positive
%        else
%            params.eventCode(1,FF) = 7; % negative - changed from 8 to 7 so that we can look at regressors togeter
%        end
%        FF = FF +1;
%    end
% end
%------------------------------------------------
if Ops.useNeighbor
    error( 'output of GLM_ExDat was updated, and neighboring  neuron anlalysis requires updates')
%     params.eventCode(1,FF) = 4;
    %     distance_in_microns = find_cell_distance(ExtractedData,0,D.rows,bl);
    %     [B,I] = sort(distance_in_microns,1);
    %     numnearcell = round(length(D.rows).*(localCellPer./100));
    %     neighbors = I(2:numnearcell+1,:); % index that matches how block is set up
    % %     ztraces = zscore(block.df_f,0,2);
    %     for i = 1:length(params.cellnum)
    %         LN = neighbors(:,params.cellIDX(i));
    %         LNtraces =  ztraces(LN,:);
    %         if useNeighborMean % average activity on neighboring cells
    %             if i ==1
    %                 FF = FF+1;
    %             end
    %             params.Cont(i).traces(1,:) = mean(LNtraces,1);
    %             params.Cont(i).eventLabel{1,:} = 'Neighboring Cell Activity';
    %
    %         else
    %             for j = 1:size(LNtraces,1) % individual activity of neighboring cells
    %                 if i ==1
    %                     FF = FF+1;
    %                 end
    %                 s2pcellnum = num2str(params.cellnum(LN(j)));
    %                 params.Cont(i).traces(j,:) = LNtraces(j,:);
    %                 params.Cont(i).eventLabel{j,:} = strcat('Cell number_',s2pcellnum );
    %
    %             end
    %         end
    %     end
end
%------------------------------------------------

if Ops.useAllCells 
    for ii = 1:size(ztraces,1)-1
            params.eventLabel{1,FF} = 'All Neurons in FOV';
            params.eventCode(1,FF) = 4;
            FF = FF +1;
    end
    
end

%------------------------------------------------
% set params.eventcode and params.event label
if Ops.useCorrCell
    if Ops.CorrMean;
        params.eventLabel{1,FF}  = 'pos corr neuron';
        params.eventCode(1,FF) = 7;
        FF = FF+1;
        params.eventLabel{1,FF}  = 'neg corr neuron';
        params.eventCode(1,FF) = 7;
        FF = FF +1;
    elseif Ops.CorrResid
         for ii = 1:size(ztraces,1)-1
            params.eventLabel{1,FF} = 'Residual traces';
            params.eventCode(1,FF) = 5;
            FF = FF+1;
         end
    end
end

%------------------------------------------------
if Ops.Residuals 
    for ii = 1:size(ztraces,1)-1
            params.eventLabel{1,FF} = 'Residual traces';
            params.eventCode(1,FF) = 5;
            FF = FF+1;
    end
end

% if Ops.CorrResid
%     if Ops.CorrMean
%         for ii = 1:2
%             params.eventLabel{1,FF} = 'Residual traces';
%             params.eventCode(1,FF) = 5;
%             FF = FF+1;
%         end
%     else
%         for ii = 1:size(ztraces,1)-1
%             params.eventLabel{1,FF} = 'Residual traces';
%             params.eventCode(1,FF) = 5;
%             FF = FF+1;
%         end
%     end
% end



end