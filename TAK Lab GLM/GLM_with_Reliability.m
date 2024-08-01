% GLM as a function of Reliability:
%------------------
% load Ridge data from GLM
%--------------------
% load data from ExtractedData
%---------------------
Summary = ExtractedData.Summary.All;
 
 
% loop through the blocks - put the CC from GLM into
% ExtractedData.([Locomotor filter]).SummaryData. 
TempCC = nan(height(Summary),1);
for i = 1:height(Summary)
    bname = Summary.Block(i); % block of interest
    bglm = find(strcmp(bname,RidgeData.BlockName));
    %match cell numbers (GLM may not have run on all cells in data set -
    %this will catch that).
    cellnum = Summary.CellNumber(i); % cell num in ExtractedData row
    if ismember(cellnum, RidgeData.params(bglm).cellnum) %check if in Ridge Data
        glmcellrow = find( RidgeData.params(bglm).cellnum == cellnum);
        TempCC(i) = RidgeData.CorrCoeff{bglm,1}(1,glmcellrow);
    end
end
Summary.GLM_CC = TempCC;
ExtractedData.Summary.All = Summary;


figure; scatter(Summary.R1, Summary.GLM_CC);
xlabel('Reliability(R1)')
ylabel('GLM performance')

% sort by cell type:
            Ops = struct;
            Ops.Loco                = 'All'; %All, Running, or NotRunning
            Ops.ResponseType        = 'activated'; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
            Ops.RF_Type             = ''; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
            Ops.IsResponsive        = 1; %1 for responsive (and reliable), 0 for not, 2 for both
            Ops.sortbyGCAMP         = 0; %0 for groups, 1 for gcamp, 2 to combine

   

            SubsetData = simple_subset_ExtractedData(ExtractedData, Ops); %This function filters the data based on the above ops


            Ops = SubsetData.Ops;
%             StimInfo = SubsetData.StimInfo;
            Summary = SubsetData.Summary;
            CellDataByStim = SubsetData.CellDataByStim;
%             StimAnalysis = SubsetData.StimAnalysis;
            Groups = Summary.Group;
            GroupList = unique(Groups);
            colors = {'k','b','r'};
  
            figure;
            for G = 1:length(GroupList)
                groupcells = find(strcmp(Summary.Group, GroupList(G)));
                if isempty(groupcells)
                    continue
                end
                
                if strcmp(GroupList(G), 'PYR')
                    color = 'k';
                elseif strcmp(GroupList(G), 'NDNF')
                    color = 'r';
                elseif strcmp(GroupList(G),'VIP')
                    color = 'b';
                else error('group not included')
                end
                
                scatter(Summary.R2(groupcells),Summary.GLM_CC(groupcells),color); hold on
               
                [Fit,s] = polyfit(Summary.R2(groupcells),Summary.GLM_CC(groupcells),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line
                [f,delta] = polyval(Fit,Summary.R2(groupcells),s); 
%                 yfit = Fit(1)*Summary.R2(groupcells)+Fit(2);
                plot(Summary.R2(groupcells),f,color);
%                 plot(Summary.R2(groupcells),f+2*delta,color,Summary.R2(groupcells),f-2*delta,color); % CI
                
            end
            xlabel('Reliability(R2)')
            ylabel('GLM performance (CorrCoeff)')
            titlename = strcat(ExtractedData.StimInfo.StimType,{' '},Ops.ResponseType,{' '}, ' GLM vs Reliability');
            title(titlename)
            
            
%             
%              x = 1:10; 
%     y1 = x + randn(1,10); 
%     scatter(x,y1,25,'b','*') 
%     P = polyfit(x,y1,1);
%     yfit = P(1)*x+P(2);
%     hold on;
%     plot(x,yfit,'r-.');

%% Leave one out GLM data sorted by cell type/Response type

for G = 1:length(GroupList)
    GroupList(G)
    count = 0;
    for b = 1:height(RidgeData)
        cb = find(strcmp(RidgeData.BlockName(b),Summary.Block));
        if ~isempty(cb)
            if strcmp(Summary.Group(cb(1)),GroupList(G))
                
                count = count+1;
                tempCells = Summary.CellNumber(cb);
                RidgeIdx = find(ismember(RidgeData.params(b).cellnum,tempCells));
                TempRidgeData.BlockName(count) = RidgeData.BlockName(b);
                TempRidgeData.TestTrace{count} = RidgeData.TestTrace{b}(RidgeIdx,:);
                TempRidgeData.params(count) = RidgeData.params(b);
                TempRidgeData.Filter{count} = RidgeData.Filter{b}(RidgeIdx,:);
                TempRidgeData.maskLabel{count} = RidgeData.maskLabel{b};
                TempRidgeData.mask{count} = RidgeData.mask{b};
                TempRidgeData.MaskTraces{count} = RidgeData.MaskTraces{b}(RidgeIdx,:,:);
                TempRidgeData.CorrCoeff{count} = RidgeData.CorrCoeff{b}(1,RidgeIdx);
            end %group list == block group
            clear cb
        end % cb empty
  
    end
    GroupRidge{G} = TempRidgeData;
    clear TempRidgeData
end

for G = 1:length(GroupList)
    GLM_leave_one_out(GroupRidge{G},0);
    temptitle = strcat(ExtractedData.StimInfo.StimType,{' '}, GroupList(G));
    title(temptitle)
end



