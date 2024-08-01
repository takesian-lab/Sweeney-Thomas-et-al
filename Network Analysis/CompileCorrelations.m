%% CompileCorrelations

% Tool to compile output tables from Network Analysis into single
% spreadsheet

% Adds a few columns to identify + and - correlations for easy
% incorporation into statistical program (JMP, excel)

% If user specifies CompileOps.avg_across_stim = 1, will append data for
% each matched pair averaged across all stim (stim type = 'AllStim')


%% Setup Compile Ops Variables

close all
clear

%Show user dialog box to prompt them to select Ops.m file
% see default template for ops files in: ...GitHub\Calcium-Imaging-Analysis\Network Analysis\Network Analysis Ops Files    

message = sprintf('Select a Network Ops file, check settings, and hit space bar');
    uiwait(msgbox(message));
    
[Network_OpsFile, Network_OpsPath] = uigetfile;

cd(Network_OpsPath)     %Go to file location
edit(Network_OpsFile)   %Open the file for the user to edit and pause while they do son
pause;                  %Wait for user to hit enter in command line
run(Network_OpsFile)    %Run the Ops.m file once they do so


Compiled_Table = [];

%% Load tables and vertically combine
for i=1:length(CompileOps.stimTypes)
    D = dir(fullfile(CompileOps.data_path,'*.mat'));
    N = {D.name};
    s = CompileOps.stimTypes{1,i}; % pulls out files in Extract Data folder containing stimulus names
    name = strcat('Correlations_Table_', s,'.mat');
    X = strcmp(N,name); %contains(N,name); %MET replaced contains with strcmp 8/3/23 so that string must be an exact match
    
    if ~any(X) %Skip stim type if not found in directory
        warning(['No Table file for stim type ' s ' found. Skipping...'])
        continue;
    elseif sum(X) > 1 %Only load first file if more than one found
        warning(['More than one Table file for stim type ' s ' found. Using first one...'])
        first_ind = find(X,1,'first');
        X(:) = 0;
        X(first_ind) = 1;
    end
    
    cd(CompileOps.data_path);
    load(N{X});
    Compiled_Table = vertcat(Compiled_Table,Correlations_Matrix_All);
end

%% Add columns to table indicating + and - correlations

for r = 1:6
    if r == 1
        T = Compiled_Table.TraceMaxCorrelations; SigT = Compiled_Table.TraceMaxZTest;
    elseif r == 2
        T = Compiled_Table.TraceZeroCorrelations; SigT = Compiled_Table.TraceZeroZTest;
    elseif r == 3
        T = Compiled_Table.NoiseCorrelations; SigT = Compiled_Table.NoiseZTest;
    elseif r == 4
        T = Compiled_Table.SignalCorrelations; SigT = Compiled_Table.SignalZTest;
    elseif r == 5
        T = Compiled_Table.NoiseXCorrelations; SigT = Compiled_Table.NoiseXZTest;
    elseif r == 6
        T = Compiled_Table.SignalXCorrelations; SigT = Compiled_Table.SignalXZTest;
    end

    % columns for noise corr
    pos_neg = nan(length(T),1);
    neg_index = find(T<0);
    pos_neg(neg_index) = -1;
    pos_index = find(T>0);
    pos_neg(pos_index) = 1;

    sig_index = find(SigT==1);

    sig = nan(length(T),1);
    sig(intersect(neg_index,sig_index))=-1;
    sig(intersect(pos_index,sig_index))=1;

    if r == 1
        Compiled_Table.PosOrNegTraceMaxCorr = pos_neg;
        Compiled_Table.SigPosOrNegTraceMaxCorr = sig;
    elseif r == 2
        Compiled_Table.PosOrNegTraceZeroCorr = pos_neg;
        Compiled_Table.SigPosOrNegTraceZeroCorr = sig;
    elseif r == 3
        Compiled_Table.PosOrNegNoiseCorr = pos_neg;
        Compiled_Table.SigPosOrNegNoiseCorr = sig;
    elseif r == 4
        Compiled_Table.PosOrNegSignalCorr = pos_neg;
        Compiled_Table.SigPosOrNegSignalCorr = sig;
    elseif r == 5
        Compiled_Table.PosOrNegNoiseXCorr = pos_neg;
        Compiled_Table.SigPosOrNegNoiseXCorr = sig;
    elseif r == 6
        Compiled_Table.PosOrNegSignalXCorr = pos_neg;
        Compiled_Table.SigPosOrNegSignalXCorr = sig;
    end

end
%% Add pair match label

%Create unique matched pair labels
Cell1_Match_Sorted = Compiled_Table.Cell1_Match;
Cell2_Match_Sorted = Compiled_Table.Cell2_Match;
C2first = Cell1_Match_Sorted > Cell2_Match_Sorted; %wherever cel11 ID > cell2 ID, swap labels
Cell1_Match_Sorted(C2first) = Compiled_Table.Cell2_Match(C2first);
Cell2_Match_Sorted(C2first) = Compiled_Table.Cell1_Match(C2first);

%Combine match labels into a unique string and store
CombineMatch = strings(size(Cell1_Match_Sorted)); 
for i = 1:length(CombineMatch)
    CombineMatch(i) = strcat(num2str(Cell1_Match_Sorted(i)), {'_'}, num2str(Cell2_Match_Sorted(i)));
end
Compiled_Table.CombineMatch = CombineMatch;

%% Save as csv and matlab file

save('CompiledTable.mat', 'Compiled_Table') ;
writetable(Compiled_Table,'CompiledTable.csv','Delimiter',',');

%% Find unique matched pairs and average across stim

if CompileOps.avg_across_stim

    tic

    %Remove Match cell NaNs
    nanind = (isnan(Compiled_Table.Cell1_Match) + isnan(Compiled_Table.Cell2_Match)) > 0;
    Compiled_Table(nanind,:) = [];

    %Preallocate table
    Loco = unique(Compiled_Table.LocoType); %Sort by Locomotion
    pairs_by_loco = cell(1,length(Loco)); %find all matched numbers per loco condition
    nRows = 0;
    for L = 1:length(Loco)
        pairs_by_loco{L} = unique(Compiled_Table.CombineMatch(strcmp(Compiled_Table.LocoType,Loco{L})));
        nRows = nRows + length(pairs_by_loco{L});
    end
    columnnames = Compiled_Table.Properties.VariableNames;
    columntypes = varfun(@class,Compiled_Table,'OutputFormat','cell');
    sz = [nRows length(columnnames)];
    T = table('Size',sz,'VariableTypes',columntypes,'VariableNames',columnnames);

    count = 1; %combined pair count

    for L = 1:length(Loco)

        currentPairs = pairs_by_loco{L};

        for p = 1:length(currentPairs)

            ind_pair = intersect(find(strcmp(Compiled_Table.CombineMatch, currentPairs(p))),find(strcmp(Compiled_Table.LocoType,Loco{L})));

            if isempty(ind_pair)
                continue;
            end

            % if pair exists
            T.MouseName(count)              = Compiled_Table.MouseName(ind_pair(1));
            T.BlockName(count)              = Compiled_Table.BlockName(ind_pair(1));
            T.FieldofView(count)            = Compiled_Table.FieldofView(ind_pair(1));
            T.StimName(count)               = "AverageAllStim";
            T.CellType(count)               = Compiled_Table.CellType(ind_pair(1));
            T.LocoType(count)               = Compiled_Table.LocoType(ind_pair(1));
            T.Cell1_s2P(count)              = Compiled_Table.Cell1_s2P(ind_pair(1)); %C1/C2 correspond to first found pair
            T.Cell2_s2P(count)              = Compiled_Table.Cell2_s2P(ind_pair(1)); %C1/C2 correspond to first found pair
            T.Cell1_Match(count)            = Compiled_Table.Cell1_Match(ind_pair(1)); %C1/C2 correspond to first found pair
            T.Cell2_Match(count)            = Compiled_Table.Cell2_Match(ind_pair(1)); %C1/C2 correspond to first found pair
            T.Cell1_RedCell(count)          = Compiled_Table.Cell1_RedCell(ind_pair(1)); %C1/C2 correspond to first found pair
            T.Cell2_RedCell(count)          = Compiled_Table.Cell2_RedCell(ind_pair(1)); %C1/C2 correspond to first found pair
            T.Distance(count)               = Compiled_Table.Distance(ind_pair(1));
            T.TraceMaxCorrelations(count)   = mean(Compiled_Table.TraceMaxCorrelations(ind_pair),'omitnan');
            T.TraceZeroCorrelations(count)  = mean(Compiled_Table.TraceZeroCorrelations(ind_pair),'omitnan');
            T.TraceMaxZTest(count)          = max(Compiled_Table.TraceMaxZTest(ind_pair));
            T.TraceZeroZTest(count)         = max(Compiled_Table.TraceZeroZTest(ind_pair));
            T.NoiseCorrelations(count)      = mean(Compiled_Table.NoiseCorrelations(ind_pair),'omitnan');
            T.SignalCorrelations(count)     = mean(Compiled_Table.SignalCorrelations(ind_pair),'omitnan');
            T.NoiseZTest(count)             = max(Compiled_Table.NoiseZTest(ind_pair));
            T.SignalZTest(count)            = max(Compiled_Table.SignalZTest(ind_pair));
            T.NoiseXCorrelations(count)     = mean(Compiled_Table.NoiseXCorrelations(ind_pair),'omitnan');
            T.SignalXCorrelations(count)    = mean(Compiled_Table.SignalXCorrelations(ind_pair),'omitnan');
            T.NoiseXZTest(count)            = max(Compiled_Table.NoiseXZTest(ind_pair));
            T.SignalXZTest(count)           = max(Compiled_Table.SignalXZTest(ind_pair));
            
            %Significant pos/negative
            T.PosOrNegTraceMaxCorr(count)       = sign(T.TraceMaxCorrelations(count));
            T.SigPosOrNegTraceMaxCorr(count)    = nan;
            T.PosOrNegTraceZeroCorr(count)      = sign(T.TraceZeroCorrelations(count));
            T.SigPosOrNegTraceZeroCorr(count)   = nan;
            T.PosOrNegNoiseCorr(count)          = sign(T.NoiseCorrelations(count));
            T.SigPosOrNegNoiseCorr(count)       = nan;
            T.PosOrNegSignalCorr(count)         = sign(T.SignalCorrelations(count));
            T.SigPosOrNegSignalCorr(count)      = nan;
            T.PosOrNegNoiseXCorr(count)         = sign(T.NoiseXCorrelations(count));
            T.SigPosOrNegNoiseXCorr(count)      = nan;
            T.PosOrNegSignalXCorr(count)        = sign(T.SignalXCorrelations(count));
            T.SigPosOrNegSignalXCorr(count)     = nan;

            if T.TraceMaxZTest(count) == 1;     T.SigPosOrNegTraceMaxCorr(count)    = T.PosOrNegTraceMaxCorr(count);    end
            if T.TraceZeroZTest(count) == 1;    T.SigPosOrNegTraceZeroCorr(count)   = T.PosOrNegTraceZeroCorr(count);   end
            if T.NoiseZTest(count) == 1;        T.SigPosOrNegNoiseCorr(count)       = T.PosOrNegNoiseCorr(count);       end
            if T.SignalZTest(count) == 1;       T.SigPosOrNegSignalCorr(count)      = T.PosOrNegSignalCorr(count);      end
            if T.NoiseXZTest(count) == 1;       T.SigPosOrNegNoiseXCorr(count)      = T.PosOrNegNoiseXCorr(count);      end
            if T.SignalXZTest(count) == 1;      T.SigPosOrNegSignalXCorr(count)     = T.PosOrNegSignalXCorr(count);     end

            %Match label
            T.CombineMatch(count)               = Compiled_Table.CombineMatch(ind_pair(1));

            count = count+1;
        end
    end

    Compiled_Table_AvgStim = T;

    %% Save as csv and matlab file
    save('CompiledTable_AverageStim.mat', 'Compiled_Table_AvgStim') ;
    writetable(Compiled_Table_AvgStim,'CompiledTable_AverageStim.csv','Delimiter',',');

    toc
end

%%
disp('done analyzing network...')