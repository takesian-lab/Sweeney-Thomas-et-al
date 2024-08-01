% Network Analysis Pipeline
% Written by Anne Takesian, modified code from Maryse Thomas, Carolyn Sweeney, Lucas Vattino & Rahul Brito
% version 8: Jan 2023 to make compatible with simple extracted data and
% model after Maryse's XCorr for faster trace analyses

% This pipeline analyzes the spontaneous, signal and noise correlations
% from extracted data

% Inputs - loads the extracted data file, blocks and matching file basedF
% on user-defined path information for a given block

% Output -
% Several graphs are saved to the save_path, organized in foldersF
% (block, stimType)
% A correlations table and spreadsheet

%% Setup NetworkOps Variables

close all
clear

%Show user dialog box to prompt them to select Ops.m file
% see default template for ops files in: ...GitHub\Calcium-Imaging-Analysis\Network Analysis\Network Analysis Ops Files   

message = sprintf('Select a Network Ops file, check settings, and hit space bar');
    uiwait(msgbox(message));
    
[Network_OpsFile, Network_OpsPath] = uigetfile;

cd(Network_OpsPath)     %Go to file location
edit(Network_OpsFile)   %Open the file for the user to edit and pause while they do so
pause;                  %Wait for user to hit enter in command line
run(Network_OpsFile)    %Run the Ops.m file once they do so

disp('Loading files...')

%% Load extracted data
for i=1:length(NetworkOps.stimTypes)
    D = dir(fullfile(NetworkOps.data_path,'*.mat'));
    N = {D.name};
    s = NetworkOps.stimTypes{1,i}; % pulls out files in Extract Data folder containing stimulus names
    name = strcat('ExtractedData_', s, '_');
    X = contains(N,name);
    
    if ~any(X) %Skip stim type if not found in directory
        warning(['No ExtractedData file for stim type ' s ' found. Skipping...'])
        continue;
    elseif sum(X) > 1 %Only load first file if more than one found
        warning(['More than one ExtractedData file for stim type ' s ' found. Using first one...'])
        first_ind = find(X,1,'first');
        X(:) = 0;
        X(first_ind) = 1;
    end
    
    cd(NetworkOps.data_path);
    load(N{X});

    % check if matched data are in ExtractedData
    if ~isempty(NetworkOps.matchFile_path)
        ExtractedData.CellList = add_matches_to_CellList(ExtractedData.CellList, NetworkOps.matchFile_path); % Maryse updated to use add_matches_to_CellList 4/3/2023
    else
        %If MatchedRow column doesn't exist, fill with NaNs
        if ~any(strcmp(ExtractedData.CellList.Properties.VariableNames, 'MatchedRow'))
            ExtractedData.CellList.MatchedRow = nan(size(ExtractedData.CellList,1),1);
            disp('No match data found. Filling with NaNs')
        end
    end
    Data.(s)=ExtractedData;
end

if ~exist('Data','var')
    error('No ExtractedData loaded');
end

%% Initialize correlations table
Correlations_Matrix_All = table; % table with all spont, noise and signal correlations data
MotorCorrelations_Matrix_All = table; % table with all motor correlations data
Correlations = struct; % structure with block correlations and info

%% Analyze Correlations - loop through all stim types
for i=1:length(NetworkOps.stimTypes)
    s = NetworkOps.stimTypes{1,i};

    % Create data structure for files corresponding to stim_protocol
    disp('Loading blocks to include...')
    blockData = Data.(s).BlockInfo;

    % Loop through each mouse and block
    mouse_names = unique(blockData.MouseID); % get mouse names

    for m = 1:size(mouse_names,1) % loop through mice
        mouse_name = mouse_names{m}; % mouse name
        current_rows = find(blockData.MouseID == mouse_name); % find rows with the mouse name
        block_num = numel(unique(blockData.Block(current_rows))); % number of blocks to analyze/mouse
        block_names =unique(blockData.Block(current_rows)); % get block names for given mouse

        for b = 1:block_num % loop through blocks
            fig_count = 1; % start figure count
            block_name = block_names{b}; % block name
            disp(strcat({'Analyzing '}, mouse_name, {' '}, block_name)); % display block name in output
            block_row = find(blockData.Block == block_name); % block row number
            cell_type = blockData.Group(block_row); % define cell type
            FOV = blockData.FOV(block_row); % define FOV for block
            all_cell_rows = find(Data.(s).CellList.Block == block_name); % rows of all cells within block 
            loco_type = fieldnames(Data.(s).CellData); % All, Running, NotRunning

            % find distances for all cells
            distances_all = find_cell_distance(Data.(s),'Rows',all_cell_rows,'BlockRow',block_row);

            % find cells across locomotion conditions
            all_ind = find(strcmp(loco_type, 'All'));
            running_ind = find(strcmp(loco_type, 'Running')); 
            notrunning_ind = find(strcmp(loco_type, 'NotRunning'));
            
            for n=1:length(loco_type)
                if NetworkOps.use_responsive_cells && ~contains(s,'pontaneous')  % if finding responsive cells & not spont block
                    tempRF = Data.(s).CellData.([loco_type{n}]).RF;
                    R_cellIDX = find(tempRF==1);
                    cell_rows{n} = intersect(R_cellIDX,all_cell_rows);
                    responsive_cells_block{n} = (cell_rows{n}-all_cell_rows(1))+1;
                else
                    cell_rows{n} =  all_cell_rows;
                    responsive_cells_block{n} = (cell_rows{n}-all_cell_rows(1))+1;
                end
                distances{n} = find_cell_distance(Data.(s),'Rows',cell_rows{n},'BlockRow',block_row);
            end

            % find times within trace with running for each analyzed (responsive) neuron
            activity = Data.(s).TrialData(block_row).Full_Loco; % loco activity in the block
            activity(activity < 0.8) = 0; %Correct for noise floor
            
            %all traces
            if ~isempty(all_ind)
                traces{all_ind} = Data.(s).TrialData(block_row).F7(responsive_cells_block{1},:); % all traces
            end
            
            %running
            if ~isempty(running_ind)
                [temp,traces{running_ind}] = find_loco_traces(Data.(s).TrialData(block_row).F7(responsive_cells_block{2},:), NetworkOps.locoThresh, activity); % running
            end
            
            %not running
            if ~isempty(notrunning_ind)
                [traces{notrunning_ind},temp] = find_loco_traces(Data.(s).TrialData(block_row).F7(responsive_cells_block{3},:), NetworkOps.locoThresh, activity); % no running
            end
            
            %% TRACE CORRELATIONS: find pair-wise cross correlations of entire traces
            data = {}; %Initialize data: contains info about which corr results to analyze

            if NetworkOps.Xcorr.ComputeXcorr
                data = [data, {'XCorrTraces'; '.cc_max'; '.cc_max_z'}, {'XCorrTraces'; '.cc_zero'; '.cc_zero_z'}]; %Concatenate for each subsequent analysis type
                NetworkOps.XCorr.fs = ExtractedData.StimInfo.fs; % frame rate
                for n=1:length(loco_type) % 1 or 3 depending on loco filter used in extract data
                    if ~isempty(traces{n})
                        [Correlations.XCorrMat{n}, Correlations.XCorrLag{n}, Correlations.XCorrTraces{n}, ~] = tak_compute_xcorr(traces{n}, traces{n}, NetworkOps.XCorr.maxlag_in_s, NetworkOps.XCorr.fs, 'Zscore', NetworkOps.XCorr.zscore, 'nShuffles', NetworkOps.XCorr.shuffles, 'pvalue', NetworkOps.XCorr.p_value);
                        %Maryse wrote above script to standardize across different analyses
                        %[Correlations.XCorrTraces{n}, Correlations.XCorrMat{n}, Correlations.XCorrLag{n}] = simple_xcorr_network(traces{n}, NetworkOps.XCorr);
                    end
                end
            else
                disp('Skipping Trace Correlation Analysis...');
            end

            %% LOCO CORRELATIONS: if motor activity, correlate traces with motor activity
            if NetworkOps.MotorCorr.ComputeMotorCorr
                
                PlotMotorFigure = 1; %Keep track of when the motor figure gets plotted so we don't plot it 4 times
                NetworkOps.MotorCorr.fs = ExtractedData.StimInfo.fs; % frame rate
                
                [hasLoco, hasXY, hasZcorr] = deal(false);
                columns = fields(Data.(s).TrialData);
                if any(ismember(columns, 'Loco')); hasLoco = ~isempty(Data.(s).TrialData(block_row).Loco); end
                if any(ismember(columns, 'Full_xoff')); hasXY = ~isempty(Data.(s).TrialData(block_row).Full_xoff); end
                if any(ismember(columns, 'Full_zcorr')); hasZcorr = ~isempty(Data.(s).TrialData(block_row).Full_zcorr); end

                if ~hasLoco
                    disp('Skipping Motor Correlation Analysis...');
                else
                    data = [data, {'MotorCorr'; '.cc_max'; '.cc_max_z'}, {'MotorCorr'; '.cc_zero'; '.cc_zero_z'}]; %Concatenate for each subsequent analysis type
                    [figures{1,fig_count}, Correlations.MotorCorr{1}, Correlations.MotorCorrMat{1}, Correlations.MotorCorrLag{1}] = ...
                        simple_motor_correlation(traces{all_ind}, Data.(s), NetworkOps.MotorCorr, block_row, 'TraceToAnalyze', 'Loco', 'PlotFigure', PlotMotorFigure);
                    fig_count = fig_count + 1;
                    figureNames(fig_count) = {'Motor'};
                    PlotMotorFigure = 0; %Set to 0 once it has been plotted for the first time
                end       

                % correlate with brain movement (suite2p x/y correction)
                if ~hasXY
                    disp('Skipping X/Y Movement Correlation Analysis...');
                else
                    data = [data, {'XYCorr'; '.cc_max'; '.cc_max_z'}, {'XYCorr'; '.cc_zero'; '.cc_zero_z'}]; 
                    [figures{1,fig_count}, Correlations.XYCorr{1}, Correlations.XYCorrMat{1}, Correlations.XYCorrLag{1}] = ...
                        simple_motor_correlation(traces{all_ind}, Data.(s), NetworkOps.MotorCorr, block_row, 'TraceToAnalyze', 'XY', 'PlotFigure', PlotMotorFigure);
                    fig_count = fig_count + 1;
                    figureNames(fig_count) = {'XY'};
                    PlotMotorFigure = 0; %Set to 0 once it has been plotted for the first time
                end  
                
                % correlate with residual of loco & x/y movement (possible control for running)
                if ~hasLoco || ~hasXY
                    disp('Skipping Loco Residuals Movement Correlation Analysis...');
                else
                    data = [data, {'ResidualsCorr'; '.cc_max'; '.cc_max_z'}, {'ResidualsCorr'; '.cc_zero'; '.cc_zero_z'}]; 
                    [figures{1,fig_count}, Correlations.ResidualsCorr{1}, Correlations.ResidualsCorrMat{1}, Correlations.ResidualsCorrLag{1}] = ...
                        simple_motor_correlation(traces{all_ind}, Data.(s), NetworkOps.MotorCorr, block_row, 'TraceToAnalyze', 'Residuals', 'PlotFigure', PlotMotorFigure);
                    fig_count = fig_count + 1;
                    figureNames(fig_count) = {'LocoResiduals'};
                    PlotMotorFigure = 0; %Set to 0 once it has been plotted for the first time
                end  
                
                % correlate with brain movement (suite2p zcorr)
                if ~hasZcorr
                    disp('Skipping Z Movement Correlation Analysis...');
                else
                    data = [data, {'ZCorr'; '.cc_max'; '.cc_max_z'}, {'ZCorr'; '.cc_zero'; '.cc_zero_z'}]; 
                    [figures{1,fig_count}, Correlations.ZCorr{1}, Correlations.ZCorrMat{1}, Correlations.ZCorrLag{1}] = ...
                        simple_motor_correlation(traces{all_ind}, Data.(s), NetworkOps.MotorCorr, block_row, 'TraceToAnalyze', 'zcorr', 'PlotFigure', PlotMotorFigure);
                    fig_count = fig_count + 1;
                    figureNames(fig_count) = {'Zcorr'};
                    PlotMotorFigure = 0; %Set to 0 once it has been plotted for the first time
                end   
            else
                disp('Skipping Loco Correlation Analysis...');
            end
            
            %% NOISE AND SIGNAL CORRELATIONS: find pairwise noise and signal correlations
            if ~contains(s,'pontaneous') && NetworkOps.NoiseCorr.ComputeNoiseCorr % if not a spontaneous block
                data = [data, {'NoiseCorr'; '.noiseCorr'; '.noise_z'}, {'NoiseCorr'; '.signalCorr'; '.signal_z'}]; 
                NetworkOps.NoiseCorr.fs = ExtractedData.StimInfo.fs; % frame rate
                for  n=1:length(loco_type) % 1 or 3 depending on loco filter used in extract data
                    if ~isempty(cell_rows{n}) && ~isempty(traces{n})
                    [Correlations.NoiseCorr{n}, blanks]...
                        = simple_noise_corr_network(Data.(s),loco_type{n},cell_rows{n},NetworkOps.NoiseCorr);
                    end
                end
            else
                disp('Skipping Noise and Signal Correlation Analysis...');
            end

             %% NOISE AND SIGNAL CORRELATIONS USING XCORR: find pairwise noise and signal correlations using xcorr (max across lags)
            if ~contains(s,'pontaneous') && NetworkOps.NoiseCorr.ComputeNoiseXCorr % if not a spontaneous block
                data = [data, {'NoiseXCorr'; '.noiseCorr'; '.noise_z'}, {'NoiseXCorr'; '.signalCorr'; '.signal_z'}]; 
                NetworkOps.NoiseCorr.fs = ExtractedData.StimInfo.fs; % frame rate
                for  n=1:length(loco_type) % 1 or 3 depending on loco filter used in extract data
                    if ~isempty(cell_rows{n}) && ~isempty(traces{n})
                    [Correlations.NoiseXCorr{n}, Xblanks]...
                        = simple_noise_Xcorr_network(Data.(s),loco_type{n},cell_rows{n},NetworkOps.NoiseCorr);
                    end
                end
            else
                disp('Skipping Noise and Signal XCorrelation Analysis...');
            end
            
            %% PLOT FIGURES
            for r = 1:size(data,2) % loop through correlation type to analyze
                figure_name = strrep([data{1,r} data{2,r}],'.','_');

                % Plot histograms of Correlations of Block
                % for All Comparisons, Only Non-Locomotion Trials, and Only Locomotion Trials
                % histograms of all noise, signal, significant noise, and significant signal correlations
                color = ['k', 'b', 'c'];
                locolegend = {'All', 'Loco','No Loco'};
                for m = 1:2 % loop through all and significant correlations
                    figures{1,fig_count} = figure;
                    figure(figures{1,fig_count});
                    if ismember(data(1,r),{'MotorCorr', 'XYCorr', 'ResidualsCorr', 'ZCorr'})
                        nCells=1;
                    else
                        nCells = length(loco_type);
                    end
                    
                    for n = 1:nCells
                        if ~isempty(cell_rows{n}) && ~isempty(traces{n})
                            data_temp = eval(string(strcat('Correlations.',data{1,r}, '{1,n}',data{2,r})));
                            if m == 2
                                % indices of significant correlated pairs
                                sig_indices = find(eval(string(strcat('Correlations.', data(1,r), '{1,n}', data(3,r)))) ==1);
                                % significant correlated pairs only
                                data_temp = data_temp(sig_indices); % only sig correlated pairs
                            end

                            if NetworkOps.histogram_line == 1   % if opted for figure with outlined histograms
                                histogram(data_temp,NetworkOps.bins, 'Normalization',NetworkOps.histogram_type,'DisplayStyle','stairs','LineWidth', 2, 'EdgeColor',color(n)); hold on
                            else % standard histogram
                                histogram(data_temp,NetworkOps.bins, 'Normalization',NetworkOps.histogram_type,'EdgeColor',color(n),'FaceColor', color(n)); hold on;
                            end
                            if m==1; title(strrep([data{1,r} data{2,r} ' histogram'], '_',' '));
                            else title(strrep([figure_name ' SigOnly_histogram'],'_',' '));end
                            xlabel('correlations'); ylabel(NetworkOps.histogram_type);
                        end
                    end
                    legend(locolegend(1:nCells));
                    if m==1; figureNames(fig_count) = {[figure_name '_Histogram']};
                    else figureNames(fig_count) = {[figure_name '_SigOnly_Histogram']}; end
                    fig_count = fig_count + 1;
                end

                % Plot Correlation as a Function of Distance for All Pairs in Block noise correlations by distance
                if ~ismember(data(1,r),{'MotorCorr', 'XYCorr', 'ZCorr', 'ResidualsCorr'})
                    for  n = 1:length(loco_type)
                        if ~isempty(cell_rows{n}) && ~isempty(traces{n})
                            data_temp = eval(string(strcat('Correlations.',data(1,r), '{1,n}',data(2,r))));
                            AA = distances{n};

                            AA(find(isnan(data_temp)))=nan;
                            
                            pairwise_distances{n} = AA(find(tril(AA,-1))); %Reshape and keep values below the diagonal, remove nans

                            sig_indices = find(eval(string(strcat('Correlations.', data(1,r), '{1,n}', data(3,r)))) ==1); % indices of significant correlated pairs
                            sig_data = data_temp(sig_indices);
                            indices_neg_correlations = find(sig_data<0);
                            indices_pos_correlations= find(sig_data>0);

                            Neg_Correlations = sig_data(indices_neg_correlations);
                            Pos_Correlations = sig_data(indices_pos_correlations);
                            pairwise_distances_sig{n} = pairwise_distances{n}(sig_indices) ; % distances of significantly correlated

                            distance_neg_noise = pairwise_distances_sig{n}(indices_neg_correlations);
                            distance_pos_noise = pairwise_distances_sig{n}(indices_pos_correlations);

                            %If no significant values to plot, continue (added by Maryse 4/3/2023)
                            if ~any(AA)
                                continue
                            end
                            
                            figures{1,fig_count} = figure;
                            figure(figures{1,fig_count});
                            scatter(pairwise_distances{n},data_temp,'k'); hold on;
                            scatter(distance_pos_noise, Pos_Correlations,'r'); hold on; % plot sig correlated in red
                            scatter(distance_neg_noise, Neg_Correlations,'c'); % plot sig correlated in bluem
                            figureNames(fig_count) = {[figure_name '_CorrbyDistance_' loco_type{n}]};

                            xlabel('distance (microns)'); ylabel('correlations');
                            title(strrep([figure_name ' by Distance '  loco_type{n}],'_',' '))
                            legend('All', 'Sig Positive', 'Sig Negative')
                            fig_count = fig_count + 1;
                        end
                    end
                end
            end

            %% SAVE RESULTS: Make list of c1 and c2 pairwise comparisons (both S2P and matched cell numbers)
            for n=1:length(loco_type)
                if ~isempty(cell_rows{n}) && ~isempty(traces{n})
                c1_list = []; sz = size(cell_rows{n},1); for mm = 1:sz; new_list = ones(sz-mm,1)*mm; c1_list=[c1_list; new_list]; end
                c2_list = []; sz = size(cell_rows{n},1); for mm = 1:sz; new_list = (mm+1:sz)'; c2_list=[c2_list; new_list]; end
                c1_list_s2p = nan(size(c1_list)) ;
                c2_list_s2p = nan(size(c2_list));
                c1_list_matchrow = nan(size(c1_list));
                c2_list_matchrow =  nan(size(c2_list));
                c1_list_redcell = nan(size(c1_list));
                c2_list_redcell = nan(size(c2_list));
                cr = cell_rows{n};

                for ff = 1:length(cr)
                    f1 = find(c1_list ==ff);
                    f2 = find(c2_list ==ff);
                    snum = Data.(s).CellList.CellNumber(cr(ff));
                    f11 = repmat(snum,[length(f1),1]);
                    f22 = repmat(snum,[length(f2),1]);
                    c1_list_s2p(f1) = f11;
                    c2_list_s2p(f2) = f22;
                    
                    redcell = Data.(s).CellList.RedCell(cr(ff));
                    f11 = repmat(redcell,[length(f1),1]);
                    f22 = repmat(redcell,[length(f2),1]);
                    c1_list_redcell(f1) = f11;
                    c2_list_redcell(f2) = f22;
                    
                    if ismember('MatchedRow', fieldnames(Data.(s).CellList))
                        mnum = Data.(s).CellList.MatchedRow(cr(ff));
                        f11 = repmat(mnum,[length(f1),1]);
                        f22 = repmat(mnum,[length(f2),1]);
                        c1_list_matchrow(f1) = f11 ;
                        c2_list_matchrow(f2) = f22 ;
                    else
                        c1_list_matchrow(f1) = nan ;
                        c2_list_matchrow(f2) = nan ;
                    end
                end
                Correlations.CellNumbers.S2P_Pairs{n} = [c1_list_s2p c2_list_s2p];
                Correlations.CellNumbers.Match_Pairs{n} = [c1_list_matchrow c2_list_matchrow];
                Correlations.CellNumbers.RedCell{n} = [c1_list_redcell c2_list_redcell];
                end
            end

            % Establish Save Paths for Figures and Data
            mouse = string(mouse_name); %(block.setup.mousename);
            stimName = string(s);
            block_filename = string(block_name);

            cd(NetworkOps.save_path)
            [~, ~, ~] = mkdir(mouse,stimName); cd(mouse);
            [~, ~, ~] = mkdir(stimName,block_filename);
            figure_path = fullfile(NetworkOps.save_path, mouse, stimName,block_filename);

            Correlations.BlockType = s;
            Correlations.CellRows = cell_rows;
            Correlations.BlockRow = block_row;

            for n=1:length(loco_type)
                if ~isempty(cell_rows{n}) && ~isempty(traces{n})
                    Correlations.CellNumbers.S2P{n} = Data.(s).CellList.CellNumber(cell_rows{n});
                    Correlations.CellNumbers.Matched{n} = Data.(s).CellList.MatchedRow(cell_rows{n});
                    Correlations.XCirc{n} = Data.(s).CellList.xcirc(cell_rows{n});
                    Correlations.YCirc{n} = Data.(s).CellList.ycirc(cell_rows{n});
                    Correlations.RedCell{n} = Data.(s).CellList.RedCell(cell_rows{n});
                end
            end

            Correlations.BlockName = convertCharsToStrings(block_name);
            Correlations.conv_factorX = Data.(s).BlockInfo.conv_factorX(block_row);
            Correlations.conv_factorY = Data.(s).BlockInfo.conv_factorY(block_row);
            Correlations.refImg = Data.(s).BlockInfo.refImg(block_row);
            Correlations.meanImg = Data.(s).BlockInfo.meanImgE(block_row);

            % Visualize Block
            if contains(s,'pontaneous') && NetworkOps.Xcorr.ComputeXcorr || (~contains(s,'pontaneous') && NetworkOps.Xcorr.ComputeXcorr  && ~NetworkOps.NoiseCorr.ComputeNoiseCorr) % if a spontaneous block
                for n = 1: length(loco_type)
                    if ~isempty(cell_rows{n}) && ~isempty(traces{n})
                    figures{1,fig_count+1} = figure; figures{1,fig_count+2} = figure;
                    [figures{1,fig_count+1},figures{1,fig_count+2}] = visualize_network(Correlations,NetworkOps.Visualize.selected_cell,n,NetworkOps.Visualize);
                    figureNames(fig_count+1) = {['Network_TraceMax' loco_type{n}]};
                    figureNames(fig_count+2) = {['Network_TraceZero' loco_type{n}]};
                    fig_count = fig_count +2;
                    end
                end
            elseif ~contains(s,'pontaneous') && NetworkOps.Xcorr.ComputeXcorr && NetworkOps.NoiseCorr.ComputeNoiseCorr
                for n = 1: length(loco_type)
                    if ~isempty(cell_rows{n}) && ~isempty(traces{n})
                    figures{1,fig_count+1} = figure; figures{1,fig_count+2} = figure; figures{1,fig_count+3} = figure; figures{1,fig_count+4} = figure;
                    [figures{1,fig_count+1},figures{1,fig_count+2},figures{1,fig_count+3},figures{1,fig_count+4}] = visualize_network(Correlations,NetworkOps.Visualize.selected_cell,n,NetworkOps.Visualize);
                    figureNames(fig_count+1) = {['Network_TraceMax' loco_type{n}]};
                    figureNames(fig_count+2) = {['Network_TraceZero' loco_type{n}]};
                    figureNames(fig_count+3) = {['Network_Noise' loco_type{n}]};
                    figureNames(fig_count+4) = {['Network_Signal' loco_type{n}]};
                    fig_count = fig_count +4;
                    end
                end
            end

            % Save Figures
            if NetworkOps.save_figures
                cd(figure_path)
                for f = 1:size(figures,2)
                    if isempty(figures{1,f}); continue; end
                    h_name = strjoin([block_filename '_Fig' num2str(f) '_' figureNames{f}]);
                    saveas(figures{1,f}, h_name, 'fig');
                end

                if isempty(Data.(s).TrialData(block_row).Loco)
                    h_name = strjoin([block_filename '_Locomotion']);
                    saveas(fig1, h_name, 'fig');
                end

                save('Correlations', 'Correlations', '-v7.3');
                close all;
                clear figures; 
            end

            % Build Table with Cell Data with all Relevant Correlations
            %loco_type = {'AllLoc', 'YesLoc','NoLoc'}; %Maryse commented so that loco_type can be defined by fieldnames instead (above)
            Cell_s2P_list = Data.(s).CellList.CellNumber(cell_rows{1});
            Matched_list = Data.(s).CellList.MatchedRow(cell_rows{1});
            RedCell_list = Data.(s).CellList.RedCell(cell_rows{1});

            for n=1:length(loco_type)
                if ~isempty(cell_rows{n}) && ~isempty(traces{n})
                    % create tables with correlations for easy plotting

                    length_corr = length(Correlations.CellNumbers.S2P_Pairs{n}(:,1));
                    %account for cases where there is only 1 responsive cell in FOV
                    if length_corr == 0
                        length_corr = 1;
                    end
                    MouseName = repmat(convertCharsToStrings(mouse_name),[length_corr,1]); %1
                    BlockName = repmat(convertCharsToStrings(block_name),[length_corr,1]); %2
                    FieldofView = repmat(convertCharsToStrings(FOV),[length_corr,1]); %3
                    StimName = repmat(convertCharsToStrings(s),[length_corr,1]);  %4
                    CellType = repmat(convertCharsToStrings(cell_type),[length_corr,1]); %5
                    LocoType = repmat(convertCharsToStrings(loco_type{n}),[length_corr,1]); %6
                    Cell1_s2P =  Correlations.CellNumbers.S2P_Pairs{n}(:,1); %7
                    Cell2_s2P = Correlations.CellNumbers.S2P_Pairs{n}(:,2); %9
                    Cell1_Match = Correlations.CellNumbers.Match_Pairs{n}(:,1);
                    Cell2_Match = Correlations.CellNumbers.Match_Pairs{n}(:,2);
                    Cell1_RedCell = Correlations.CellNumbers.RedCell{n}(:,1);
                    Cell2_RedCell = Correlations.CellNumbers.RedCell{n}(:,2);
                    Distance = pairwise_distances{n};

                    if NetworkOps.Xcorr.ComputeXcorr
                        TraceMaxCorrelations = Correlations.XCorrTraces{1,n}.cc_max;
                        TraceMaxCorrelationsLag = Correlations.XCorrTraces{1,n}.lag_max;
                        TraceMinCorrelations = Correlations.XCorrTraces{1,n}.cc_min;
                        TraceMinCorrelationsLag = Correlations.XCorrTraces{1,n}.lag_min;
                        TraceZeroCorrelations = Correlations.XCorrTraces{1,n}.cc_zero;
                        TraceMaxZTest = Correlations.XCorrTraces{1,n}.cc_min_z;
                        TraceMinZTest = Correlations.XCorrTraces{1,n}.cc_min_z;
                        TraceZeroZTest = Correlations.XCorrTraces{1,n}.cc_zero_z;
                    else
                        TraceMaxCorrelations = nan(length_corr,1);
                        TraceMaxCorrelationsLag = nan(length_corr,1);
                        TraceMinCorrelations = nan(length_corr,1);
                        TraceMinCorrelationsLag = nan(length_corr,1);
                        TraceZeroCorrelations = nan(length_corr,1);
                        TraceMaxZTest = nan(length_corr,1);
                        TraceMinZTest = nan(length_corr,1);
                        TraceZeroZTest = nan(length_corr,1);
                    end
                    
                    if ~contains(s,'pontaneous') && NetworkOps.NoiseCorr.ComputeNoiseCorr 
                        NoiseCorrelations = Correlations.NoiseCorr{1,n}.noiseCorr;
                        SignalCorrelations = Correlations.NoiseCorr{1,n}.signalCorr;
                        NoiseZTest = Correlations.NoiseCorr{1,n}.noise_z;
                        SignalZTest = Correlations.NoiseCorr{1,n}.signal_z;
                    else 
                        NoiseCorrelations = nan(length_corr,1);
                        SignalCorrelations = nan(length_corr,1);
                        NoiseZTest = nan(length_corr,1);
                        SignalZTest = nan(length_corr,1);
                    end

                     if ~contains(s,'pontaneous') && NetworkOps.NoiseCorr.ComputeNoiseXCorr % if compute NoiseCorr with XCorr
                        NoiseXCorrelations = Correlations.NoiseXCorr{1,n}.noiseCorr;
                        SignalXCorrelations = Correlations.NoiseXCorr{1,n}.signalCorr;
                        NoiseXZTest = Correlations.NoiseXCorr{1,n}.noise_z;
                        SignalXZTest = Correlations.NoiseXCorr{1,n}.signal_z;
                    else 
                        NoiseXCorrelations = nan(length_corr,1);
                        SignalXCorrelations = nan(length_corr,1);
                        NoiseXZTest = nan(length_corr,1);
                        SignalXZTest = nan(length_corr,1);
                    end

                    if length_corr == 1
                        %Set all empty variables to NaN if only 1 cell
                        [Cell1_s2P, Cell2_s2P, Cell1_Match, Cell2_Match, Cell1_RedCell, Cell2_RedCell, Distance,...
                        TraceMaxCorrelations, TraceMaxCorrelationsLag, TraceMinCorrelations, TraceMinCorrelationsLag, TraceZeroCorrelations, TraceMaxZTest, TraceMinZTest, TraceZeroZTest, ...
                        NoiseCorrelations, SignalCorrelations, NoiseZTest, SignalZTest] = deal(nan);
                    end
                    
                    Correlations_Matrix = table(MouseName, BlockName, FieldofView, StimName, CellType, LocoType,...
                        Cell1_s2P, Cell2_s2P, Cell1_Match, Cell2_Match, Cell1_RedCell, Cell2_RedCell, Distance,...
                        TraceMaxCorrelations, TraceMaxCorrelationsLag, TraceMinCorrelations,  TraceMinCorrelationsLag, TraceZeroCorrelations, TraceMaxZTest, TraceMinZTest, TraceZeroZTest, ...
                        NoiseCorrelations, SignalCorrelations, NoiseZTest, SignalZTest,NoiseXCorrelations, SignalXCorrelations, NoiseXZTest, SignalXZTest);

                    Correlations_Matrix_All = [Correlations_Matrix_All; Correlations_Matrix];
                end
            end

            %Save Motor, XY, and Z correlations in Motor table
            if NetworkOps.MotorCorr.ComputeMotorCorr
                if ~isempty(Data.(s).TrialData(block_row).Loco) || ~isempty(Data.(s).TrialData(block_row).Full_xoff) || ~isempty(Data.(s).TrialData(block_row).Full_zcorr)
                    if ~isempty(traces{all_ind})
                        nCells = size(traces{all_ind},1);

                        columnsToKeep = {'cc_zero', 'cc_zero_z', 'cc_zero_sign', 'cc_max', 'lag_max', 'cc_max_z', 'cc_min', 'lag_min', 'cc_min_z', 'cc_sign'};
                        renameColumns = {'Correlations', 'ZTest', 'Sign', 'MaxCorrelations', 'MaxCorrelationsLag', 'MaxZTest', 'MinCorrelations', 'MinCorrelationsLag', 'MinZTest', 'MaxMinSign'};

                        Motor_Matrix = table;
                        Motor_Matrix.MotorMouseName = repmat(convertCharsToStrings(mouse_name),[nCells,1]);
                        Motor_Matrix.MotorBlockName = repmat(convertCharsToStrings(block_name),[nCells,1]);
                        Motor_Matrix.MotorFieldofView = repmat(convertCharsToStrings(num2str(FOV)),[nCells,1]);
                        Motor_Matrix.MotorStimName = repmat(convertCharsToStrings(s),[nCells,1]);
                        Motor_Matrix.MotorCellType = repmat(convertCharsToStrings(cell_type),[nCells,1]); 
                        Motor_Matrix.Cell_s2P_list = Cell_s2P_list;
                        Motor_Matrix.Matched_list = Matched_list;
                        Motor_Matrix.RedCell_list = RedCell_list;

                        if hasLoco 
                            %Get percentage of time mouse was running and average running speed
                            loco_percent = (sum(activity > 0)/length(activity))*100;
                            loco_speed = mean(activity(activity > 0));

                            Motor_Matrix.MotorPercent = repmat(loco_percent, [nCells,1]);
                            Motor_Matrix.MotorSpeed = repmat(loco_speed, [nCells,1]);
                            for c = 1:length(columnsToKeep); Motor_Matrix.(strcat('Motor',renameColumns{c})) = Correlations.MotorCorr{1,1}.(columnsToKeep{c}); end
                        else
                            Motor_Matrix.MotorPercent = nan(nCells,1);
                            Motor_Matrix.MotorSpeed = nan(nCells,1);
                            for c = 1:length(columnsToKeep); Motor_Matrix.(strcat('Motor',renameColumns{c})) = nan(nCells,1); end
                        end

                        if hasXY
                            for c = 1:length(columnsToKeep); Motor_Matrix.(strcat('XY',renameColumns{c})) = Correlations.XYCorr{1,1}.(columnsToKeep{c}); end
                        else
                            for c = 1:length(columnsToKeep); Motor_Matrix.(strcat('XY',renameColumns{c})) = nan(nCells,1); end
                        end

                        if hasXY && hasLoco
                            for c = 1:length(columnsToKeep); Motor_Matrix.(strcat('Residuals',renameColumns{c})) = Correlations.ResidualsCorr{1,1}.(columnsToKeep{c}); end
                        else
                            for c = 1:length(columnsToKeep); Motor_Matrix.(strcat('Residuals',renameColumns{c})) = nan(nCells,1); end
                        end

                        if hasZcorr
                            for c = 1:length(columnsToKeep); Motor_Matrix.(strcat('Z',renameColumns{c})) = Correlations.ZCorr{1,1}.(columnsToKeep{c}); end
                        else
                            for c = 1:length(columnsToKeep); Motor_Matrix.(strcat('Z',renameColumns{c})) = nan(nCells,1); end
                        end

                        MotorCorrelations_Matrix_All = [MotorCorrelations_Matrix_All; Motor_Matrix];
                    end
                end
            end

            % alter block name to be compatible as field name
            new_name = makeFieldNameFromBlock(block_name); %MET created function 6/25 (in Utils folder)
            Correlations_Summary.(mouse_name).(stimName).(new_name) = Correlations;
            close all; clear Correlations;
        end
    end

end

%Correlations.Correlation_Data.([s]) = Correlation_Data;
cd(NetworkOps.save_path)

% save inter-cell correlations table
writetable(Correlations_Matrix_All,char(strcat('Correlations_Table_', NetworkOps.stimTypes, '.csv')),'Delimiter',',');
save(char(strcat('Correlations_Table_', NetworkOps.stimTypes, '.mat')), 'Correlations_Matrix_All', '-v7.3');
Correlations_Summary.NetworkOps = NetworkOps;
save(char(strcat('Correlations_Summary_', NetworkOps.stimTypes)), 'Correlations_Summary', '-v7.3');

% TODO: Spontaneous analysis needs to be re-evaluated first
% if only spontaneous data processed, add the spontaneous analysis from ExtractedData (N transients, transient rate, etc) to motor table
% if isequal(NetworkOps.stimTypes, {'Spontaneous'})
%     MotorCorrelations_Matrix_All = add_spontaneous_data_to_network_analysis(ExtractedData, MotorCorrelations_Matrix_All);
% end

% save motor correlations table
if NetworkOps.MotorCorr.ComputeMotorCorr
    writetable(MotorCorrelations_Matrix_All,char(strcat('MotorCorrelations_Table_', NetworkOps.stimTypes,'.csv')),'Delimiter',',');
    save(char(strcat('MotorCorrelations_Table_',NetworkOps.stimTypes,'.mat')), 'MotorCorrelations_Matrix_All') ;
end
disp('done analyzing network...')

%% save summary figures for pairwise cell correlations
if NetworkOps.Xcorr.ComputeXcorr || NetworkOps.NoiseCorr.ComputeNoiseCorr
    cd(NetworkOps.save_path);
    [~, ~, ~] = mkdir(NetworkOps.save_path,char(strcat('Summary_PairwiseCells_', NetworkOps.stimTypes)));

    for cell = 1:length(NetworkOps.Summary.all_celltypes) % loop through cell types
        NetworkOps.Summary.celltype = NetworkOps.Summary.all_celltypes{cell};
        cd(NetworkOps.save_path);
        [~, ~, ~] = mkdir(char(strcat('Summary_PairwiseCells_', NetworkOps.stimTypes)), NetworkOps.Summary.celltype);
        cd(NetworkOps.save_path)
        summary_path = fullfile(NetworkOps.save_path, char(strcat('Summary_PairwiseCells_', NetworkOps.stimTypes)),NetworkOps.Summary.celltype);
        PlotAllCorrelations(Correlations_Matrix_All,summary_path,NetworkOps.Summary);
    end

    cd(NetworkOps.save_path);
    [~, ~, ~] = mkdir(char(strcat('Summary_PairwiseCells_', NetworkOps.stimTypes)), 'AllCells');
    summary_path = fullfile(NetworkOps.save_path, char(strcat('Summary_PairwiseCells_', NetworkOps.stimTypes)),'AllCells');
    NetworkOps.Summary.celltype = [];
    PlotAllCorrelations(Correlations_Matrix_All,summary_path,NetworkOps.Summary);
end

%% save summary figures for motor correlations
if NetworkOps.MotorCorr.ComputeMotorCorr
    cd(NetworkOps.save_path);
    [~, ~, ~] = mkdir(NetworkOps.save_path,char(strcat('Summary_Motor_', NetworkOps.stimTypes)));

    for cell = 1:length(NetworkOps.Summary.all_celltypes) % loop through cell types
        NetworkOps.Summary.celltype = NetworkOps.Summary.all_celltypes{cell};
        cd(NetworkOps.save_path);
        [~, ~, ~] = mkdir(char(strcat('Summary_Motor_', NetworkOps.stimTypes)), NetworkOps.Summary.celltype);
        cd(NetworkOps.save_path)
        summary_path = fullfile(NetworkOps.save_path, char(strcat('Summary_Motor_', NetworkOps.stimTypes)),NetworkOps.Summary.celltype);
        PlotAllMotorCorrelations(MotorCorrelations_Matrix_All,summary_path,NetworkOps.Summary);
    end

    cd(NetworkOps.save_path);
    [~, ~, ~] = mkdir(char(strcat('Summary_Motor_', NetworkOps.stimTypes,'_AllCells')));
    summary_path = fullfile(NetworkOps.save_path, char(strcat('Summary_Motor_', NetworkOps.stimTypes,'_AllCells')));
    NetworkOps.Summary.celltype = [];
    PlotAllMotorCorrelations(MotorCorrelations_Matrix_All,summary_path,NetworkOps.Summary);
end

disp('All done :)')