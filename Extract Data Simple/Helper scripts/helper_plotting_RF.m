%% helper_quick visualization of all RF:::

%% User options

%Load ExtractedData first

save_figures = 1;
sortbyGCAMP = 0;
save_path = pwd;
PercentOpaque = 0.3; %Set opacity of RF mask (1 = no mask)

%Script will create one folder for each combination of these elements
PLOTBINARY = [0 1]; %[0 1]
RESPONSETYPE = {'activated', 'prolonged', 'suppressed', 'excitatory'}; %, '', 'excitatory'};
RFTYPE = {''}; %{'mixed', 'excitatory', 'inhibitory'};

%% Loop through all filtering options
cd(save_path)

for b = 1:length(PLOTBINARY)
    for r = 1:length(RESPONSETYPE)
        for f = 1:length(RFTYPE)
            
            if save_figures
                folder_name = strcat('RF_', RESPONSETYPE{r}, '_', RFTYPE{f});
                if PLOTBINARY(b) == 1
                    folder_name = [folder_name '_binary'];
                end
                [~, ~, ~] = mkdir(save_path,folder_name);
                cd([save_path '/' folder_name])
            end
            
            plot_binary_RF = PLOTBINARY(b);

            %Options for filtering data
            Ops = struct;
            Ops.Loco                = 'All'; %All, Running, or NotRunning
            Ops.ResponseType        = RESPONSETYPE{r}; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
            Ops.RF_Type             = RFTYPE{f}; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
            Ops.IsResponsive        = 1; %1 for responsive (and reliable), 0 for not, 2 for both
            Ops.sortbyGCAMP         = sortbyGCAMP; %0 for groups, 1 for gcamp, 2 to combine

            %% Setup

            SubsetData = simple_subset_ExtractedData(ExtractedData, Ops); %This function filters the data based on the above ops

            Ops = SubsetData.Ops;
            StimInfo = SubsetData.StimInfo;
            Summary = SubsetData.Summary;
            CellDataByStim = SubsetData.CellDataByStim;
            StimAnalysis = SubsetData.StimAnalysis;
            Groups = Summary.Group;
            GroupList = unique(Groups);

            %Setup for RF Mapping
            RF_Map = ExtractedData.StimInfo.RF_Map;
            nanmat = nan(max(RF_Map));
            RF_ind = sub2ind(size(nanmat),RF_Map(:,1),RF_Map(:,2));

            %Continuous variables to sort RFs by
            variables = {'Sparseness', 'Binary_Sparseness', 'R2', 'ISI', 'dPrime'}; %'BF', 'BF_I', 'BWBF_I', 'BWInt'};

            %% Plot one tiled figure for each group

            for G = 1:length(GroupList)
                for v = 1:length(variables)
                    
                    if ~any(ismember(Summary.Properties.VariableNames, variables{v}))
                        continue
                    end

                    %Sort data based on variable
                    [sorted_var,order] = sort(Summary.([variables{v}]));

                    FigHandle = figure;
                    t = tiledlayout('flow');

                    title(t,strcat(GroupList(G), {' sorted by '} , variables{v}))

                    %Go through all cells in Summary
                    for c = 1:length(Summary.CellNumber)
                        %Only plot the ones corresponding to the current group
                        if strcmp(Summary.Group(order(c)),GroupList(G))

                            %Get RF data
                            PeakData = StimAnalysis.Response(order(c),:);
                            [vmat, IsRF] = deal(nanmat);
                            vmat(RF_ind) = PeakData;
                            IsRF(RF_ind) = CellDataByStim(order(c)).RF;
                            RFmask = ones(size(nanmat))*PercentOpaque;
                            RFmask(IsRF == 1) = 1;

                            if plot_binary_RF
                                vmat = IsRF;
                            end

                            %Plot
                            nexttile
                            h = imagesc(vmat);
                            set(h, 'AlphaData', RFmask)
                            axis off
                            colormap(gca, bluewhitered(256))
                            title(num2str(round(sorted_var(c),2)), 'FontSize', 8)

                            %drawnow
                        end
                    end % length c

                    %% Save figures

                    if save_figures
                        FigSavePath = strcat(pwd, '/', GroupList(G), '_sortedby_', variables{v});
                        fullscreen(FigHandle); % Make all plots fullscreen before saving
                        set(0, 'CurrentFigure', FigHandle);
                        savefig(strcat(FigSavePath, '.fig'));
                        saveas(FigHandle,strcat(FigSavePath, '.png'));
                        close(FigHandle)
                    end

                end % variables list
            end% group list
        end
    end
end

disp('Done plotting')