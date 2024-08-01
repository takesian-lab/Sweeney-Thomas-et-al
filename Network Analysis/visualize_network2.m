function [fig1, fig2, fig3, fig4] = visualize_network2(Correlations,selected_cell,loco,Ops)
% Modified code from Maryse Thomas, Carolyn Sweeney and Anne Takesian

% This function previews the network data output from Network Analysis. 

% Argument(s): 
% Correlations - Correlations from Network_Analysis_Pipeline 
% selected_cell - cell chosen by user to visualize correlations 
% loco = 1 all, 2 no loco, 3 loco
% Returns:
% figures showing the correlation of all types available (trace @ max lag, trace @
% zero lag, noise correlation, signal correlation) with selected cell in
% given field of view (block)

%% Setup
if isempty(Ops)
    Ops.plot_graph = 1;
    Ops.sig_only = 0; % only plot significant correlations with selected cell
    Ops.use_matched_cell = 0; % block cells = 0, matched cells =1
    Ops.ROI = 'mask';
end
loco_type = {'all loco', 'no loco', 'loco'};

%Return empty if figure could not be made
fig1 = []; fig2 = []; fig3 = []; fig4 = [];

%% Show Correlations in FOV
if Ops.plot_graph
   
    % if spont block, less figures
    if strcmp(Correlations.BlockType, 'Spontaneous') || strcmp(Correlations.BlockType, 'spontaneous')
        n=2;
    else
        n=4;
    end

    % define cell numbers for blocks
    if Ops.use_matched_cell == 0
        cells = Correlations.CellNumbers.S2P{loco};
        pairs = Correlations.CellNumbers.S2P_Pairs{loco};
    else 
        cells = Correlations.CellNumbers.Matched{loco};
        pairs = Correlations.CellNumbers.Match_Pairs{loco};
    end

    for i=1:n % 1 = Trace Max, 2 = Trace Zero, 3 = Noise, 4 = Signal
        if i == 1
            if ~isfield(Correlations, 'XCorrTraces')
                continue;
            end
            CorrelationData = Correlations.XCorrTraces{1,loco}.cc_max;
            ZMatrix = Correlations.XCorrTraces{1,loco}.cc_max_z;
        elseif i == 2
            if ~isfield(Correlations, 'XCorrTraces')
                continue;
            end
            CorrelationData = Correlations.XCorrTraces{1,loco}.cc_zero;
            ZMatrix = Correlations.XCorrTraces{1,loco}.cc_zero_z;
        elseif i == 3
            if ~isfield(Correlations, 'NoiseCorr')
                continue;
            end
            CorrelationData = Correlations.NoiseCorr{1,loco}.noiseCorr;
            ZMatrix = Correlations.NoiseCorr{1,loco}.noise_z;
        elseif i == 4
            if ~isfield(Correlations, 'NoiseCorr')
                continue;
            end
            CorrelationData = Correlations.NoiseCorr{1,loco}.signalCorr;
            ZMatrix = Correlations.NoiseCorr{1,loco}.signal_z;
        end

        % Find indices of all correlations of selected cell
        exampleCell = find(cells == selected_cell); 
        exampleCell_ind = [find(pairs(:,1)==selected_cell);find(pairs(:,2)==selected_cell)]; 

        figure('units','normalized','outerposition',[0 0 1 1]) % initialize figure
        ax0=axes; hold on; 
%        for m = 1:2
            if i == 1
                fig1 = figure; hold on
            elseif i == 2
                fig2 = figure; hold on
            elseif i == 3
                fig3 = figure; hold on
            elseif i == 4
                fig4 = figure; hold on;
            end

        %    if m == 1 % left plot is reference image
                textColor = 'w';
                ax = gca;
                image(Correlations.refImg{1,1});
%             elseif m == 2 % right plot is mean image
%                 textColor = 'c'; 
%                 ax = gca;
%                 image(Correlations.meanImg{1,1});
%             end

            if Ops.sig_only
                if i == 1
                   title(['Significant Trace Max Correlations with ' num2str(selected_cell)])
                elseif i == 2
                    title(['Significant Trace Zero Correlations with ' num2str(selected_cell)])
                elseif i == 3
                     title(['Significant Noise Correlations with ' num2str(selected_cell)])
                elseif i == 4
                    title(['Significant Signal Correlations with ' num2str(selected_cell)])
                end
            else
                if i == 1
                    title(['Trace Max Correlations with ' num2str(selected_cell)])
                elseif i == 2
                    title(['Trace Zero Correlations with ' num2str(selected_cell)])
                elseif i == 3
                    title(['Noise Correlations with ' num2str(selected_cell)])
                elseif i == 4
                    title(['Signal Correlations with ' num2str(selected_cell)])
                end
            end

            axis square
            xlim([0 512])
            ylim([0 512])
            colormap('bone')
            set(gca,'YDir','reverse')
            set(gca,'YTick',[0:100:500])
            set(gca,'YTickLabel',round([0:100:500]*Correlations.conv_factorX))
            set(gca,'XTick',[0:100:500])
            set(gca,'XTickLabel',round([0:100:500]*Correlations.conv_factorX))
            xlabel('Microns')
            ylabel('Microns')
  

           maxCorrelation = ceil(max(CorrelationData(exampleCell_ind))*100);
           minCorrelation = ceil(min(CorrelationData(exampleCell_ind))*100);

           if strcmp(Ops.ROI,'mask')
               imask=Correlations.meanImg{1,1}*0;
               imask2=Correlations.meanImg{1,1}*0;
               alpha = zeros(size(Correlations.meanImg{1,1}));
               alpha2 = zeros(size(Correlations.meanImg{1,1}));
           end

            % Plot cell color as a function of correlation
            for a = 1:length(cells)
                pairCell = cells(a); 
                pairCell_ind = [find(pairs(:,1)==pairCell);find(pairs(:,2)==pairCell)];
                
                if strcmp(Ops.ROI,'circle')
                    xc = Correlations.XCirc{loco}(a);
                    xcirc = xc{1,1};
                    yc = Correlations.YCirc{loco}(a);
                    ycirc = yc{1,1};

                    if pairCell == selected_cell
                        C = [1,1,1]; % White
                        alpha_value = 1;
                        subplot(1,2,m);
                        plot(xcirc,ycirc,'Linewidth', 5, 'Color', C);
                        alpha(alpha_value);
                        text(max(xcirc),max(ycirc),['Starter' num2str(pairCell)], 'Color', textColor); %Label with cell number
                    else
                        ind_pair = intersect(exampleCell_ind, pairCell_ind);
                        correlation = ceil(abs(CorrelationData(ind_pair)*100));
                        if ~isnan(correlation)

                            if CorrelationData(ind_pair)>0
                                colours = flipud(hot(maxCorrelation)); %Make list of hot colors for + associations
                            else
                                colours = flipud(cool(minCorrelation)); %Make list of cool colors for - associations
                            end

                            if Ops.sig_only % if only plotting significant correlations
                                if ZMatrix(ind_pair) == 1 % if significant correlation
                                    C = colours(correlation,:);
                                else
                                    C = [0,0,0]; %black circles for non-significant
                                end
                            else
                                C = colours(correlation,:); % plot all correlations in color
                            end

                        else
                            C = [0,0,0]; % don't plot neurons that don't have noise correlations
                        end

                        subplot(1,2,m);
                        plot(xcirc,ycirc,'Linewidth', 3, 'Color', C);
                        text(max(xcirc),max(ycirc),num2str(pairCell), 'Color', textColor); %Label with cell number
                    end

                elseif strcmp(Ops.ROI,'mask')
                    xpoints = double(cell2mat(Correlations.XPix{loco}(a)));
                    ypoints = double(cell2mat(Correlations.YPix{loco}(a)));
               %     ax1=axes; %linkprop([ops.Ax ax1],'Position');
               %     ax2=axes; %linkprop([ops.Ax ax2],'Position');

                    %                     xc = Correlations.XCirc{loco}(a);
                    %                     xcirc = xc{1,1};
                    %                     yc = Correlations.YCirc{loco}(a);
                    %                     ycirc = yc{1,1};
                    

                    if pairCell == selected_cell
                        %C = [1,1,1]; % White
                        imask2((xpoints)*512+ypoints)=1; %Make mask white
                        alpha(imask2>0)=1; %White mask will be less transparent than the others
                      %  alpha(imask2==0)=0;
                     %   text(max(xpoints),max(ypoints),['Starter' num2str(pairCell)], 'Color', textColor); %Label with cell number
                    else
                        ind_pair = intersect(exampleCell_ind, pairCell_ind);
                        correlation = ceil(abs(CorrelationData(ind_pair))*100);
                        if ~isnan(correlation)

%                             if CorrelationData(ind_pair)>0
%                                 colours = flipud(hot(maxCorrelation)); %Make list of hot colors for + associations
%                             else
%                                 colours = flipud(cool(minCorrelation)); %Make list of cool colors for - associations
%                             end

                            if Ops.sig_only % if only plotting significant correlations
                                if ZMatrix(ind_pair) == 1 % if significant correlation
                                    imask((xpoints)*512+ypoints)= correlation;
                                else
                                    imask2((xpoints)*512+ypoints)= -1; %rand(1)*.7; %Make mask a random color %black circles for non-significant
                                end
                            else
                                imask((xpoints)*512+ypoints)= correlation; % plot all correlations in color
                            end

                        else
                            imask2((xpoints)*512+ypoints)= -1; % don't plot neurons that don't have noise correlations
                        end
                        alpha(imask>0)=0.6;
                        alpha(imask==0)=0;
                        alpha2(imask2>0) = 1;
                        alpha2(imask2==0)=0;


                      %  text(max(xpoints),max(ypoints),num2str(pairCell), 'Color', textColor); %Label with cell number
                    end
                end
            end
            ax1=axes; linkprop([ax0 ax1],'Position');
            ax2=axes; linkprop([ax0 ax2],'Position');
            im2 = imagesc(ax2,imask,'alphadata',alpha);
            im1 = imagesc(ax1,imask2,'alphadata',alpha2);

            % axis 1 for nan 
            h = colorbar(ax1,'TickLabels','none');
            set(h,'visible','off');

            FOV_color_map = load('CFColormap.mat');
            FOV_color_map = FOV_color_map.mymap;
            colormap(ax2,FOV_color_map);
            colorbar % my colormap
            colormap(ax1,"bone"); caxis(ax1,[-1 1]);
            
            if ~isnan(minCorrelation) && ~isnan(maxCorrelation) && maxCorrelation>minCorrelation
                caxis(ax2,[0 40]);
                %caxis(ax2,[minCorrelation maxCorrelation]);
            end

            linkprop([ax0 ax2],'Position');
            linkprop([ax0 ax1],'Position');
            set([ax1,ax2],'Position',[.17 .11 .685 .815]);
            ax1.Position = ax2.Position;
            ax2.Visible = 'off';
            ax1.Visible = 'off';
            %ax1.ColorbarVisible = 'off';
            axis(ax1,'square');
            axis(ax2,'square');
            xlim([0 512])
            ylim([0 512])
            set(gca,'YDir','reverse');


            A = Correlations.BlockName; newA = strrep(A,'_', ' ');
            if isequal(version('-release'),'2021a') || isequal(version('-release'),'2021b')
                sgtitle(strjoin([newA, ' ', loco_type{loco}]));

                if Ops.use_matched_cell
                    subtitle('matched cells')
                else
                    subtitle('block cells')
                end
            else
                suptitle(strjoin([newA, ' ', loco_type{loco}]));
                if Ops.use_matched_cell
                    subtitle('matched cells')
                else
                    subtitle('block cells')
                end
            end
        end
    end
end


