%% helper_make_all_ExtractedData_the_same_ActivityType
% Use this if you want to plot all of your ExtractedData cells together,
% regardless of activity type (i.e. mix activated, prolonged, and suppressed)

%% load ExtractedData
orig_ExtractedData = ExtractedData; %Store just in case

%% COMBINE ACTIVATED AND PROLONGED

Loco = ExtractedData.Ops.loco_filter;

for L = 1:length(Loco)

    %Edit CellData
    isProlonged = ExtractedData.CellData.(Loco{L}).ResponseType == "prolonged";
    ExtractedData.CellData.(Loco{L}).ResponseType(isProlonged) = "activated";
    
    %Edit summary
    isProlonged = ExtractedData.Summary.(Loco{L}).ResponseType == "prolonged";
    ExtractedData.Summary.(Loco{L}).ResponseType(isProlonged) = "activated";
end

%% MAKE ALL ACTIVATED

Loco = ExtractedData.Ops.loco_filter;

for L = 1:length(Loco)

    %Edit CellData
    isRF = ExtractedData.CellData.(Loco{L}).RF == 1;
    ExtractedData.CellData.(Loco{L}).RF_Type(isRF) = "excitatory";
    ExtractedData.CellData.(Loco{L}).ResponseType(isRF) = "activated";
    
    %Edit summary
    isRF = ExtractedData.Summary.(Loco{L}).RF == 1;
    ExtractedData.Summary.(Loco{L}).RF_Type(isRF) = "excitatory";
    ExtractedData.Summary.(Loco{L}).ResponseType(isRF) = "activated";
end

%% Plot stuff

PlotOps = struct;
PlotOps.Loco                 = 'All'; %All, Running, or NotRunning
PlotOps.ResponseType         = ''; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
PlotOps.RF_Type              = ''; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
PlotOps.IsResponsive         = 2; %1 for responsive (and reliable), 0 for not, 2 for both
PlotOps.FOV                  = ''; %'' -> no filtering, or add FOV name here to only keep that FOV (e.g. L1), put ~ in front to keep everything excep that FOV
PlotOps.Group                = ''; %'' -> no filtering, or add group name here to only keep that group, put ~ in front to keep everything excep that Group
PlotOps.sortbyGCAMP          = 0; %0 for groups, 1 for gcamp, 2 to combine
PlotOps.sortbyCondition      = 0; %0 to ignore, 1 for first part, 2 for second part [e.g. Passive_60dB]
PlotOps.sortbyRedCell        = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
PlotOps.plotCorrelations	 = 0; %0 = don't plot dependent variable correlations, 1 = plot
PlotOps.smooth_rasters       = 1; %0 = don't smooth, 1 = smooth
PlotOps.plotShadedErrorBars  = 0; %0 = don't plot error bars, 1 = do
PlotOps.FibPhotChannel       = 2; %1 = blue, 2 = green, 3 = red, 4 = all of them
PlotOps.SuppressOutput       = 0; %0 to print Ops to command line, 1 to suppress
PlotOps.ByMouse              = 1; %Sort by unique mice. Currently only utilized in FreqDisc (stimProtocol =7)

plot_simple_ExtractedData(ExtractedData, PlotOps);    
%plot_simple_ExtractedData_AllBehavior(ExtractedData, PlotOps)