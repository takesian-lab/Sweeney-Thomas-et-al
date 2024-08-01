function plot_simple_ExtractedData_FreqDisc(ExtractedData, UserOps)
% Plot behavior data (loco, whisker, pupil, etc)

%% Subset data from ExtractedData (check parameters at the top of function for options!)

%Return if not behavior data
if ~ismember(ExtractedData.StimInfo.StimProtocol, [21]) %Repetative stim behavior (ABCXYZ)
    return;
end

%I RECOMMEND NOT CHANGING THESE AND USING USEROPS INSTEAD
if nargin < 2
    UserOps = struct;
    UserOps.Loco                = 'All'; %All, Running, or NotRunning
    UserOps.ResponseType        = ''; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
    UserOps.RF_Type             = ''; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
    UserOps.IsResponsive        = 1; %1 for responsive (and reliable), 0 for not, 2 for both
    UserOps.FOV                 = ''; %'' -> no filtering, or add FOV name here to only keep that FOV (e.g. L1)
    UserOps.sortbyGCAMP         = 0; %0 for groups, 1 for gcamp, 2 to combine
    UserOps.sortbyCondition     = 0; %0 to ignore, 1 for first part, 2 for second part [e.g. Passive_60dB]
    UserOps.sortbyRedCell       = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
    UserOps.smooth_rasters      = 1; %0 = don't smooth, 1 = smooth
    UserOps.FibPhotChannel      = 2; %1 = blue, 2 = green, 3 = red, 4 = all of them
    UserOps.SuppressOutput      = 0; %0 to print Ops to command line, 1 to suppress
    UserOps.ByMouse             = 1; %Sort by unique mice. Currently only utilized in FreqDisc (stimProtocol =7)
    UserOps.ByDay               = 1; % sort by days
    UserOps.plotShadedErrorBars = 1; %%0 = don't plot, 1 = plot shaded error bars
    UserOps.NumberDays         = 3; % how many days of data did we collect 
end

disp('Plotting ABC/XYZ data')

SubsetData = simple_subset_ExtractedData(ExtractedData, UserOps); %
Ops = SubsetData.SubsetOps;
BlockInfo = SubsetData.BlockInfo;
TrialData = SubsetData.TrialData;
StimInfo = SubsetData.StimInfo;
Summary = SubsetData.Summary;
GroupList = Summary.Group;
Groups = unique(GroupList);
MouseList = BlockInfo.MouseID;
Mice = unique(MouseList);
% BehaviorData = ExtractedData.BehaviorData; 
Days = unique(SubsetData.Summary.Date);


if isempty(TrialData)
    disp('No behavior data found to plot')
    return
end

%Make groups label without underscores for plots
GroupsLabel = Groups;
for g = 1:length(GroupsLabel)
    GroupsLabel{g} = regexprep(Groups{g},'_',' ','emptymatch');
end

%Useful for plotting
nBaselineFrames = StimInfo.nBaselineFrames;
xInFrames = nBaselineFrames + StimInfo.after_inFrames;
xInSeconds = xInFrames/StimInfo.fs;
xstart = -nBaselineFrames./StimInfo.fs;
graphx = linspace(xstart, xInSeconds, xInFrames); %changes trial length (x) in 180 frames to 6 seconds
xtick = 0:nBaselineFrames:xInFrames;
xticklabels = xtick./StimInfo.fs;

if  UserOps.FibPhotChannel < 4
    channel = UserOps.FibPhotChannel;
else
    channel = 1:3;
end

%% Add Day labels to ExtractedData
Day1 = zeros(height(ExtractedData.Summary.All.Group),1);
Day2 = zeros(height(ExtractedData.Summary.All.Group),1);
Day3 = zeros(height(ExtractedData.Summary.All.Group),1);
Day = zeros(height(ExtractedData.Summary.All.Group),1);
A_R = nan(height(ExtractedData.Summary.All.Group),1);
B_R = nan(height(ExtractedData.Summary.All.Group),1);
C1_R = nan(height(ExtractedData.Summary.All.Group),1);
C2_R = nan(height(ExtractedData.Summary.All.Group),1);
X_R = nan(height(ExtractedData.Summary.All.Group),1);
Y_R = nan(height(ExtractedData.Summary.All.Group),1);
Z1_R = nan(height(ExtractedData.Summary.All.Group),1);
Z2_R = nan(height(ExtractedData.Summary.All.Group),1);


A_peak = nan(height(ExtractedData.Summary.All.Group),1);
B_peak = nan(height(ExtractedData.Summary.All.Group),1);
C1_peak = nan(height(ExtractedData.Summary.All.Group),1);
C2_peak = nan(height(ExtractedData.Summary.All.Group),1);
X_peak = nan(height(ExtractedData.Summary.All.Group),1);
Y_peak = nan(height(ExtractedData.Summary.All.Group),1);
Z1_peak = nan(height(ExtractedData.Summary.All.Group),1);
Z2_peak = nan(height(ExtractedData.Summary.All.Group),1);

allMice = unique(ExtractedData.Summary.All.MouseID);
% allDays = unique(ExtractedData.Summary.All.Date);
for M = 1:length(allMice)
   ix1 = find(strcmp(ExtractedData.CellList.MouseID,allMice(M))); %rows that correspond to mice
   %find dates for this mouse.there should be 3. sort data by date. 
   mouseDays = sort(unique(ExtractedData.CellList.Date(ix1)));
  for j = 1:length(ix1)
      if strcmp(ExtractedData.CellList.Date(ix1(j)),mouseDays(1))
          Day1(ix1(j)) = 1;
          Day(ix1(j)) = 1;
      elseif strcmp(ExtractedData.CellList.Date(ix1(j)),mouseDays(2))
            Day2(ix1(j)) = 1;
            Day(ix1(j)) = 2;
      elseif strcmp(ExtractedData.CellList.Date(ix1(j)),mouseDays(3))
            Day3(ix1(j)) = 1;
            Day(ix1(j)) = 3;
      else
          error('day no here')
      end
      A_R(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).ReliabilityData(1).R;

  B_R(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).ReliabilityData(2).R;
   C1_R(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).ReliabilityData(3).R;
    C2_R(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).ReliabilityData(4).R;
     X_R(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).ReliabilityData(5).R;
      Y_R(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).ReliabilityData(6).R;
       Z1_R(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).ReliabilityData(7).R;
        Z2_R(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).ReliabilityData(8).R;
        
        
    A_peak(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).PeakData.Peak(1);
    B_peak(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).PeakData.Peak(2);
    C1_peak(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).PeakData.Peak(3);
     C2_peak(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).PeakData.Peak(4);
       X_peak(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).PeakData.Peak(5);
         Y_peak(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).PeakData.Peak(6);
           Z1_peak(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).PeakData.Peak(7);
             Z2_peak(ix1(j)) = ExtractedData.CellDataByStim.All(ix1(j)).PeakData.Peak(8);
             
     
        
        
  end
    
end

Loco = {'All','Running','NotRunning'};

for L=1:length(Loco)
ExtractedData.Summary.([Loco{L}]).Day1 = Day1;
ExtractedData.Summary.([Loco{L}]).Day2 = Day2;
ExtractedData.Summary.([Loco{L}]).Day3 = Day3;
ExtractedData.Summary.([Loco{L}]).Day = Day;

ExtractedData.CellData.([Loco{L}]).Day1 = Day1;
ExtractedData.CellData.([Loco{L}]).Day2 = Day2;
ExtractedData.CellData.([Loco{L}]).Day3 = Day3;
ExtractedData.Summary.([Loco{L}]).Day = Day;
end
ExtractedData.CellList.Day1 = Day1;
ExtractedData.CellList.Day2 = Day2;
ExtractedData.CellList.Day3 = Day3;
ExtractedData.Summary.([Loco{L}]).Day = Day;

ExtractedData.Summary.All.A_R = A_R;
ExtractedData.Summary.All.B_R = B_R;
ExtractedData.Summary.All.C1_R = C1_R;
ExtractedData.Summary.All.C2_R = C2_R;
ExtractedData.Summary.All.X_R = X_R;
ExtractedData.Summary.All.Y_R = Y_R;
ExtractedData.Summary.All.Z1_R = Z1_R;
ExtractedData.Summary.All.Z2_R = Z2_R;

ExtractedData.Summary.All.A_peak = A_peak;
ExtractedData.Summary.All.B_peak = B_peak;
ExtractedData.Summary.All.C1_peak = C1_peak;
ExtractedData.Summary.All.C2_peak = C2_peak;
ExtractedData.Summary.All.X_peak = X_peak;
ExtractedData.Summary.All.Y_peak = Y_R;
ExtractedData.Summary.All.Z1_peak = Z1_peak;
ExtractedData.Summary.All.Z2_peak = Z2_peak;






%% temporary plotting:
% average all stim across all mice/all days:
AllStim = [];
for j = 1:length(SubsetData.CellDataByStim)
    AllStim(:,:,j) = (SubsetData.CellDataByStim(j).StimTracesAveraged);
end
T= tiledlayout(2,4);

for ii = 1:size(AllStim,1)
%     if ii <5
%         xc = ii; yc = 1;
%     else
%         xc = ii-4; yc = 2;
%     end
    t1 = squeeze(AllStim(ii,:,:));
    t2 = smooth(mean(t1,2));
    t3 = smooth(std(t1,0,2)./size(t1,2));
    nexttile
    shadedErrorBar(graphx,t2,t3);
    title(SubsetData.StimInfo.V1(ii))
end
title(T,'Mean All Days')
% xlabel(t,'Distance (mm)')
% ylabel(t,'Size (mm)')


%% subset by days
if    UserOps.NumberDays == 3; 
Day1 = []; Day2 = []; Day3 = [];
numDays = 3;
else
    error('protocol changed. not 3 days long. update numDays')
end
for M = 1:length(Mice)
   ix1 = find(strcmp(SubsetData.CellList.MouseID,Mice(M)));
   dayUQ = sort(unique(SubsetData.CellList.Date(ix1)));
   for ii = 1:numDays
       keep = ix1(find(strcmp(SubsetData.CellList.Date(ix1),dayUQ(ii))));
       for jj = 1:length(keep)
           if ii ==1
               Day1 = cat(3,Day1,SubsetData.CellDataByStim(keep(jj)).StimTracesAveraged);
           elseif ii==2
               Day2 = cat(3,Day2,SubsetData.CellDataByStim(keep(jj)).StimTracesAveraged);
           else
               Day3 = cat(3,Day3,SubsetData.CellDataByStim(keep(jj)).StimTracesAveraged);
           end
       end
   end
    
end


%% plot output by days:
TT = tiledlayout(2,4)
for ii = 1:size(Day1,1)
   d1 = squeeze(Day1(ii,:,:));
   d2 = squeeze(Day2(ii,:,:));
   d3 = squeeze(Day3(ii,:,:));
   
   dm1 = smooth(mean(d1,2));
   dm2 = smooth(mean(d2,2));
   dm3 = smooth(mean(d3,2));
   
    ds1 = smooth(std(d1,0,2)./size(d1,2));
    ds2 = smooth(std(d2,0,2)./size(d2,2));
    ds3 = smooth(std(d3,0,2)./size(d3,2));
    
    nexttile
%     shadedErrorBar(1:size(d1,1),dm1,ds1);
plot(graphx,dm1)
    title(SubsetData.StimInfo.V1(ii))
    hold on
    plot(graphx,dm2)
    plot(graphx,dm3)
    legend ('day1','day2','day3')
    ylim([0 0.8])
%        nexttile
%       shadedErrorBar(1:size(d2,1),dm2,ds2);
%     title(strcat('Day2_',SubsetData.StimInfo.V1(ii)))
    
%        nexttile
%       shadedErrorBar(1:size(d3,1),dm3,ds3);
%     title(strcat('Day3_',SubsetData.StimInfo.V1(ii)))
end

%% make table of different measures across days:


%% Matching across Days:
matches = unique(ExtractedData.Summary.All.MatchedRow);
matchmat = zeros(length(matches),3);
matchedIX = matchmat;
matchedR1 = matchmat;
matchedR2 = matchmat;
matchedpeak = matchmat;
matchedRF = matchmat
for mr = 1:length(matches)
    m1 = find(ExtractedData.Summary.All.MatchedRow == matches(mr));
    for m2 = 1:length(m1)
        daycolumn = ExtractedData.Summary.All.Day(m1(m2));
        matchedIX(mr,daycolumn) = m1(m2);
        matchedR1(mr,daycolumn) = ExtractedData.Summary.All.R1(m1(m2));
        matchedR2(mr,daycolumn) = ExtractedData.Summary.All.R2(m1(m2));
        matchedpeak(mr, daycolumn)= ExtractedData.Summary.All.Peak(m1(m2));
        matchedRF(mr, daycolumn)= ExtractedData.Summary.All.RF(m1(m2));
    end
end


k = sum(matchedRF,2);
keep = find(k==3);
Y = matchedR1(keep,:);
for c = 1:length(Y)
    for cc = 1:3
        if Y(c,1) == 0
            Ynorm(c,cc) = Y(c,cc)./0.000001;
        else 
    Ynorm(c,cc) = Y(c,cc)./Y(c,1);
        end
    end
end


for mr = 1:length(Y)
    plot(1:3,Y(mr,:)); hold on
end
xlim([0 4])

%% R by stim by day for JMP



end