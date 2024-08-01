%% Pupillometry analysis after deep lab cut
% note to carolyn - you are going to have to go back and figure out the
% frame time/numbers. Kens analysis appears to make a index of frames, not
% an absolute value for frames. Also, what are teh tiems?! Check them with
% the fibphot data to see if they all make sense!!!!!


setup = block.setup;
tosca_path = setup.Tosca_path;
cd(tosca_path)
stim_protocol = setup.stim_protocol;
loco_filter = 1;

%----------------------------------------------------------
if stim_protocol ==10
    dB = block.parameters.variable1;
   
elseif stim_protocol ==3 || 15 % Fm sweep and ripple
    dB =  block.parameters.variable2;
    varstim = block.parameters.variable1;
end
%------------------
% when did the stimulus occur?
st = block.Sound_Time;

%blank trials
stimIDX = find(dB>0);       blankIDX = find(dB<=0);

% loco trials
locoTrY = find(block.active_trials==1);
locoTrN = find(block.active_trials==0);

% remove loco trials from stim -  loco_filter : 0= dont sort by motor activity, 1= only non motor, 2= only motor
if loco_filter ==1 LF = locoTrN; elseif loco_filter ==2  LF = locoTrY; else  LF = 1:length(st); end

stimIDX = intersect(stimIDX,LF);        stimTime = st(stimIDX);
blankIDX = intersect(blankIDX,LF);      blankTime = st(blankIDX);
%---------------------------------------------------------


%%
% Do you want to analyze all points in the DLC output. If no, you can get
% rid of some here.
excludepoints = 1;
columnIDX = 1; % this is teh eye center. Im not using it!

filedir = dir;
filenames = {filedir(:).name};
csvFiles = filenames(endsWith(filenames, '.csv')); %find all .csv files

for i = 1:length(csvFiles)
    PF = csvFiles{i};
    M = readmatrix(PF);
    xs = M(:,2:3:end);
ys = M(:,3:3:end);
frames = M(:,1);


if excludepoints ==1
    xs(:,columnIDX) = [];
    ys(:,columnIDX) = [];
end
    

for j = 1:length(frames)
    [xc,yc,R,a] = circfit(xs(j,:),ys(j,:));
    xcenter(j) = xc;
    yc(j) = yc;
    Radius(j) = R;
end
  PupilSize{i} = Radius;
  clear frames Radius xs ys M xc yx R a
end

%%
%pupil size stim trials vs blank trials

Pupilstim = PupilSize(stimIDX);
Pupilblank = PupilSize(blankIDX);

% % temporary mean/plot
% meanstim = mean(Pupilstim{:}(:,1:300),1);
% meanblanks = mean(Pupilblank{:}(:,1:300),1);


plot(PupilSize{1,1})


