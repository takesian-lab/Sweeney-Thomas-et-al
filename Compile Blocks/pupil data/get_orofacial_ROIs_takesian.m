function ROIdata = get_orofacial_ROIs_takesian(avifiles, filename, save_ROI_figure, plot_figure)
%Adapted from the Polley lab - August 2022

%this script takes all of the videos in one run of one session, then
%calculates the absolute difference between every frame of every video,
%utilizing the roi from first video of data set, set by user 

%% Selecting the ROI-- manually for the first video, then apply to all

%Naming and reading the video file
K=VideoReader(avifiles{1});

% get the first ROI
% here you need to preload the facial_ROI.mat to reuse the same ROIs
% use h = imrect(gca); to generate new ROIs and drag
% and place the rectangle
% h=imrect(gca, crop1);

ROI_finalized = 0;
while ROI_finalized == 0
    figure
    im1 = (read(K,K.NumFrames));
    imshow(im1(:,:,1));
    figtitle = char(filename);
    figtitle = figtitle(9:end); %remove 'ROIData_'
    title(figtitle)
    h=imrect(gca);
    ROI_finalized = input('Draw ROI. Resize box as needed and hit enter when done.');
    %ROI_finalized = input('Save ROI? 1 for yes, 0 for no: '); %Give user a chance to correct mistakes
end

% get the coordinates
Crop=getPosition(h);
crop1(1,1)=Crop(1,1);
crop1(1,2)=Crop(1,2);
crop1(1,3)=Crop(1,3);
crop1(1,4)=Crop(1,4);

% % optionally get a second ROI
% % h1=imrect(gca, crop2);
% h1=imrect(gca);
%
% % get the coordinates
% Crop_new=getPosition(h1);
% crop2(1,1)=Crop_new(1,1);
% crop2(1,2)=Crop_new(1,2);
% crop2(1,3)=Crop_new(1,3);
% crop2(1,4)=Crop_new(1,4);

%% save a reference figure with ROI in Tosca folder

if save_ROI_figure
    fig1 = figure; hold on
    im1_gray = rgb2gray(im1);
    imshow(im1_gray);
    rectangle('Position',crop1,'EdgeColor','w','LineWidth',2,'LineStyle','--'); hold on; % add the roi
    title(figtitle)
    saveas(fig1,filename,'fig');
    saveas(fig1,filename,'jpg');
    close all
    % rectangle('Position',crop2,'EdgeColor','r','LineWidth',2,'LineStyle','--'); hold on; % add the roi
    % export_fig('ROI',  '-png')
    % export_fig('ROI_fullface',  '-png')
end

%% Applying the crop (crop1) to the rest of the videos and their frames
% This is the time-consuming part of the function

roi = cell(1,length(avifiles));

%parallelize reading the video files
parfor i = 1:length(avifiles)
% %     If you get a codec error, try running this try/catch in a regular for loop
%     try
        v = VideoReader(avifiles{i});
%     catch
%         error(['Problem with reading video file ' avifiles{i} '. Try deleting and reconverting'])
%     end
    numframes = v.NumFrames;
    for ii = 1:numframes
        im = read(v,ii);
        im_gray = rgb2gray(im);
        %cropping to the ROI manually created in the first run
        roi{i}(:,:,ii) = single(imcrop(im_gray, crop1)); %convert to single so we can use NaNs later
        %roi_anterior{i}(:,:,ii) = imcrop(im_gray, crop2);
    end
end

%% compute diff between frames without concatenating everything

mean_diff_face = cell(1,length(roi)); %initialize a space for mean diff

for i = 1:length(roi)
    
    %check if any pixels are saturated and replace with nan
    roi{i}(roi{i} == 255) = nan;
    roi{i}(roi{i} == 0) = nan;
        
    diff_face = abs(diff(roi{i}, 1, 3));
    if (i < length(roi))
        left = roi{i}(:,:,end); %last frame
        right = roi{i+1}(:,:,1); %first frame of next video
        inbetween_diff = abs(right - left);
        diff_face = cat(3,diff_face,inbetween_diff); %concat with all differences from before
    end
    mean_diff_face{i} = mean(squeeze(mean(diff_face,1)),1,'omitnan'); %taking the average of the concatonated differences
end

Mtrace_face = []; %initialize a space for mean diff for every file
for i = 1:length(roi)
    Mtrace_face = [Mtrace_face,mean_diff_face{i}];
end

%% Saving all of the calculated differences

Mtrace_face = [0 Mtrace_face]; %adding zero to the start of the array
Mtrace_face_split = cell(1,length(avifiles));
framesinAvi = [];
for i = 1:length(avifiles)
    framesinAvi = [framesinAvi size(roi{i},3)];
end
cumframesAvi = [0 cumsum(framesinAvi)];
for i = 1:length(avifiles)
    Mtrace_face_split{i} = Mtrace_face(cumframesAvi(i)+(1:framesinAvi(i)));
end

ROIdata = struct;
ROIdata.Trace = Mtrace_face;
ROIdata.TraceByTrial = Mtrace_face_split;
ROIdata.Crop = crop1;

%% Plot figures

if plot_figure
    
    %Plot example for methodology
    v = 1; %CHANGE # TO PICK RANDOM VIDEO
    f = 1; %CHANGE # TO PICK RANDOM FRAME
    
    im = read(VideoReader(avifiles{v}),f); 
    im_gray = rgb2gray(im);
    im_crop = single(imcrop(im_gray, crop1));
     
    figure; hold on
    colormap('gray')
    
    subplot(1,3,1)
    image(im_crop)
    title('Raw grayscale image')
    axis off
    
    subplot(1,3,2)
    image(roi{v}(:,:,f))
    title('Saturated pixels removed')
    axis off
    
    subplot(1,3,3)
    diff_face = abs(diff(roi{v}, 1, 3));
    image(diff_face(:,:,f))
    title('Diff from next frame')
    axis off
    
    %Plot concatenated trials
    figure; hold on
    title('Whisker pad')
    estimated_x = (1:length(Mtrace_face))*1/50;
    xlim([1 estimated_x(end)])
    xlabel('Time (s)')
    plot(estimated_x, Mtrace_face)
    vline(cumframesAvi+1)
    for i = 1:length(cumframesAvi)
        text(cumframesAvi(i), nanmedian(Mtrace_face)*1.5, num2str(i))
    end  
end