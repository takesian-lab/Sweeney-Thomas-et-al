 %% Batch enhance videos to prepare for DeepLabCut analysis
% Authors: Polley Lab modified by Selorm Quarshie and Maryse Thomas

%% Parameters to modify

redo = 0; %0 to skip videos that have already been enhanced, 1 to redo

disp(['Save over previously enhanced videos = ' num2str(redo)])

%% Select path

disp('Select Tosca or Session folder');
ToscaPath = uigetdir(path, 'Select Tosca or Session folder');
disp(ToscaPath);

%Check whether folder is a Tosca root folder or Session folder
[~,name,~] = fileparts(ToscaPath);

if strcmp(name(1:5), 'Tosca')
    rootfolder = 'Tosca';
elseif strcmp(name(1:7), 'Session')
    rootfolder = 'Session';
else
    error('Selected path is not a Tosca or Session folder, try again.')
end

cd(ToscaPath)

%% Get list of Session folders to enhance

switch rootfolder
    case 'Tosca'
        D = struct2table(dir);
        D = D(D.isdir,:); %Only keep folders
        path = string(unique(D.folder));
        sessionList = sort_nat(string(D.name(contains(D.name, 'Session'))));
        sessionNames = sessionList; %For displaying progress messages below
        for s = 1:length(sessionList)
            sessionList(s) = strcat(path, '\', sessionList(s));
        end
        
    case 'Session'
        sessionList = string(pwd);
        sessionNames = string(name);
end

if isempty(sessionList)
    disp('No sessions found to process, check root folder')
    return
end

%% Make Brighter AVIs

savingDir = "enhanced videos"; %name for folder where videos will be saved - RECOMMENDED DO NOT CHANGE
addgain = 0; %arbitrary number (e.g. 40). Higher = brighter
framerate = 50;

for s = 1:length(sessionList)
    
    %cd to session folder
    disp(sessionNames(s));
    avifolder = sessionList(s);
    cd(avifolder);
    
    %Get list of avis in directory
    files = struct2table(dir('*.avi'));
    filenames = string(files.name);
    
    %Skip if there are no avis
    if isempty(files)
        continue;
    end
    
    %Check if enhanced videos already exist
    if ~redo
        if isfolder(savingDir)
            disp('Found enhanced videos folder')
            
            %cd to enhanced video folder
            cd(strcat(avifolder, '/', savingDir))
            
            %Get list of files that had previously been enhanced
            existing_files = struct2table(dir('*.avi'));                
            existing_filenames = string(existing_files.name);
            
            previouslyEnhanced = false(size(filenames));
            if ~isempty(existing_filenames)
                for f = 1:length(existing_filenames)
                    F = char(existing_filenames(f));
                    F = F(10:end);
                    ind = find(strcmp(F, filenames));
                    if ~isempty(ind)
                        previouslyEnhanced(ind) = 1;
                    end
                end
            end
            
            %Remove from list of files to process
            filenames(previouslyEnhanced) = [];
            
            %Skip if no files left to enhance
            if isempty(filenames)
                continue
            end
            
            %cd back to main folder
            cd(avifolder);
        end
    end

    %Make enhanced videos folder
    [~,~,~] = mkdir(savingDir);
    
    tic                                                                 %starts elapsed timer

    for i = 1:length(filenames)

        filename = filenames(i);
        v = fullfile(avifolder, filename);
        v_temp = fullfile(sprintf('temp_%s',filename));                 %Creates temporary video file
        copyfile(filename, v_temp);

        video = VideoReader(v_temp);                                    %loads temporary video file
        numframes = video.NumFrames;                                    %number of total frames for each video file

        writeavi = VideoWriter(sprintf('brighter_%s',filename)); 
        writeavi.FrameRate = framerate;                                 %set your frame rate here
        open(writeavi);                                                 %start writing video

            for ii = 1:numframes
                im = read(video,ii);
                im_gray = rgb2gray(im);                                 %change video to grayscale
                im_hist = histeq(im_gray);                              %add im_hist for additional contrast
                im_hist = im_hist + addgain;                            %add gain
                writeVideo(writeavi,im_hist);
            end 

        close(writeavi)                                                 %finishes avi video writing 
        temp_file = video.Name;

        clear video
        delete(fullfile(avifolder, temp_file))                          %deletes the temporary video file
        movefile(fullfile(avifolder, writeavi.Filename), savingDir, "f" )
        
        disp(strcat({'Saved '}, filename));

    end

    toc                                                                   %ends elapsed timer
end

disp('All done!')