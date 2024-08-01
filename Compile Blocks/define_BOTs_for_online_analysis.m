function [block] = define_BOTs_for_online_analysis(block)
% This function uses PrairieView BOTs as ROIs
%
% Argument(s):
%   block (struct)
%
% Returns:
%   block (struct)
%
% Notes:
%
% TODO:
% Search 'TODO'

%% Skip this function if not performing online analysis

keep_red_ROIs = 0; %0 to delete all red (ch1) BOTs, 1 to keep

if ismissing(block.setup.suite2p_path) || ~isequal(block.setup.suite2p_path, block.setup.block_path)
    return
end

%Accomodate multiplane data
if isfield(block, 'MultiplaneData')
    multiplaneData = true;
    nPlanes = block.setup.XML.nPlanes;
else
    multiplaneData = false;
    nPlanes = 1;
end

if ismissing(string(block.setup.BOT_filename)) && ~multiplaneData
    %block.setup.BOT_filename is set to nan when there is more than one BOT
    %file (e.g. multiplane data or T-Series trials).
    error('New datatype found. Troubleshoot live BOT analyasis')
end

disp('Performing live BOT analysis...');

if ~keep_red_ROIs
    disp('Deleting red channel BOTs')
end

%% SINGLE PLANE DATA: Load BOT file and store ROIs in block like suite2p data
 
if nPlanes == 1 
    
    %Load BOT csv from path stored in block
    filename = strcat(block.setup.block_path, '/', block.setup.BOT_filename);
    BOT_data = csvread(filename,1,0);
    timestamp = BOT_data(:,1)-BOT_data(1,1);
    ROIs = BOT_data(:,2:end)';
    N = size(ROIs,1);
    
    %Load XML from path stored in block
    disp('Reading XML file')
    xml_filename = strcat(block.setup.block_path, '/', block.setup.XML.filename);
    X = readstruct(xml_filename, 'AttributeSuffix', '');
    Regions = struct2table(X.Sequence(1).PVBOTs.Region);

    %Blank images for visualization - if Widefield, load reference image
    if ~ismissing(block.setup.VR_filename)

image = refImage_forBOTs(block);
        block.img.meanImg = image;
        block.img.refImg = image;
        block.img.max_proj = image;
        block.img.meanImgE = image;
        block.img.Vcorr = image;
    else
        block.img.meanImg = zeros(512,512);
        block.img.refImg = zeros(512,512);
    block.img.max_proj = zeros(512,512);
    block.img.meanImgE = zeros(512,512);
    block.img.Vcorr = zeros(512,512);
    end

    %Cell and neuropil data
    block.cell_number = 1:N;
    block.stat = {};
    block.ops.badframes = zeros(1,length(timestamp));
    block.ops.xoff = zeros(1,length(timestamp));
    block.ops.yoff = zeros(1,length(timestamp));
    block.F = ROIs;
    block.Fneu = zeros(size(ROIs));
    block.spks = zeros(size(ROIs));
    block.F7 = ROIs; %Neuropil corrected trace
    block.df_f = (block.F7 - mean(block.F7, 2, 'omitnan'))./mean(block.F7, 2, 'omitnan'); %DF/F = (total-mean)/mean
    
    %Fill stat and redcell with BOT region info
    block.redcell = zeros(N,1);
    for r = 1:size(Regions,1)
         if ~ismissing(block.setup.VR_filename)
             % for widefield. the BOT is slightly cropped compared to the
             % actual image. Adding 20 pixels fixes it well enough.
        x = Regions.x(r)+20;
        y = Regions.y(r)+20;
         else
                x = Regions.x(r);
        y = Regions.y(r);
         end
        width = Regions.width(r);
        height = Regions.height(r);
        
        %generate an ellipse
        t = -pi:0.063:pi; %0.63 gives 100 points
        xcirc = x + (width./2)*cos(t);
        ycirc = y + (height./2)*sin(t);
        %figure; scatter(xcirc, ycirc)

        %store in stat
        block.stat{1,r}.med(2) = x;
        block.stat{1,r}.med(1) = y;
        block.stat{1,r}.ycirc = ycirc;
        block.stat{1,r}.xcirc = xcirc;
        block.stat{1,r}.ypix = [];
        block.stat{1,r}.xpix = [];
        
        %update redcell
        channel = Regions.channel(r);
        if channel == 1
            block.redcell(r) = 1;
        end
    end
    
    %Remove red channel regions (optional)
    if ~keep_red_ROIs
        block.cell_number(block.redcell==1) = [];
        block.F(block.redcell==1,:) = [];
        block.Fneu(block.redcell==1,:) = [];
        block.spks(block.redcell==1,:) = [];
        block.F7(block.redcell==1,:) = [];
        block.df_f(block.redcell==1,:) = [];
        block.stat(block.redcell==1) = [];
        block.redcell(block.redcell==1) = [];
    end
    
    %Compute distance between BOT regions
    block.distance_in_microns = find_cell_distance(block);

end

%% MULTIPLANE DATA: Load multiple BOT csvs and store ROIs in block

if nPlanes > 1 
    
    %cd to block path to load csvs
    cd(block.setup.block_path)
    filedir = dir;
    filenames = {filedir(:).name};
    BOT_files = filenames(endsWith(filenames, 'botData.csv'));

    %Load and concatenate data across BOTs, separating by plane
    disp(['Loading ' num2str(length(BOT_files)) ' BOT files'])
    BOT_data_per_plane = cell(1,nPlanes);
    for b = 1:length(BOT_files)
        temp_BOT_data = csvread(BOT_files{b}, 1,0);
        for n = 1:nPlanes
            if b == 1
                BOT_data_per_plane{n} = temp_BOT_data(n,:);
            else
                BOT_data_per_plane{n} = [BOT_data_per_plane{n}; temp_BOT_data(n,:)];
            end
        end
    end   

    %Load XML from path stored in block
    disp('Reading XML file')
    xml_filename = strcat(block.setup.block_path, '/', block.setup.XML.filename);
    X = readstruct(xml_filename, 'AttributeSuffix', '');
    Regions = struct2table(X.Sequence(1).PVBOTs.Region);
    
    for n = 0:nPlanes-1
        currentPlane = strcat('plane', num2str(n));
        ROIs = BOT_data_per_plane{n+1}(:,2:end)';
        N = size(ROIs,1);
        T = size(ROIs,2);
        timestamp = block.timestamp.(currentPlane);
        
        %Make sure ROI duration matches timestamp duration
        %They can be mismatched if the number of frames is not divisible by the number of planes (there will be extra rows in the very last BOT)
        %This is usually fixed by compile_blocks but since we are just previewing data now we will simply pad the final mismatched frames with NaNs
        if length(timestamp) ~= T
            newT = length(timestamp);
            frameDiff = newT - T;
            paddedFrames = nan(N,frameDiff);
            ROIs = [ROIs, paddedFrames];
        end
        
        %Images for visualization
        block.img.(currentPlane).meanImg = zeros(512,512);
        block.img.(currentPlane).refImg = zeros(512,512);
        block.img.(currentPlane).max_proj = zeros(512,512);
        block.img.(currentPlane).meanImgE = zeros(512,512);
        block.img.(currentPlane).Vcorr = zeros(512,512);

        %Cell and neuropil data
        block.cell_number.(currentPlane) = 1:N;
        block.stat.(currentPlane) = {};
        block.ops.badframes.(currentPlane) = zeros(1,length(timestamp));
        block.ops.xoff.(currentPlane) = zeros(1,length(timestamp));
        block.ops.yoff.(currentPlane) = zeros(1,length(timestamp));
        block.F.(currentPlane) = ROIs;
        block.Fneu.(currentPlane) = zeros(size(ROIs));
        block.spks.(currentPlane) = zeros(size(ROIs));
        block.F7.(currentPlane) = ROIs; %Neuropil corrected trace
        block.df_f.(currentPlane) = (block.F7.(currentPlane) - mean(block.F7.(currentPlane), 2, 'omitnan'))./mean(block.F7.(currentPlane), 2, 'omitnan'); %DF/F = (total-mean)/mean
        
        %Special for multiplane data
        block.ops.nplanes = nPlanes;
        
        %Fill stat and redcell with BOT region info [the BOTs will be the same in both planes]
        block.redcell.(currentPlane) = zeros(N,1);
        for r = 1:size(Regions,1)
            x = Regions.x(r);
            y = Regions.y(r);
            width = Regions.width(r);
            height = Regions.height(r);

            %generate an ellipse
            t = -pi:0.063:pi; %0.63 gives 100 points
            xcirc = x + width*cos(t);
            ycirc = y + height*sin(t);
            %figure; scatter(xcirc, ycirc)

            %store in stat
            block.stat.(currentPlane){1,r}.med(2) = x;
            block.stat.(currentPlane){1,r}.med(1) = y;
            block.stat.(currentPlane){1,r}.ycirc = ycirc;
            block.stat.(currentPlane){1,r}.xcirc = xcirc;
            block.stat.(currentPlane){1,r}.ypix = [];
            block.stat.(currentPlane){1,r}.xpix = [];
            
            %update redcell
            channel = Regions.channel(r);
            if channel == 1
                block.redcell.(currentPlane)(r) = 1;
            end
        end
        
        %Remove red channel regions (optional)
        if ~keep_red_ROIs
            block.cell_number.(currentPlane)(block.redcell==1) = [];
            block.F.(currentPlane)(block.redcell==1,:) = [];
            block.Fneu.(currentPlane)(block.redcell==1,:) = [];
            block.spks.(currentPlane)(block.redcell==1,:) = [];
            block.F7.(currentPlane)(block.redcell==1,:) = [];
            block.df_f.(currentPlane)(block.redcell==1,:) = [];
            block.stat.(currentPlane)(block.redcell==1) = [];
            block.redcell.(currentPlane)(block.redcell==1) = [];
        end
    
        %Compute distance between BOT regions
        block.distance_in_microns.(currentPlane) = find_cell_distance(block, 'Plane', n);
    end
end

end
