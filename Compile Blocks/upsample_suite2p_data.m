function block = upsample_suite2p_data(block,fs)
%% Upsample suite2p data to framerate fs
% y = resample(x,tx,fs) resamples x with timestamp tx at the uniform sample rate fs using a polyphase antialiasing filter

% Variables that will be upsampled:
% -timestamp
% -F
% -Fneu
% -spks
% -F_chan2
% -Fneu_chan2
% -zcorr
% -xoff
% -yoff

% Variables that will be recomputed based on upsampled data:
% -F7
% -df_f
% -F7_chan2
% -df_f_chan2
% -aligned_stim

%% Initialize variables

%Update framerate in block
fs_orig = block.setup.framerate;
block.setup.framerate_orig = fs_orig;
block.setup.framerate = fs;

%Check that original framerate is lower than fs
if fs <= block.setup.framerate_orig
    warning(['Original framerate is already higher than or equal to ' num2str(fs)])
end

%Check that fs is an integer
if mod(fs,1) ~= 0
    error('fs should be an integer')
end

%Accomodate multiplane data
if isfield(block, 'MultiplaneData')
    multiplaneData = true;
    nPlanes = block.setup.XML.nPlanes;
else
    multiplaneData = false;
    nPlanes = 1;
end


%% Perform upsampling

for n = 1:nPlanes
    
    if multiplaneData
        planeName = strcat('plane', num2str(n - 1));
        tx = block.timestamp.(planeName);
        [block.F.(planeName),ty] = resample_mat(block.F.(planeName),tx,fs);
        [block.Fneu.(planeName),~] = resample_mat(block.Fneu.(planeName),tx,fs);
        [block.spks.(planeName),~] = resample_mat(block.spks.(planeName),tx,fs);
        block.ops.xoff.(planeName) = resample_mat(block.ops.xoff.(planeName),tx,fs); %3/6/23 MET changed from: interp1(tx,single(block.ops.xoff.(planeName)),ty,'linear')'; %Linear interpolate xoff
        block.ops.yoff.(planeName) = resample_mat(block.ops.yoff.(planeName),tx,fs); %3/6/23 MET changed from: interp1(tx,single(block.ops.yoff.(planeName)),ty,'linear')'; %Linear interpolate yoff
        block.F7.(planeName) = block.F.(planeName) - block.setup.constant.neucoeff*block.Fneu.(planeName); %Neuropil corrected trace
        block.df_f.(planeName) = (block.F7.(planeName) - mean(block.F7.(planeName), 2, 'omitnan'))./mean(block.F7.(planeName), 2, 'omitnan'); %DF/F = (total-mean)/mean
        block.timestamp.(planeName) = ty;
        if isfield(block, 'F_chan2')
            [block.F_chan2.(planeName),~] = resample_mat(block.F_chan2.(planeName),tx,fs);
            [block.Fneu_chan2.(planeName),~] = resample_mat(block.Fneu_chan2.(planeName),tx,fs);
            block.F7_chan2.(planeName) = block.F_chan2.(planeName) - block.setup.constant.neucoeff*block.Fneu_chan2.(planeName); %Neuropil corrected trace
            block.df_f_chan2.(planeName) = (block.F7_chan2.(planeName) - mean(block.F7_chan2.(planeName), 2, 'omitnan'))./mean(block.F7_chan2.(planeName), 2, 'omitnan'); %DF/F = (total-mean)/mean
        end
        if isfield(block, 'zcorr')
            block.zcorr.(planeName) = resample_mat(block.ops.zcorr.(planeName),tx,fs); %3/6/23 MET changed from: interp1(tx,block.zcorr.(planeName)',ty,'linear')'; %Linear interpolate zcorr
        end
    else
        tx = block.timestamp;
        [block.F,ty] = resample_mat(block.F,tx,fs);
        [block.Fneu,~] = resample_mat(block.Fneu,tx,fs);
        [block.spks,~] = resample_mat(block.spks,tx,fs);
        block.ops.xoff = resample_mat(block.ops.xoff,tx,fs); %3/6/23 MET changed from: interp1(tx,single(block.ops.xoff),ty,'linear')'; %Linear interpolate xoff
        block.ops.yoff = resample_mat(block.ops.yoff,tx,fs); %3/6/23 MET changed from: interp1(tx,single(block.ops.yoff),ty,'linear')'; %Linear interpolate yoff
        block.F7 = block.F - block.setup.constant.neucoeff*block.Fneu; %Neuropil corrected trace
        block.df_f = (block.F7 - mean(block.F7, 2, 'omitnan'))./mean(block.F7, 2, 'omitnan'); %DF/F = (total-mean)/mean
        block.timestamp = ty;
        if isfield(block, 'F_chan2')
            [block.F_chan2,~] = resample_mat(block.F_chan2,tx,fs);
            [block.Fneu_chan2,~] = resample_mat(block.Fneu_chan2,tx,fs);
            block.F7_chan2 = block.F_chan2 - block.setup.constant.neucoeff*block.Fneu_chan2; %Neuropil corrected trace
            block.df_f_chan2 = (block.F7_chan2 - mean(block.F7_chan2, 2, 'omitnan'))./mean(block.F7_chan2, 2, 'omitnan'); %DF/F = (total-mean)/mean
        end
        if isfield(block, 'zcorr')
            block.zcorr = resample_mat(block.zcorr,tx,fs); %3/6/23 MET changed from: interp1(tx,block.zcorr',ty,'linear')'; %Linear interpolate zcorr
        end
    end 
end

%% Redo align to stim if needed

if isfield(block, 'aligned_stim')
   block = align_to_stim(block);
end

end

function [y,ty] = resample_mat(x,tx,fs)
%As far as I can tell, resample cannot be used on a matrix, so I will do that here

if isempty(x)
    y = x;
    ty = tx;
else
    x = double(x); %Convert suite2p variables to double before using resample

    for i = 1:size(x,1)
        if i == 1
            [y,ty] = resample(x(i,:),tx,fs);
            y = [y; nan(size(x,1)-1,length(ty))]; %Preallocate once we know how long ty is
        else
            y(i,:) = resample(x(i,:),tx,fs);
        end
    end

    y = single(y); %Convert back to single

    plot_figure = 0;
    if plot_figure
        figure

        cellrow = 1; 

        subplot(2,1,1); hold on
        %stem(tx,x(cellrow,:))
        plot(tx,x(cellrow,:))
        stairs(tx,x(cellrow,:),'--')
        title('Original signal')
        set(gca, 'FontSize', 12)

        subplot(2,1,2); hold on
        %stem(ty,y(cellrow,:))
        plot(ty,y(cellrow,:))
        stairs(ty,y(cellrow,:),'--')
        title(['Upsampled signal ' num2str(fs) 'Hz'])
        xlabel('Seconds')
        set(gca, 'FontSize', 12) 
    end
end

end

