
function meanImage = refImage_forBOTs(block)
blockname = block.setup.block_name;
folder = convertStringsToChars(block.setup.block_path); %convertStringsToChars(setup.BOTpath);
cd(folder);
d = dir([folder '*/*.ome.tif']);%extract tiffs

if contains(d(1).name,'Ch1')
    chan = 'Ch1';
    if contains(d(end).name,'Ch2')
        error(strcat('Collected both Ch1 and Ch2 for widefield for block ',block.setup.block_name));
    end
elseif contains(d(1,1).name,'Ch2')
    chan = 'Ch2';
end


% grab image 1 to check if the data are single or multipage tiffs:
num_BOT = sprintf('%06d',1);
fig_name = strcat(block.setup.block_name, '_Cycle00001_', chan, '_', num_BOT, '.ome.tif'); %Concatenate strings in () horizontally as the fig name
fig_name = convertStringsToChars(fig_name);
Y = imfinfo(fig_name);
 if length(Y) >1
     for i = 1:20
         image(:,:) = imread(fig_name,i);
         Full_Tile_Matrix(:,:,i) = image;
     end
 else
 

     for k=1:20 % load first 20 images to average. Will need to be updated for multipage tiff (April 2024)
         num_BOT = sprintf('%06d',k);
         fig_name = strcat(block.setup.block_name, '_Cycle00001_', chan, '_', num_BOT, '.ome.tif'); %Concatenate strings in () horizontally as the fig name
         fig_name = convertStringsToChars(fig_name);
         image(:,:,:) = imread(fig_name);
         Full_Tile_Matrix(:,:,k) = image;
     end
end
meanImage = mat2gray(mean(Full_Tile_Matrix,3));


end

