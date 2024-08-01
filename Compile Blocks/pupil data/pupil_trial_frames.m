function  [framenum, frametime] = pupil_trial_frames(FrameData,start)




framenum = [];
frametime = [];
for i = 1:length(FrameData)
    framenum = [framenum; FrameData(i).frames];
    zero_time = FrameData(i).tframe(:) - start;
    frametime = [frametime; zero_time];
end

end