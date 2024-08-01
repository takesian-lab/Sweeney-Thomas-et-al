
% If getting the index of one data set (not divided up as in 20% randomly
% removed) set Per1 to one (keep 100% data) and use resulting output from
% tIDX1/tTrials1
% if creating a hold-out set of 20%, set Per1 to 0.2. tIDX1 is the hold out
% tIDX2 is the remaining. tTrials1 and 2 correspond respectively. 

% numTrial: number of total trials to parse out
% Per1    : percent of traces to leave out
% trigTime: times of events/trials
% t       : size of output index

function [tIDX1, tIDX2, tTrials1, tTrials2] = RandTrials_GLM(numTrial,Per1,trigTime,t)


    
    %which trials are picked for test/train
    numTest = round(Per1*numTrial); 
    drawNum = randsample(numTrial,numTrial,'false'); %each time draw another set without replacement
    tTrials1 = drawNum(1:numTest);
    
    if Per1<0.95
    tTrials2 = drawNum(numTest+1:end);
    else 
         tTrials2 = NaN;
    end
    
    %what imaging frames do those test/train/validation trials correspond to
   tIDX1=false(t,1);
    for j=1:numel(tTrials1)
        idx1 = trigTime(tTrials1(j),1);
        idx2 = trigTime(tTrials1(j),2);
        tIDX1(idx1:idx2)=true;
    end
    
      if Per1<0.95  
tIDX2=false(t,1);
    for j=1:numel(tTrials2)
        idx1 = trigTime(tTrials2(j),1);
        idx2 = trigTime(tTrials2(j),2);
       tIDX2(idx1:idx2)=true;
    end
      else
          tIDX2 = NaN;
      end
end