% GLM for Ca imaging data
% V3 
% By Carolyn Sweeney
% see Pillow Lab tutorial on apollo
% 


%% Get your ops file and make sure that the settings are correct
 message = sprintf('Select an GLM_Ops file, check settings, and hit space bar');
        uiwait(msgbox(message));
  
    
[GLM_OpsFile, GLM_OpsPath] = uigetfile;
cd(GLM_OpsPath)
edit(GLM_OpsFile)
    pause;
run(GLM_OpsFile)
SetupData.Ops = Ops;
 %% Load extracted data
for i=1:length(Ops.stimTypes)
    D1 = dir(fullfile(data_path,'*.mat'));
    N = {D1.name};
    s = Ops.stimTypes{1,i};
    s1 = strcat(s,'_'); 
    name = strcat('ExtractedData_', s1);
    X = contains(N,name);
    cd(data_path);
    load(N{X});
    Data.([s])=ExtractedData;
end
clear ExtractedData

% load the compiled_table
if Ops.useCorrCell
    cd(NetworkAnalysisPath)
    load('CompiledTable')
end


%% Get baseline data:

% figure out which rows correspond to each block
s = Ops.stimTypes{1};
blist = {Data.([s]).TrialData.Block}; % list of blocks
bs = Data.([s]).StimInfo.baseline; % baseline
as = Data.([s]).StimInfo.after; % time after sound onset 
fs = Data.([s]).StimInfo.fs; % frame rate
bsfr = bs*fs; % frames of baseline data
ntfilt = Ops.window*fs; % number of frames to shift the data for regressor matrix
RidgeData  = table;% store the data from Ridge

%% loop through the blocks
tic
for bl = 1:length(blist)
    blist{bl}
    RidgeData.BlockName(bl) = blist(bl);
    RidgeData.Group(bl) = Data.([s]).BlockInfo.Group(bl);
    RidgeData.MouseID(bl) = Data.([s]).BlockInfo.MouseID(bl);
    t = Data.([s]).TrialData(bl).Timestamp';
    [D , params, ztraces, ResidTraces] = gSetup(bl,Data,s,blist,Ops,bsfr,t); % stimulus data/times/indicies
    
        
    if isempty(params.cellnum) % if there are no active cells (or other cells that meet experimental criteria)
        continue
    else
    
     for BR = 1:Ops.numBoots  % todo: do this 100-500x ==================================

         %-----------------------------------------
         % set test and training data sets aside:
         %-----------------------------------------
         if Ops.keepContig
             [testIdx, trainIdx1] =  GLM_mContigTr(t,FracTest); 
         else
             [testIdx, trainIdx1,testTrials, trainTrialsAll] = RandTrials_GLM(size(params.trigTime,1),...
                 Ops.FracTest,params.trigTime, size(Data.([s]).TrialData(bl).Timestamp',1));
         end
         
         if Ops.PredictMotor
             ztraces1 = ztraces';
         else
             ztraces1 = ztraces(params.cellIDX,:); % Ca traces for cells to analyze in GLM
         end
        
        if Ops.useAllCells | Ops.Residuals
              ztraces2 = ztraces; % Ca traces to be used as regressors
              sz = size(ztraces2);
        else
              sz = size(ztraces1);
        end
  

        
        %-----------------------------------------
        % make regressor matrix:
        % if you selected Ops.activeCells, they are not included. Rather,
        % they are added in during validation/testing because each
        % regressorMat will be different for each neuron. 
        %-----------------------------------------
       [XmatOnes, mask, maskLabel,params] = GLM_makeRegMat(params,Ops,ntfilt,sz,t);
        RidgeData.params(bl) = params; % save the parameters for the block
       
       % -----------------------------------------
       % make design mat for test data: we will make the training one below
       % with the validation. 
       % -----------------------------------------
        Xtest = XmatOnes(find(testIdx),:);
        
        
        %% ===  Ridge regression======================
        for ii = 1:size(ztraces1,1) % loop through neurons in the block.
            % ----------------------------------------------
            % Training/validation
            % ----------------------------------------------
            for cv = 1:Ops.RepeatFits % 10x validation - number of tests set in GLM_simple_eventOps
                if Ops.keepContig
                    [valIdx, trainIdx2] =  GLM_mContigTr(t,FracTest);
                else
                    [valIdx, trainIdx2, valTrials, trainTrials] = RandTrials_GLM(numel(trainTrialsAll),Ops.FracTest,...
                        params.trigTime(trainTrialsAll,:),size(t,1));
                end
                
                if Ops.useAllCells 
                    otherZtrace = ztraces2';
                    otherZtrace(:,ii) = []; % remove current neuron from the regressors
                    Xtrain = [XmatOnes(trainIdx2,:),otherZtrace(trainIdx2,:)];
                    Xval  = [XmatOnes(valIdx,:),otherZtrace(valIdx,:)];
                    if Ops.Residuals
                        otherZtrace = ResidTraces';
                        otherZtrace(:,ii) = []; % remove current neuron from the regressors
                        Xtrain = [Xtrain,otherZtrace(trainIdx2,:)];
                        Xval  = [Xval,otherZtrace(valIdx,:)];
                    end
                elseif Ops.useCorrCell
                    [CorrTraces, params, mask, maskLabel] = GLM_CorrCellReg(Compiled_Table,blist{bl},params.cellnum(ii),ztraces,Ops.CorrType,Ops,...
                        Data.([s]),ResidTraces,params,mask,maskLabel);
                    if size(CorrTraces,1)<size(CorrTraces,2)
                        CorrTraces = CorrTraces';
                    end
                    Xtrain = [XmatOnes(trainIdx2,:),CorrTraces(trainIdx2,:)];
                    Xval  = [XmatOnes(valIdx,:),CorrTraces(valIdx,:)];
                else
                    Xtrain = XmatOnes(trainIdx2,:);
                    Xval = XmatOnes(valIdx,:);
                    if Ops.Residuals
                        otherZtrace = ResidTraces';
                        otherZtrace(:,ii) = []; % remove current neuron from the regressors
                        Xtrain = [Xtrain,otherZtrace(trainIdx2,:)];
                        Xval  = [Xval,otherZtrace(valIdx,:)];
                    end
                end
                
                
                
                 
                
                
                
                
                spstrain = ztraces1(ii,trainIdx2); %Ca data for training - this is a single neuron trace
                spsval = ztraces1(ii,valIdx); % Ca data for validation - single neuron trace
                XXtr = Xtrain'*Xtrain; % Precompute some quantities (X'X and X'*y) for training and val data
                XYtr = Xtrain'*spstrain';  % design mat training data * Ca trace training data
                
                Imat = eye(size(Xtrain,2)); % identity matrix of size of filter + const
                Imat(1,1) = 0; % don't apply penalty to constant coeff
                
                for L = 1:numel(Ops.testLambda) % different lambdas to test
                    % "The MAP estimate for the LG model parameters has a closed form, making it
                    % simple and fast to compute: w_hat = (X'*X + lambda*I)^{-1} * X^T*Y "
                    w = (XXtr+(Ops.testLambda(L)*Imat)) \ XYtr;
                    
                    tempfilt = w(2:end); % filter, w(1) = constant
                    predRidge = w(1) + Xval(:,2:end)*tempfilt; % constant from model + (model * validation matrix)
                    
                    tempCC =  corrcoef(predRidge,spsval);% correlation between prediction and val data
                    CC(cv,L) = tempCC(1,2); % we will use these to figure out best lambda value to use on hold-out data
                end %lambda loop
            end % cross val loop
            
            % -----------------------------------
            % to find the "best lambda" take the mean CC for each lambda that
            % we tested and find the highest one
            % -----------------------------------
            [bestlamCC, lamIdx]  = max(nanmean(CC,1));
            
            % take all training data and use best lambda to fit the model (put train and val back together):
%             if Ops.keepContig
%                     [valIdx, trainIdx2] =  GLM_mContigTr(t,FracTest);
%                 else
%             [trainIdxAll, ~, ~,~] = RandTrials_GLM(numel(trainTrialsAll),1,params.trigTime(trainTrialsAll,:),size(t,1));
%             end
            
   
            if Ops.useAllCells
                otherZtrace = ztraces2';
                otherZtrace(:,params.cellIDX(ii)) = [];
                XtrainAll = [XmatOnes(trainIdx1,:),otherZtrace(trainIdx1,:)];
                if Ops.Residuals % if using all cells AND residuals, add the residual traces onto the all cell traces made above
                    otherZtrace = ResidTraces';
                    otherZtrace(:,params.cellIDX(ii)) = [];
                    XtrainAll = [XtrainAll,otherZtrace(trainIdx1,:)];
                end
            elseif Ops.useCorrCell; % if you chose Ops.CorrResid, the residuals were already taken into account
                XtrainAll = [XmatOnes(trainIdx1,:),CorrTraces(trainIdx1,:)];
            else
                XtrainAll = XmatOnes(trainIdx1,:);
                if Ops.Residuals % if youre using all residual traces but not the Ca activity from the traces
                    otherZtrace = ResidTraces';
                    otherZtrace(:,params.cellIDX(ii)) = [];
                    XtrainAll = [XtrainAll,otherZtrace(trainIdx1,:)];
                end
                
            end
            
            spstrainAll = ztraces1(ii,trainIdx1); %Ca data for training - this is a single neuron trace
            XXtr = XtrainAll'*XtrainAll; % compute some quantities (X'X and X'*y) for All training data (train + val)
            XYtr = XtrainAll'*spstrainAll';  % design mat training data * Ca trace training data
            
            Imat = eye(size(XtrainAll,2)); % identity matrix of size of filter + const
            Imat(1,1) = 0; % don't apply penalty to constant coeff
            
            w = (XXtr+Ops.testLambda(lamIdx)*Imat) \ XYtr; %: w_hat = (X'*X + lambda*I)^{-1} * X^T*Y "
            
            AllBestLam(ii) = Ops.testLambda(lamIdx); % labmda to use in model
            AllBestFilt(:,ii) = w; % filter (constant = w(1))
            
            if Ops.useCorrCell
                CorrTraceTest{ii} = CorrTraces(testIdx,:);
            end
            
        end
        toc
        % store the data:
        RidgeData.TestTrace{bl,BR} = single(ztraces1(:,testIdx)); % Ca trace for testing
        RidgeData.XmatTest{bl,BR}  = Xtest; % regressor mat for testing
        RidgeData.Filter{bl,BR} = AllBestFilt; % final model for neuron
        RidgeData.Lambda{bl,BR} = AllBestLam; % "best" lambda from validation
        RidgeData.mask{bl,BR} = mask;
        RidgeData.maskLabel{bl,BR} = maskLabel;
        if Ops.useAllCells
             RidgeData.TestTraceAll{bl,BR} = single(ztraces2(:,testIdx)); % Ca trace for testing
        end

        if Ops.Residuals
            RidgeData.TestResidualAll{bl,BR} = single(ResidTraces(:,testIdx));
        end
        if Ops.useCorrCell
            RidgeData.TestCorrAll{bl,BR} =  CorrTraceTest;
            clear CorrTraceTest
        end
    
        clear Xtest  AllBestFilt AllBestLam  RidgeEst mask maskLabel 
        
        %     toc
    end % for checking if block is empty
  
    end % numboots
      clear ztraces ztraces1 
     if Ops.useAllCells
           clear ztraces2
     end
      if Ops.Residuals
           clear ResidTraces
     end
end %loop through blocks

%% how does the model perform?
% make a histogram of each cell's corrcoeff for predicted vs. measured
% using the hold out data and model
ccAll = [];

for bl =  1:length(RidgeData.BlockName)
 
    if isempty(RidgeData.TestTrace{bl})
        continue
    else
        for BR = 1:Ops.numBoots
          
            testTrace = RidgeData.TestTrace{bl,BR}; % Ca data to test
            XtestTemp = RidgeData.XmatTest{bl,BR}; % regressor mat for test;
            lambda = RidgeData.Lambda{bl,BR}; % lambda chosen from validation/training
            Ws =  RidgeData.Filter{bl,BR}; % filter from validation/training
            
            % loop through cells in FOV :
            for ii =  1:size(testTrace,1)
                
                if Ops.useAllCells == 1
                    otherTraces = RidgeData.TestTraceAll{bl}';
                    otherTraces(:,RidgeData.params(bl).cellIDX(ii)) = [];
                    Xtest = [XtestTemp,otherTraces];
                    if Ops.Residuals
                        otherTraces = RidgeData.TestResidualAll{bl}';
                        otherTraces(:,RidgeData.params(bl).cellIDX(ii)) = [];
                        Xtest = [Xtest,otherTraces];
                    end
                elseif Ops.useCorrCell
                    AllCTr =RidgeData.TestCorrAll{bl};
                    CorrTr = AllCTr{1,ii};
                    Xtest = [XtestTemp,CorrTr];
                else
                    Xtest = XtestTemp;
                    if Ops.Residuals
                        otherTraces = RidgeData.TestResidualAll{bl}';
                        otherTraces(:,RidgeData.params(bl).cellIDX(ii)) = [];
                        Xtest = [Xtest,otherTraces];
                    end
                end
                
                
                T = testTrace(ii,:);
                
                const = Ws(1,ii);
                mask = RidgeData.mask{bl,BR};
                maskLabel = RidgeData.maskLabel{bl,BR};
                start = 1;
                pred = 0;
                % put each Beta in iteratively
                for j = 1:size(mask,2)
                    maskend = mask(j)+start;
                    predtemp =  Xtest(:,start+1:maskend)*Ws(start+1:maskend,ii);
                    start = maskend;
                    pred = pred + predtemp;
                end
               
                pred = pred + const; % put the constant in...
                R = corrcoef(pred,T');
                ccFOV(ii) = R(1,2);
            end
            
            
            RidgeData.CorrCoeff{bl,BR} = ccFOV;
%             ccAll = [ccAll,ccFOV];
            clear ccFOV
        end
        %         figure;
        %         histogram(ccFOV);
        %         title(RidgeData.BlockName{bl})
        %         xlabel('CC')
        %         ylabel('number of cells')
        %         drawnow % print figure before loop done
        %         clear ccFOV  predtemp pred
    end
    
    
end
% figure;
% histogram(ccAll);
% xlabel('CC')
% ylabel('number of cells')
% title('All cells')

%% Plot traces with CC>0.25
% TODO: combine with section above.

for bl = 1:length(RidgeData.BlockName)
     if isempty(RidgeData.TestTrace{bl})
        continue
    else
    figure('units','normalized','outerposition',[0 0 1 1]);
    TL = tiledlayout('flow');
    title(TL,RidgeData.BlockName{bl});
    xlabel(TL,'frames');
    ylabel(TL,'zscore');
    ccblock = RidgeData.CorrCoeff{bl};
    % find neurons that meet a threshold;
    ThrF  = ccblock>Ops.threshFits;
    testTrace = RidgeData.TestTrace{bl}; % Ca data to test
    XtestTemp = RidgeData.XmatTest{bl}; % regressor mat for test;
    lambda = RidgeData.Lambda{bl}; % lambda chosen from validation
    Ws =  RidgeData.Filter{bl}; % filter from validation
    mask = RidgeData.mask{bl};
    maskLabel = RidgeData.maskLabel{bl};
    
    % loop through cells in FOV :
    for ii =  1:size(testTrace,1)
        if Ops.useAllCells == 1
            otherTraces = RidgeData.TestTraceAll{bl}';
            otherTraces(:,RidgeData.params(bl).cellIDX(ii)) = [];
            Xtest = [XtestTemp,otherTraces];
            if Ops.Residuals
                otherTraces = RidgeData.TestResidualAll{bl}';
                otherTraces(:,RidgeData.params(bl).cellIDX(ii)) = [];
                Xtest = [Xtest,otherTraces];
            end
        elseif Ops.useCorrCell
            AllCTr =RidgeData.TestCorrAll{bl};
            CorrTr = AllCTr{ii};
            Xtest = [XtestTemp,CorrTr];
        else
            Xtest = XtestTemp;
            if Ops.Residuals
                otherTraces = RidgeData.TestResidualAll{bl}';
                otherTraces(:,RidgeData.params(bl).cellIDX(ii)) = [];
                Xtest = [Xtest,otherTraces];
            end
        end
        
        T = testTrace(ii,:);
        
        const = Ws(1,ii);
        start = 1; % this will grab from the second column -- column 1 is the constant!
        pred = 0;
        % put each Beta in iteratively
        for j = 1:size(mask,2)
            maskend = mask(j)+start;
            predtemp =  Xtest(:,start+1:maskend)*Ws(start+1:maskend,ii);
            start = maskend;
            pred = pred + predtemp;
            predset(j,:) = predtemp;
        end
        pred = pred + const; % put the constant in...
        
        MaskTraces(ii,:,:) = predset;
        
        nexttile
        if  ThrF(ii)
            linecolor_1 = [0 0 1];
            linecolor_2 = [0 1 0];
        else
            linecolor_1 = [.5 .5 .5];
            linecolor_2 = [0.3010, 0.7450, 0.9330];
        end
        %           figure;
        plot(smooth(T,10),'Color',linecolor_1); hold on
        plot(smooth(pred,10),'Color' ,linecolor_2);
        title(num2str(RidgeData.params(bl).cellnum(ii)))
        drawnow
        
    end
    RidgeData.MaskTraces{bl} = MaskTraces;
    clear MaskTraces predset
     end
end



%% leave one out CDF math:: Work in Progress::::::::::
% sort Ridge data by cell type:
% for bl = 1:length(blist)
%     RidgeData.Group(bl) = Data.([s]).BlockInfo.Group(bl);
%     RidgeData.MouseID(bl) = Data.([s]).BlockInfo.MouseID(bl);
% end

Groups = unique(RidgeData.Group);
for G = 1:length(Groups)
    gIX = find(strcmp(RidgeData.Group,Groups(G)));
    RD = RidgeData(gIX,:);
GLM_leave_one_out(RD,0)
title(Groups(G))
end


%% save figs?
FolderName = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\GLM\Example Features Nov2023\AllResidualsInFOV_allStim\Figures\AllMatched Cells Sound vs all Betas';
% FolderName = save_path;
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = num2str(get(FigHandle, 'Number'));
    set(0, 'CurrentFigure', FigHandle);
    %   saveas(FigHandle, strcat(FigName, '.png'));
    saveas(FigHandle, fullfile(FolderName,strcat(FigName, '.fig'))); % specify the full path
    saveas(FigHandle, fullfile(FolderName,strcat(FigName, '.png'))); % specify the full path
end

close all
%%
FolderName = 'Z:\Carolyn\Fiber Photometry\Behavior\GLM_VIP mice';
cd(FolderName)
savename = string('RidgeRegression_FreqDiscVIP_001.mat');
save(savename, 'RidgeData', '-v7.3');

savename= string('Ops_FreqDiscVIP_001.mat');
save(savename, 'Ops', '-v7.3');




