% use corrcell for regressor in matrix
function [outputTraces, params, mask, maskLabel] = GLM_CorrCellReg(Compiled_Table,block,cellOI,ztraces,CorrType,Ops,ExtractedData,ResidTraces,params,mask,maskLabel)
% load Compiled_Table
% match the block numbers
% match the cell numbers
% find the significant ones
% average them together (both positive and negative correlations)
% make sure to set up parameters and event times for the variable.
% need to have a variable with ALL cells from block because some correlated
% cells might be excluded from the

%  these variables will be inputs to the function:
% cells = params.cellnum;
% cellOI = cells(ii);
% block = blist{bl};

% CorrType = 'Trace'; %'Noise' 'Signal', 'Trace'


% find rows that correspond to the block that you're modeling:
blockrow = find(strcmp(string(block),Compiled_Table.BlockName));
exd_row = find(strcmp(string(block),ExtractedData.CellList.Block));
tempExDat = ExtractedData.CellList(exd_row,:);


% find matched row that corresponds with cellOI within blocksrows:
c1 = find(Compiled_Table.Cell1_s2P(blockrow) == cellOI);
c2 = find(Compiled_Table.Cell2_s2P(blockrow) == cellOI);

if Ops.CorrResid & ~Ops.CorrMean
    outputTraces = zeros(length(exd_row),size(ztraces,2));
end

if isempty(c1)  & isempty(c2)& ~ Ops.CorrResid
    outputTraces = zeros(2,size(ztraces,2));
else
    
    %get matched row as identifier:
    if ~isempty(c1)
        matchedCell = Compiled_Table.Cell1_Match(blockrow(c1(1)));
    else
        matchedCell = Compiled_Table.Cell2_Match(blockrow(c2(1)));
    end
    
    % find rows that correspond to the stim type
    if ~strcmp(Ops.Corrstim, 'ExceptCurrent')
        stimrow = find(Compiled_Table.StimName == Ops.Corrstim);
    else
        stimrow = find(Compiled_Table.StimName ~= Ops.stimTypes);
    end
    
    % find cellOI in stimrow:
    c3 = find(Compiled_Table.Cell1_Match(stimrow) == matchedCell);
    c4 = find(Compiled_Table.Cell2_Match(stimrow) == matchedCell);
    
    crow = [stimrow(c3);stimrow(c4)];
    if ~isempty(crow)
        if strcmp(string(CorrType),string('Noise'))
            colName = 'SigPosOrNegNoiseCorr';
        elseif strcmp(string(CorrType),string('Trace'))
            colName = 'SigPosOrNegTraceMaxCorr';
        elseif strcmp(string(CorrType),string('Signal'))
            colName = 'SigPosOrNegSignalCorr';
        else
            error('corr type not found')
        end
        posRow = crow(find(Compiled_Table.([colName])(crow)>0));
        negRow = crow(find(Compiled_Table.([colName])(crow)<0));
        
        
        % posRow and negRow can sometimes represent the same cell more than
        % once (i.e. a NeuronA can have noisecorr with NeuronB across
        % multiple stim types). Since we are only interested in the
        % identity of that neuron, here we find duplicate cell numbers and
        % remove them from the data. Otherwise, these cells will be added
        % to the matrix mulitple times below and be overrepresented in the
        % regressors matrix.
        if ~isempty(posRow)
            for p = 1:length(posRow)
                if Compiled_Table.Cell1_Match(posRow(p)) == matchedCell
                    pcell = Compiled_Table.Cell2_Match(posRow(p)); % identify the other match number
                elseif Compiled_Table.Cell2_Match(posRow(p)) == matchedCell
                    pcell = Compiled_Table.Cell1_Match(posRow(p));
                end
                
                % check that 'pcell' is in the current block. If you are
                % using 'ExceptCurrent' as an option, there is a
                % possibility that you are finding a significant
                % correlation for a cell that is in teh FOV, but was not
                % matched in teh current block of interest. In this case,
                % we need to skip this cell or you will get an error.
                if ~isempty(tempExDat.CellOrder(find(tempExDat.MatchedRow == pcell)))
                    pzIDX(p) = tempExDat.CellOrder(find(tempExDat.MatchedRow == pcell));
                else
                    pzIDX(p) = NaN;
                end
            end
            posIDX = unique(pzIDX);
            chk = isnan(posIDX);
            posIDX(chk) = [];
        else
            posIDX = [];
        end
        
        if ~isempty(negRow)
            for p = 1:length(negRow)
                if Compiled_Table.Cell1_Match(negRow(p)) == matchedCell
                    pcell = Compiled_Table.Cell2_Match(negRow(p));
                elseif Compiled_Table.Cell2_Match(negRow(p)) == matchedCell
                    pcell = Compiled_Table.Cell1_Match(negRow(p));
                end
                
                if ~isempty(tempExDat.CellOrder(find(tempExDat.MatchedRow == pcell)))
                    nzIDX(p) = tempExDat.CellOrder(find(tempExDat.MatchedRow == pcell));
                else
                    nzIDX(p) = NaN;
                end
            end
            negIDX = unique(nzIDX);
            chk2 = isnan(negIDX);
            negIDX(chk2) = [];
        else
            negIDX = [];
        end
        
        
        
        
        posZtraces = [];
        negZtraces = [];
        % now we need to find the cell numbers that go with those rows:
        
        for p = 1:length(posIDX)
            if Ops.CorrResid
                posZtraces = [posZtraces; ResidTraces(posIDX(p),:)];
            else
                posZtraces = [posZtraces; ztraces(posIDX(p),:)];
            end
        end
        
        for p = 1:length(negIDX)
            if Ops.CorrResid
                negZtraces = [negZtraces; ResidTraces(negIDX(p),:)];
            else
                negZtraces = [negZtraces; ztraces(negIDX(p),:)];
            end
        end
        
        
        if Ops.CorrMean
            otr1 = mean(posZtraces,1,'omitnan');
            otr2 = mean(negZtraces,1,'omitnan');
            
        else
            otr1 = posZtraces;
            otr2 = negZtraces;
        end
    else
        otr1 = [];
        otr2 = [];
        warning('missing corr traces for positive and negative')
    end
    
    
    if isempty(otr1)
        otr1 = zeros(1,length(ztraces));
    end
    if isempty(otr2)
        otr2 = zeros(1,length(ztraces));
    end
    
    % set params.eventcode and params.event label - !!!! moved to GLM_make_Regmat and GLM_simple_set_Regs!!!!!!
    % params.eventLabel{1,length(params.eventLabel)+1}  = 'pos corr neuron';
    % params.eventLabel{1,length(params.eventLabel)+1}  = 'neg corr neuron';
    %
    % for p = 1:size(otr1,1)
    % params.eventCode(1,length(params.eventCode)+1) = 7;
    % mask(length(mask)+1) = 1;
    % maskLabel{1,length(maskLabel)+1} = 'pos corr neuron';
    % end
    % for p = 1:size(otr2,1)
    % params.eventCode(1,length(params.eventCode)+1) = 7;
    % mask(length(mask)+1) = 1;
    % maskLabel{1,length(maskLabel)+1} = 'neg corr neuron';
    % end
    
    if Ops.CorrResid
        d1 = size(otr1,1);
        d2 = size(otr2,1);
        outputTraces(1:d1,:) = otr1;
        outputTraces(d1+1:d1+d2,:) = otr2;
    else
        outputTraces = [otr1;otr2];
    end
end
end





