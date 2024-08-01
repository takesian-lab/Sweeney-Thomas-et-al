
function GLM_leave_one_out(RD,plot_all)

AllMaskCC = [];

for bl = 1:length(RD.BlockName)
    ccByMask = [];
    if isempty(RD.TestTrace{bl})
        continue
    else
%         leg = RidgeData.maskLabel{bl};
%         leg{1,length(leg)+1} = 'all';
%             leg = {'Motor','Sound', 'Ca Activity FOV', 'Residual FOV', 'All'};
 leg = {'Motor','Sound', 'CorrNeuron', 'All'};
        
        %-------------------------------------
        % error fix - The V2 data generated for April 2023 were missing an
        % event code. This will set that code to 2 (sound). this is
        % temporary so that We dont need to completely rerun for the paper.
        zerocode = find(RD.params(bl).eventCode ==0);
        if ~isempty(zerocode)
            RD.params(bl).eventCode(zerocode) =2;
        end
        %---------------------------------------
        
        
        uniqueMaskCode = unique(RD.params(bl).eventCode);

        
        testTrace = RD.TestTrace{bl}; % Ca data to test
        Ws =  RD.Filter{bl}; % filter from validation
        mask = RD.mask{bl};
        maskLabel = RD.maskLabel{bl};
        ccByMask = [];
        for j = 1:length(uniqueMaskCode)%(mask,2)
            maskset = 1:size(mask,2); % all the masks in dataset; remove last one because that is the whole model.
            codeCheck = uniqueMaskCode(j);%RidgeData.params(bl).eventCode(j); % which masks should be grouped together. If codes match
            matchCodes = find(codeCheck == RD.params(bl).eventCode);
            
      
            maskset(matchCodes) = []; % leave out this code (aka regressor)
            
            % loop through cells in FOV :
            for ii =  1:size(testTrace,1) %all cells
                T = testTrace(ii,:); %
                const = Ws(1,ii);
                leaveone = 0;
                for k = 1:length(maskset) % masks with one mask code removed
                    beta = squeeze(RD.MaskTraces{bl}(ii,maskset(k),:)); 
                    leaveone = leaveone +beta;
                end
           
                leaveone = leaveone +const;
                R = corrcoef(leaveone,T');
                ccByMasktemp(ii,j) = R(1,2);
            end
        end
        ccByMask = [ccByMask;ccByMasktemp];
        clear ccByMasktemp
    end
    adddim = size(ccByMask,2)+1;
    ccByMask(:,adddim) = RD.CorrCoeff{bl};
    % make cdf:
%     if plot_all
%     figure;
% %     for L = 1:size(ccByMask,2) 
% %         cdfplot(ccByMask(:,L)); hold on
% %         title(RidgeData.BlockName{bl})
% %         
% %     end
% %     legend('motor', ' sound', 'neurons','all')
% legend(leg)
%     end
%     clear ccByMask ccByMasktemp
AllMaskCC = [AllMaskCC;ccByMask];
end 
figure;
for L = 1:size(AllMaskCC,2) 
        cdfplot(AllMaskCC(:,L)); hold on 
end
% title('All data')
legend(leg)
    
end
