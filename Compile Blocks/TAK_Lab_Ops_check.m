function [MissingOps] = TAK_Lab_Ops_check(CompBlk_OpsFile,CompBlk_OpsPath)


cd(CompBlk_OpsPath)     %Go to file location
run(CompBlk_OpsFile)    %Run the Ops.m file once they do so
UserOps = who;

keep CompBlk_OpsPath CompBlk_OpsFile UserOps
Ops_TEMPLATE = 'Ops_CompBlk_TEMPLATE'
run(Ops_TEMPLATE)

TemplateOps = who;

% compare the two ops variable lists:
count = 1;
MissingOps = [];
for i = 1:length(TemplateOps)
    a = TemplateOps{i};
   if ~strcmp(a,'UserOps') & ~strcmp(a,'Ops_TEMPLATE') % this ends up in TemplateOps by default using the who function
    if ~any(strcmp(UserOps,a))
        MissingOps{count} = a;
        count = count+1;
    end
   end

  

end

%  if ~isempty(MissingOps)
%      display(MissingOps)
%        error('User Ops file does not match template. Missing variables listed in Missing Ops')
%  end
 
end