function [value] = getXMLAttribute(XMLstruct, attribute1, attribute2, attribute3, attribute4)
%Access variables from XML struct based on attribute name
%Written by Maryse Thomas, May 2023 for use with MATLAB readstruct

if nargin < 3
    attribute2 = [];
end

if nargin < 4
    attribute3 = [];
end

if nargin < 5
    attribute4 = [];
end

%Convert struct to table
XMLtable = struct2table(XMLstruct); 

%First level of attributes
row1 = strcmp(XMLtable.key, attribute1);
value = XMLtable.value(row1);

%Second level of attributes
if ~isempty(attribute2)
    IndexedValueTable = struct2table(XMLtable.IndexedValue{row1}); %Convert substruct to table
    
    %get column names
    IVcolumns = IndexedValueTable.Properties.VariableNames;
    
    if ~any(ismember(IVcolumns, 'description'))
        columnname = 'index';
    else
        columnname = 'description';
    end
    
    row2 = strcmp(IndexedValueTable.(columnname), attribute2);
    value = IndexedValueTable.value(row2);
end

%Third level of attributes
if ~isempty(attribute3)
    SubindexedValueTable = struct2table(XMLtable.SubindexedValues{row1}); %Convert substruct to table
    
    %get column names
    IVcolumns = SubindexedValueTable.Properties.VariableNames;
    
    if ~any(ismember(IVcolumns, 'description'))
        columnname = 'index';
    else
        columnname = 'description';
    end
    
    row2 = strcmp(SubindexedValueTable.(columnname), attribute3);
    
    try
       SubSubindexedValueTable = struct2table(SubindexedValueTable.SubindexedValue(row2));
    catch
       SubSubindexedValueTable = struct2table(SubindexedValueTable.SubindexedValue{row2});
    end
    
    if height(SubSubindexedValueTable) > 1
        row3 = strcmp(SubSubindexedValueTable.description, attribute4);
        value = SubSubindexedValueTable.value(row3);
    else
        value = SubSubindexedValueTable.value;
    end

end
