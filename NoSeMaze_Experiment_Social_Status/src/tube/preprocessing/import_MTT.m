%% Functions
function [match_matrix, IDs] = import_MTT(filename)
    
    MTT_table = readtable(filename);
    IDs = MTT_table.Var1;
    
    MTT_table = table2cell(MTT_table(:,2:end));
    match_matrix = zeros(size(MTT_table));
    match_matrix(contains(MTT_table,'W')) = 1;
end