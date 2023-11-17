%% ---- QUICK EDITS TO .INP AND .DAT FILES -----


% - oldFolder is the name of the folder where the old file is stored
% - newFolder is the name of the new folder where the new file is saved
% - oldFile is the .INP or .DAT file that needs to be edited
% - newFile is the new saved file. Code will make a copy of oldFilename
% and save it as a new file after editing
% - oldExpr is the expression to change
% - newExpr is the expression that will replace expr2look

% By Andrew J. Buggee
%%

function [] = edit_INP_DAT_files(oldFolder,newFolder,oldFile,newFile,oldExpr,newExpr)

% if new folder doesnt exist, create it. 


if isfolder(newFolder) == false
    
    [status,msg,msgID] = mkdir(newFolder);
    
    
end


if ischar(newExpr)==1
    
    % open old file to edit
    oldText = fileread([oldFolder,oldFile]);
    
    % find location of expression to change
    [~,endI] = regexp(oldText,oldExpr,'match'); % find the old expression location
    
    % replace old expression with new expression
    oldText(endI:(endI+length(oldExpr)-1)) = newExpr; % insert new expression
    
    % save the edited text as a new file in the new folder
    writematrix(oldText,[newFolder,newFile],'Delimiter',' ','FileType','text','QuoteStrings',0);
    
elseif iscell(newExpr)==true
    % if there are multiple expressions to edit, the oldExpr cell has to be
    % the same length as the newExpr cell
    
    
    if length(oldExpr)==length(newExpr)
        
        % open old file to edit
        oldText = fileread([oldFolder,oldFile]);
        
        for ii = 1:length(oldExpr)
            
        % find location of expression to change
        [~,endI] = regexp(oldText,oldExpr{ii},'match'); % find the old expression location
        
        % replace old expression with new expression
        
        % check out the following example that can add strings to lines
        % (https://www.mathworks.com/matlabcentral/answers/122189-how-to-edit-a-text-file-using-matlab)
        oldText(endI:(endI+length(oldExpr{ii})-1)) = newExpr{ii}; % insert new expression
        
        end
        
        % save the edited text as a new file in the new folder
        writematrix(oldText,[newFolder,newFile],'Delimiter',' ','FileType','text','QuoteStrings','none');
        
        
    elseif length(oldExpr) ~= length(newExpr)
        
        
        error('The number of old expressions doesnt equal the number of new expressions')
        
    end
    
    
    
end