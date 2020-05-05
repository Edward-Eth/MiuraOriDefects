% Function to extract the contents of .FOLD files to MATLAB variables.
% Created by Steve Grey 25/06/2018 - University of Bristol
% Email: steven.grey@bristol.ac.uk
% NOTE: I am aware this is not the most efficient way to do what this
% code does. If any user of this code makes improvements I would appriciate
% access to any improved versions if possible.
% Edited by Edward Etheridge for speed.

function [file,frame,vertices,edges,faces,faceOrders,edgeOrders] = FOLD_reader(filename,rundocfolder)

%% Open the file and read
fid = fopen(fullfile(rundocfolder,filename));
numchar = 'placeholder';

% Initialise the counter
ii = 0;

% Run through reading every line and storing them individually
while ~feof(fid)
    ii = ii + 1;
    raw{ii} = fgets(fid);
end

%% Save data from the file in MATLAB format

% Initialise counters
count = 0;
kk = 0;

% This tells us where we are in the string
start_reading = 0;
for ii = 1:length(raw)-1
    
    % Starting off looking at the variables to which we are going to be
    % writing
    
    % This is for storing data which is a variable string
    if start_reading == 2 && isempty(strfind(raw{ii},'  ],')) ~= 0
        % Increment the counter for this particular variable. This will be
        % reset to 0 when we get to the end of this variable
        kk = kk+1;
        
        % if the data is in characters then this tells us that
        if length(strfind(raw{ii},'"')) == 2
            
            % Where are the quotes
            locs = strfind(raw{ii},'"');
            
            % Store whatever is between the quotes
            command_str = [var_name '{kk,:}=raw{ii}(locs(1)+1:locs(2)-1);'];
            
            % Evaluate the command
            eval(command_str)
        else
            % If it is not characters then is must be numbers
            
            % Initialise logic and counter vars
            break_logic = 0;
            ii2 = 0;
            
            % Run through the string looking for the end (represented by a
            % ']')
            while break_logic == 0
                
                % Increment counter
                ii2 = ii2 + 1;
                
                % If it is the last value in the string then save this and
                % break out
                if raw{ii}(ii2) == ']'
                    break_logic = 1;
                    vector_str(ii2) = raw{ii}(ii2);
                else
                    % If not save this element of the string
                    vector_str(ii2) = raw{ii}(ii2);
                end
            end
            
            % Convert the string to numbers and save
            command_str = [var_name '{kk,:}=str2num(vector_str);'];
            
            % Evaluate the command
            eval(command_str)
        end
        
    end
    
    % If we have reached the end of the variable we don't want to read
    % anything more into it and the variable counter should be reset to 0
    if isempty(strfind(raw{ii},'  ],')) == 0
        start_reading = 99;
        kk = 0;
    end
    
    % If this line has a variable name in it then we want to find this
    if isempty(strfind(raw{ii},'":')) == 0
        
        % Initialise the counters and logical tools
        start_reading = 0;
        var_logic = 1;
        var_count = 0;
        num_count = 0;
        char_count = 0;
        
        % Run through the line in question element by element
        for jj = 1:length(raw{ii})
            
            % The quote represents the start of the variable name
            if raw{ii}(jj) == '"'
                
                % Log that we are looking at the variable name
                % (start_reading = 1)
                start_reading = start_reading + 1;
                
            end
            
            % If we are looking at the variable name and it is not the
            % quotation mark still
            if start_reading == 1 && raw{ii}(jj) ~= '"'
                
                % Underscores in the name represent a new layer of a
                % structure
                if raw{ii}(jj) == '_'
                    
                    % Keep track of what layer of structure we are in
                    var_logic = var_logic + 1;
                    
                    % Initialise the count of which letter we are on
                    var_count = 0;
                else
                    % If it is not an underscore then it must be a letter
                    % we want to save somewhere
                    
                    % First thing to look at is the layer of structure we
                    % are in
                    if var_logic == 1
                        % Then increment the counter
                        var_count = var_count + 1;
                        % And save the character to the variable name
                        % string
                        var1(var_count) = raw{ii}(jj);
                    end
                    if var_logic == 2
                        var_count = var_count + 1;
                        var2(var_count) = raw{ii}(jj);
                    end
                    if var_logic == 3
                        var_count = var_count + 1;
                        var3(var_count) = raw{ii}(jj);
                    end
                    
                    % Currently only goes to 3 layers deep. Add more if
                    % statements in the same format as above if required
                    if var_logic > 3
                        error('You have gone too deep... Only structure within a structure supported.')
                    end
                end
            end
            
            % Now we have found the second quotation mark so this means
            % that we have finished the variable name
            if start_reading == 2 && raw{ii}(jj) == '"'
                
                % Construct the structure variable name
                if var_logic == 1
                    var_name = var1;
                elseif var_logic == 2
                    var_name = [var1 '.' var2];
                elseif var_logic == 3
                    var_name = [var1 '.' var2 '.' var3];
                end
                
                % Currently only goes to 3 layers deep. Add more if
                % statements in the same format as above if required
                if var_logic > 3
                    error('You have gone too deep... Only structure within a structure supported.')
                end
            end
            
            % Now we are looking at the contents if they are in-line with
            % the variable name
            if start_reading == 2 && raw{ii}(jj) ~= '"'
                
                % If it is numeric then save to a string
                if isempty(str2num(raw{ii}(jj))) ~= 1
                    num_count = num_count + 1;
                    num_char(num_count) = raw{ii}(jj); %#ok<*AGROW>
                end
                
                % The comma represents the end of the line
                if raw{ii}(jj) == ','
                    % Build the command to save the number to the right
                    % variable
                    command_str = [var_name '=' num_char ';'];
                    
                    % Evaluate
                    eval(command_str);
                    
                    % We are done with this variable now
                    start_reading = 99;
                end
                
            end
            
            % If there are this many quotation marks then it must be an
            % in-line character variable to save
            if start_reading == 3 && raw{ii}(jj) ~= '"'
                % Save the characters to a string ready for saving to the
                % correct variable
                char_count = char_count + 1;
                char_str(char_count) = raw{ii}(jj);
            end
            
            % This is the end of the contents of the variable so time to
            % execute
            if start_reading == 4 && raw{ii}(jj) == ','
                
                % Build the command
                command_str = [var_name '=' '''' char_str '''' ';'];
                
                % Evaluate
                eval(command_str);
                
                % We are finished with the writing to the variable
                start_reading = 99;
            end
        end
    end
    
    % Clear this each time incase the vector strings are of different length
    clear vector_str
    
    % If we have finished with a variable clear all the strings which could
    % be of a different length next time
    if start_reading > 3
        clear var1 var2 var3 command_str num_char var_name char_str
    end
end

% Now we're finished with everything clear anything 
clear break_logic char_count command_str count filename ii ii2 jj kk locs num_count raw start_reading var1 var2 var_count var_logic var_name

[~] =fclose(fid);
clear fid

% Check to see if all of the variables were in the file
C = who;

% There should be 7 if not find the one missing and write a 0 in it to
% return
if length(C) < 7
    
    logic_find = zeros(7,1);
    
    output_name{1} = 'file';
    output_name{2} = 'frame';
    output_name{3} = 'vertices';
    output_name{4} = 'edges';
    output_name{5} = 'faces';
    output_name{6} = 'faceOrders';
    output_name{7} = 'edgeOrders';
   
    for ii = 1:length(C)
       
        logic_find(1) = logic_find(1) + length(strfind(C{ii},output_name{1}));
        logic_find(2) = logic_find(2) + length(strfind(C{ii},output_name{2}));
        logic_find(3) = logic_find(3) + length(strfind(C{ii},output_name{3}));
        logic_find(4) = logic_find(4) + length(strfind(C{ii},output_name{4}));
        logic_find(5) = logic_find(5) + length(strfind(C{ii},output_name{5}));
        logic_find(6) = logic_find(6) + length(strfind(C{ii},output_name{6}));
        logic_find(7) = logic_find(7) + length(strfind(C{ii},output_name{7}));
        
    end
    locs = find(logic_find == 0);
    for jj = 1:length(find(logic_find == 0))
        var_name = output_name{locs(jj)};
        command_str = [var_name '=0;'];
        eval(command_str)
    end
end
