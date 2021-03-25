function [values, time, names] = parseTestFiles(numFiles, fileStruct)
%% Get number of lines in files
fid = fopen(fileStruct(1).name, 'r');
lines = 0;
while feof(fid) == 0
  tline = fgetl(fid);
  lines = lines+1;
end
fclose(fid);
%% Preallocating the arrays


% this will save the "names" of the files which are just the 
names = zeros(1, numFiles);
% values from each of the files
values = [];
time = zeros(lines - 5, numFiles);

%% Running through each file
for i = 1:numFiles
    % get the name from the input struct
    stringname = fileStruct(i).name;
    % set the name as the file name minus the "_iw.txt" or "_in.txt"
    names(1,i) = str2double(strtok(stringname,"_"));
    % open the file
    fid = fopen(fileStruct(i).name, 'r');
    
    % skip through the first 6 lines since they are all the same and start
    % on line 6
    for j = 1:6
        line = fgetl(fid);
    end
    
    % initialize the number of ms each file runs through
    ms = 0;
    % initialize the indexer
    k = 1;
    
    % run through this while there is still a line
    while line ~= -1
        % split up the lines so we can grab the 3rd occurance and place it
        % in the values array
        junk = split(line);
        values(k, i) = str2double(junk{3});
        
        % add time to the array
        time(k,i) = ms;
        
        % iterate
        ms = ms + 1;
        k = k + 1;
        line = fgetl(fid);
    end
    
    % close the file and keep running through them until there are none
    % left
    fclose(fid);
end