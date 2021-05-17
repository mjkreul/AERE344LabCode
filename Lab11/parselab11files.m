%% AerE 344 Parse Lab 11 Data
% saves time by exporting read files as csv

%initializing variables
AOA4 = zeros(3750,4,100);
AOA8 = zeros(3750,4,100);
AOA12 = zeros(3750,4,100);
AOA16 = zeros(3750,4,100);

%getting file path names
for i =1:100
    filename(i) = strcat("B00", num2str(i,'%03.f'),".dat");
end

for i = 1:100
    fid = fopen(strcat('AOA4\',filename(i)), 'r');
    for j = 1:4
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
        AOA4(k,1,i) =  str2double(junk{1});
        AOA4(k,2,i) =  str2double(junk{2});
        AOA4(k,3,i) =  str2double(junk{3});
        AOA4(k,4,i) =  str2double(junk{4});
        
        
        % iterate
        ms = ms + 1;
        k = k + 1;
        line = fgetl(fid);
    end
    
    % close the file and keep running through them until there are none
    % left
    fclose(fid);
end

for i = 1:100
    fid = fopen(strcat('AOA8\',filename(i)), 'r');
    for j = 1:4
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
        AOA8(k,1,i) =  str2double(junk{1});
        AOA8(k,2,i) =  str2double(junk{2});
        AOA8(k,3,i) =  str2double(junk{3});
        AOA8(k,4,i) =  str2double(junk{4});
        
        
        % iterate
        ms = ms + 1;
        k = k + 1;
        line = fgetl(fid);
    end
    
    % close the file and keep running through them until there are none
    % left
    fclose(fid);
end
for i = 1:100
    fid = fopen(strcat('AOA12\',filename(i)), 'r');
    for j = 1:4
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
        AOA12(k,1,i) =  str2double(junk{1});
        AOA12(k,2,i) =  str2double(junk{2});
        AOA12(k,3,i) =  str2double(junk{3});
        AOA12(k,4,i) =  str2double(junk{4});
        
        
        % iterate
        ms = ms + 1;
        k = k + 1;
        line = fgetl(fid);
    end
    
    % close the file and keep running through them until there are none
    % left
    fclose(fid);
end
for i = 1:100
    fid = fopen(strcat('AOA16\',filename(i)), 'r');
    for j = 1:4
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
        AOA16(k,1,i) =  str2double(junk{1});
        AOA16(k,2,i) =  str2double(junk{2});
        AOA16(k,3,i) =  str2double(junk{3});
        AOA16(k,4,i) =  str2double(junk{4});
        
        
        % iterate
        ms = ms + 1;
        k = k + 1;
        line = fgetl(fid);
    end
    
    % close the file and keep running through them until there are none
    % left
    fclose(fid);
end

csvwrite('AOA4.csv',AOA4)
csvwrite('AOA8.csv',AOA8)
csvwrite('AOA12.csv',AOA12)
csvwrite('AOA16.csv',AOA16)