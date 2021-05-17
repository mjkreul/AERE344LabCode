%% AerE344 Lab 7 Create matricies

% This script uses M Kreul's parseTestFiles.m used previously as a basis
% for writing .csv files containing only necessary data. Due to the amount
% of time needed to simply read the raw data, this script must only be run
% once in order to create .csv files of only the voltages that can be read
% and calculated from directly.
clear, clc


%% Preallocating matricies, naming values
distvals = ["2-5";"2-7";"2-9";"3-1";"3-3";"3-5";"3-7";"3-9";"4-1";"4-3";"4-5";"4-7";"4-9";"5-1";"5-3";"5-5"];
aoavals = ["0";"4";"8";"12"];
aoa0 = zeros(10000,16);
aoa4 = zeros(10000,16);
aoa8 = zeros(10000,16);
aoa12 = zeros(10000,16);

%% Reading the data
for m = 1:4
    for n = 1:16
        cur_aoa = aoavals(m);
        cur_dist = distvals(n);
        filename = strcat('AOA',cur_aoa,'/',cur_dist,'.txt');
        fid = fopen(filename, 'r');
%         disp(fid)
        
        % skip through the first 6 lines since they are all the same and start
        % on line 6
        for j = 1:6
            line = fgetl(fid);
        end
        
        % initialize the indexer
        k = 1;
        
        % run through this while there is still a line
        while line ~= -1
            % split up the lines so we can grab the 3rd occurance and place it
            % in the values array
            junk = split(line);
            if cur_aoa == '0'
                aoa0(k, n) = str2double(junk{3});
            elseif cur_aoa == '4'
                aoa4(k, n) = str2double(junk{3});
            elseif cur_aoa == '8'
                aoa8(k, n) = str2double(junk{3});
            else
                aoa12(k, n) = str2double(junk{3});
            end
            % iterate
            k = k + 1;
            line = fgetl(fid);
        end
        
        fclose(fid);
    end
end

%% Saving to .csv files

csvwrite('aoa0.csv',aoa0)
csvwrite('aoa4.csv',aoa4)
csvwrite('aoa8.csv',aoa8)
csvwrite('aoa12.csv',aoa12)