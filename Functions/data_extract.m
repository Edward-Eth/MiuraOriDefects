function [timeenergy] = data_extract(rundocfolder,filename)
%Function to scan through .dat files line by line to find desired data and
%return it.

% rundocfolder = 'run_method0';
% filename = 'Miura_sheet_n1_m1_a30_b30.FOLD.dat';

%filename = fullfile(rundocfolder,filename);
FID = fopen(filename);
tcount = 1;
ecount = 1;

while(feof(FID)~=1)
    data = fgetl(FID);
    if(contains(data,'TOTAL TIME COMPLETED'))
        splitdata = split(data,' ');
        timesvals(tcount) = splitdata(29); %#ok<*AGROW>
        tcount = tcount + 1;
    end
    
    if(contains(data,'TOTAL STRAIN ENERGY (STRESS POWER)'))
        splitdata = split(data,' ');
        energyvals(ecount) = splitdata(end);
        ecount = ecount + 1;
    end
    
end

if(tcount>1)
    times = zeros(size(timesvals,1),size(timesvals,2));
    times = str2double(timesvals);
    energy = zeros(size(energyvals,1),size(energyvals,2));
    energy = str2double(energyvals);
    timeenergy = [times;energy];
else
    timeenergy = 0;
end

fclose(FID);