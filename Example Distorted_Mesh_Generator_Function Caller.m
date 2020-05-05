%A script written to call the function version of the distorted mesh
%generator, in order to easily generate swept data sets.


% 0 = No Deviation
% 1 = Centre Point Deviation (LEGACY)
% 2 = All Points Deviation
% 3 = No Deviation with tears
% 4 = Deviated with tears
% 5 = No Deviation, fold stiffness varying
% 6 = Deviation with varying fold stiffness
% 7 = Tears with varying fold stiffness
% 8 = Deviation, tears and varying fold stiffness
tear = [0:0.02:0.5];
%tear = 0;
%dev = [0:0.1:3];
dev = 0;
%stiff = [[0:0.1:4],[4.5:0.5:20]];
stiff = 0;
repeats = 8;

if(length(dev) > 1)
    method = 2;
    for kk = 1:length(dev)
        Distorted_Mesh_Generator_Function(method,dev(kk),tear,stiff,repeats);
    end
elseif(length(tear) > 1)
    method = 3;
    for jj = 1:length(tear)
        Distorted_Mesh_Generator_Function(method,dev,tear(jj),stiff,repeats);
    end
elseif(length(stiff) > 1)
    method = 5;
    for ii = 1:length(stiff)
        Distorted_Mesh_Generator_Function(method,dev,tear,stiff(ii),repeats);
    end
end
