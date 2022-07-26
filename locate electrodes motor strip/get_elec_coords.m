%% Extract surface norms and E-field norms from CFM
% (c) Carys Evans, UCL
% carys.evans@ucl.ac.uk
% July 2022

%This script:
%Extracts electrode coordinates from ROAST.
%Requires modified roast scripts: 'roast2.m' and 'electrodePlacement2.m'
%Modified roast scripts must be in the roast folder
%roast2.m runs roast model until electrode placement.
%'electrodePlacement2.m outputs electrode coordinates.

%% Define file names and directories

%set directories
roastdir = 'C:\Matlab2018b\roast-3.0';    %roast directory
MRIdir = 'D:\MRI\DATA\';                                    %path to MRIs
roastMods = 'D:\ROAST\MODELS\';                             %path to roast models
datadir = 'D:\SAVE\DATA\HERE\';                             %path to save data

%get subject MRI filenames
cd(MRIdir)
subjs = dir('subject*.nii'); subjs = {subjs.name};
% subjs = extractBefore(subjs,'.')';
% subj_sort = extractAfter(subjs,'t');
% subj_sort = str2double(subj_sort);
% [subj_sort, idx] = sort(subj_sort);
% for i = 1:length(subjs)
% subjs_new(i,1) = subjs(idx(i));
% end
% subjs = subjs_new;
% clear idx subj_sort subjs_new

%roast model specifications. This will feed into roast2.m
montage = {'CP3FCZ','CPZFC3','C1FP2'};
montagespecs = [{'CP3',2,'FCz',-2};{'CPz',2,'FC3',-2};{'C1',2,'Fp2',-2}];

%output filename
savefilename = 'allPTs_electrodeCoords.mat';

%% Get electrode coordinates for each montage

for i = 1:length(subjs)
    cd(roastMods)
    %load MRI
    myvol =spm_vol(char(strcat(roastMods,subjs{i},'_ras.nii')));
    mri = spm_read_vols(myvol); 
    
    for j = 1:length(montage)
        cd(roastMods)
        
        if ~exist(char(strcat(subjs(i),'_',montage(j),'_2mA.mat')))
            continue
        else
            myroast = char(strcat(roastMods,subjs(i),'.nii'));
            cd(roastdir)
            %outputs Anode and Cathode coordinates in voxel space
            [myelectrodecoordA, myelectrodecoordC] = roast2((myroast),montagespecs(j,:), 'simulationTag', (strcat(montage{j},'_2mA')), 'elecsize', [17 2]);
         
            %Anode and Cathode coordinates in MNI space (xyz)
            myanodesxyz = myvol.mat(1:3,1:3) * myelectrodecoordA' + myvol.mat(1:3,4);
            mycathodesxyz = myvol.mat(1:3,1:3) * myelectrodecoordC' + myvol.mat(1:3,4);
            
            mysubjs(i) = subjs(i);
            
            electrodeCoord_data.(montage{j}).anodeCoords(i,:) = myelectrodecoordA;
            electrodeCoord_data.(montage{j}).anodeCoordsxyz(i,:) = myanodesxyz;
            electrodeCoord_data.(montage{j}).cathodeCoords(i,:) = myelectrodecoordC;
            electrodeCoord_data.(montage{j}).cathodeCoordsxyz(i,:) = mycathodesxyz;
            electrodeCoord_data.(montage{j}).subjs = mysubjs';
            
            clear myelectrodecoordA myelectrodecoordC myanodesxyz mycathodesxyz
        end
        
    end
    clear myvol mri
end

%save data
cd(datadir)
save(savefilename,'electrodeCoord_data');
