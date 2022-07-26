%% Define Motor Strip coordinates on GM Surface
% (c) Carys Evans, UCL 
% carys.evans@ucl.ac.uk
% July 2022

%This script:
%%Loads GM surface to allow manual selection of Motor Strip coordinates
%User defines motor strip using MNI coordinates (see line 75)
%Produces .mat file containing the medial and lateral coordinates of the
%Motor Strip (i.e. precentral gyrus)

%% Define file names and directories

%set directories
addpath(genpath('C:\Matlab2018b\spm12\'));          %spm
surfPATH = 'D:\PATH\TO\SURFACE\DATA\';              %surface path
cfmPATH = 'D:\PATH\TO\ROAST\MODELS\';               %ROAST models

%set files
cfmResultFile = 'CP3FCZ_2mA_roastResult.mat';       %ROAST result file (after subjectID)
gmfilename = 'ro_white-pial';                       %grey matter mask file
ROIfile = 'M1coord.mat';                            %ROI data                    

%desired filename
CoordFile = 'M1_motorstrip_coords.mat';

%% SUBJECT MRI.NII
%=================
%subject folder names
cd(surfPATH)
k = dir('1*'); subj={k.name}'; clear k

%subjects equivalent name in roast
%nb: only required if subject names for surface data & roast models do not match
cd(cfmPATH) 
k = dir('s*_*.nii'); k = unique(extractBetween({k.name},'subject','_'))';
k = str2double(k); k = sort(k);
for i = 1:length(k)
subj_cfm{i} = ['subject', num2str(k(i))];
end
clear i k


%% ==============================================
% CREATE GM SURFACE & EXTRACT SURFACE NORM AND E-FIELD DIRECTION
%%===============================================

for sub = 1:length(subj_cfm) %change to subj if subj_cfm not required
     
    %E-FIELD RESULT.MAT from Roast
    %=================
    cd(cfmPATH)
    cfmData = [sprintf('%s_', subj_cfm{sub}), cfmResultFile];
    if isfile(cfmData)
        load([cfmPATH,cfmData])
    else
        continue
    end
    
    %MRI .nii volume
    ef_vol = spm_vol([cfmPATH, sprintf('%s_ras.nii',subj_cfm{sub})]);
    
    
    %GREY MATTER SURFACE.GII
    %=================
    surface = [surfPATH, sprintf('%s/%s/surf/%s.gii', subj{sub},subj{sub}, pialfilename)];
    
    %generate gifti from surface file
    sl=gifti(surface);
    
    %select Motor Strip coordinates
    cd([surfPATH, sprintf('%s/%s', subj{sub},subj{sub})]);
        
    %load Coord file if created previously or plot surface & identify Motor
    %Strip coordinates. NB: define coordinates by selecting 'data tips' on
    %the figure window and clicking on desired location on GM surface.
    %Enter coordinates (rounded)
    if isfile(CoordFile)
         load(CoordFile); %check if coordinates already created previously
    else %pause;
        ROI_centre_coord_approx = importdata(ROIfile);
        ROI_centre_ind = knnsearch(sl.vertices,ROI_centre_coord_approx);
        ROI_vert_ind = [ROI_centre_ind];
        meanROI_vertices = sl.vertices(ROI_vert_ind,:);
      
        figure = plot(sl);
        hold on; scatter3(meanROI_vertices(:,1),meanROI_vertices(:,2),meanROI_vertices(:,3),'c*');

        %set(gcf,'Position', [1763 12 1293 840])
        M1medial = input('define most MEDIAL point of M1 as 3 numbers (e.g. -36 -28 56): ', 's');
        M1medial = str2num(M1medial);       
       
        M1lateral = input('define most LATERAL point of M1 as 3 numbers (e.g. -36 -28 56): ', 's');
        M1lateral = str2num(M1lateral);
        
        save(CoordFile, 'M1medial', 'M1lateral');

        close
     end
    
end