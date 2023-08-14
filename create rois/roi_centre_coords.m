%% Define ROI centre coordinates on GM surface
% (c) Carys Evans, UCL 
% carys.evans@ucl.ac.uk
% July 2022

%This script:
%Loads GM surface to allow manual selection of ROI centre
%User defines ROI centre using MNI coordinates (see line 76)
%Produces .mat file containing coordinates for centre of each ROI

%% Define file names and directories

%set directories

%addpath(genpath('C:\Matlab2018b\spm12\'));         %spm
%surfPATH = 'D:\PATH\TO\FREESURFER\SUBJECT\DATA\';  %surface path
%cfmPATH = 'D:\PATH\TO\ROAST\MODELS\';              %ROAST models

%spm12 PATH
addpath(genpath('/Users/carysevans/Documents/MATLAB/spm12/'));          %spm

%surface data PATH
disp('select subject folder containing freesurfer data (i.e. folder containing "surf" folder)')
surfPATH = uigetdir; surfPATH = [surfPATH,'/'];

%ROAST current flow model PATH
disp('select folder containing current flow model data for individual subject')
cfmPATH = uigetdir; cfmPATH = [cfmPATH,'/'];

%set files
cfmResultFile = 'CP3FCZ_2mA_roastResult.mat';       %ROAST result file (after subjectID)
gmfilename = 'ro_white-pial';                       %grey matter mask file

%desired ROI filename. Repeat for each ROI
ROIfile = 'M1coord.mat';                            %e.g. M1 bank coordinates                     

showROI = 1; %Display selected ROI centrepoint on surface? 1 = yes, 0 = no.

%%
%SUBJECT MRI.NII
%=================

%subject folder names, freesurfer
cd(surfPATH)
%k = dir('1*'); subj={k.name}'; clear k
slashIdx = strfind(surfPATH, '/');
subj = extractBetween(surfPATH, (slashIdx(1,end-1)+1),(slashIdx(1,end)-1));
%subj = cell2mat(extractBetween(surfPATH, (slashIdx(1,end-1)+1),(slashIdx(1,end)-1)));
clear slashIdx


%subjects equivalent name in roast
%nb: only required if subject file names for surface data & roast models do not match
cd(cfmPATH) 
% k = dir('s*_*.nii'); k = unique(extractBetween({k.name},'subject','_'))';
% k = str2double(k); k = sort(k);
% for i = 1:length(k)
% subj_cfm{i} = ['subject', num2str(k(i))];
% end
% clear i k

k = dir(['*',cfmResultFile]);
subj_cfm = {extractBefore(k.name,['_',cfmResultFile])};


%% ==============================================
% CREATE GM SURFACE & EXTRACT SURFACE NORM AND E-FIELD DIRECTION
%%===============================================

for sub = 1:length(subj_cfm) %change to subj if subj_cfm not required
    
    %E-FIELD RESULT.MAT from Roast
    %=================
    cd(cfmPATH)
    cfmData = [sprintf('%s_', subj_cfm{sub}), cfmResultFile];
    
    %load roast model or skip subject if it does not exist
    if isfile(cfmData)
        load([cfmPATH,cfmData])
    else
        continue
    end
    
    %MRI .nii volume (ras transformed by roast)
    if  isfile([cfmPATH, sprintf('%s_ras.nii',subj_cfm{sub})])
        ef_vol = spm_vol([cfmPATH, sprintf('%s_ras.nii',subj_cfm{sub})]);
    else
        ef_vol = spm_vol([cfmPATH, sprintf('%s.nii',subj_cfm{sub})]);
    end

    
    %GREY MATTER SURFACE.GII
    %=================
    %surface = [surfPATH, sprintf('%s/surf/%s.gii', subj{sub}, gmfilename)];
    surface = [surfPATH, sprintf('/surf/%s.gii', gmfilename)];
    
    %generate gifti image from surface file
    sl=gifti(surface);
    
    %select ROI coordinates
    %cd([surfPATH, sprintf('%s',subj{sub})]); %subject folder within surface path  
    cd([surfPATH,'surf/'); %subject folder within surface path   

    
    %load ROI file if created previously or plot surface & identify ROI coordinates
    %nb: define ROI centre by selecting 'data tips' on the figure window and
    %clicking on desired location on GM surface. Enter coordinates (rounded)
    if isfile(ROIfile) 
        ROI_centre_coord_approx = importdata(ROIfile); 
    else 
        figure = plot(sl);
        fprintf('FIND ROI.\nSelect "Data Tips" on the figure navigation pane.\nClick on the cortical area where you want your ROI coordinates to be.\nThe XYZ values are your coordinates.\n');
        ROI_centre_coord_approx = input('Type your ROI centre here. Use round numbers (e.g. -36 -28 56): ', 's');
        ROI_centre_coord_approx = str2num(ROI_centre_coord_approx);
        save(ROIfile, 'ROI_centre_coord_approx');
        close
    end    

    % find index of closest vertex for ROI centre coords on gm surface
    ROI_centre_ind = knnsearch(sl.vertices,ROI_centre_coord_approx);
    % initialise ROI_vert_ind
    ROI_vert_ind = [ROI_centre_ind];
    ROI_vertices = sl.vertices(ROI_vert_ind,:);
    meanROI_vertices = sl.vertices(ROI_centre_ind,:);
    
    %plot ROI centre to confirm location
    if showROI == 1
        figure = plot(sl);
        hold on; scatter3(meanROI_vertices(:,1),meanROI_vertices(:,2),meanROI_vertices(:,3),'c*');
        
        answer = input('Keep ROI? (Y/N): ','s');
        if answer == 'N'
            fixROI(sub) = subj_cfm(sub);
        end
    end
    
    close 
      
end

%display whether any ROIs need redoing
if exist('fixROI','var')
    disp(['redo the following ROIs: ', fixROI])
else
    disp('all ROIs complete')
    
end
