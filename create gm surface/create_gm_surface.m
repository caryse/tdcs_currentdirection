%% Combine L & R Hemisphere surfaces & Correct RAS offset
% (c) Carys Evans, UCL 
% carys.evans@ucl.ac.uk
% July 2022

%This script:
%Corrects RAS offset introduced by freesurfer for L & R hemisphere surfaces
%Combines L & R Hemisphere surfaces to produce one surface
%RAS offset surfaces are prefixed with "ro_"
%Creates grey matter surface gifti, named "ro_white-pial.gii"

%Script also requires scripts 'combine_surfaces.m', 'write_surf_gifti.m'
%these scripts can also be found here https://github.com/jbonaiuto/MEGsurfer 


%% ============================================
%STEP 1 - IN TERMINAL
%==============================================

%%OBTAIN RAS OFFSET FROM FREESURFER
%nb: may need to setup freesurfer in terminal 1st

%AUTOMATIC METHOD:
%%paste following into bash terminal (e.g Ubuntu or Terminal)
%%script does following:
%%- retrieves subject file path
%%- converts surfaces to gifti
%%- obtains ras offset and saves to txt file 'rasoffset'
%%- outputs subject directory to check correct

%%Start Freesurfer
%export FREESURFER_HOME="/path/to/freesurfer"
%export SUBJECTS_DIR=/path/to/subjects_dir/
%source $FREESURFER_HOME/SetUpFreeSurfer.sh

%%Convert surfaces and obtain ras offset
% for dir in /mnt/d/NIHdata/*/*/ %(subject file path)
% do #for loop
% cd $dir/surf
% mris_convert lh.pial lh.pial.gii
% mris_convert rh.pial rh.pial.gii
% mris_convert lh.white lh.white.gii
% mris_convert rh.white rh.white.gii
% mri_info --cras $dir/mri/orig.mgz > rasoffset.txt
% echo $dir
% done 

%MANUAL METHOD:
   %type following to get 3 element vector:
%mri_info --cras /usr/local/freesurfer/subjects/subj_id/mri/orig.mgz'

%%CONVERT SURFACES TO GIFTI 
  %type the following for each surface being converted:
% mris_convert lh.pial lh.pial.gii
% mris_convert rh.pial rh.pial.gii
% mris_convert lh.white lh.white.gii
% mris_convert rh.white rh.white.gii

%NB: surface files must exist in original freesurfer directory in order to
%accurately convert images

 
%% ============================================
%STEP 2 - IN MATLAB
%==============================================

%set directories:
freepath = 'D:\FREESURFER\OUTPUT\';              %freesurfer output 
addpath(genpath('D:\Matlab2018b\spm12'));       %spm pathway
addpath(genpath('D:\PATHWAY\FOR\SCRIPTS\'));    %scripts pathway
cfmPATH = 'D:\ROAST\MODELS\';                   %ROAST models pathway

%freesurfer surface & hemisphere filenames
hemfile = {'lh','rh'};             %hemisphere .pial.gii files
surftype = {'pial','white'};       %surface being converted. Also filename. 

%subject folders
cd(freepath)
k = dir('1*'); subj={k.name}';      %subject folders begin with '1'


%% Read in each hemisphere's gifti file and adjust for RAS offset
for sub=1:length(subj)        %Loop for each subject
    
    surfdir = [freepath, sprintf('%s/%s/surf/',subj{sub},subj{sub})];
    
    %specify RAS offset
    ras_offset=dlmread([freepath, sprintf('%s/%s/surf/rasoffset.txt',subj{sub},subj{sub})]);
    
    for i=1:size(surftype,2)
        
        for j=1:size(hemfile,2)   %Repeat for each hemisphere
            g=gifti([freepath, sprintf('%s/%s/surf/%s.%s.gii',subj{sub},subj{sub},hemfile{j},surftype{i})]);
            % Set transformation matrix to identify
            g.mat=eye(4);
            g=set_mat(g,'NIFTI_XFORM_UNKNOWN','NIFTI_XFORM_TALAIRACH');
            % Apply RAS offset
            g.vertices=g.vertices+repmat(ras_offset,size(g.vertices,1),1);
            save(g,[freepath,sprintf('%s/%s/surf/ro_%s.%s.gii',subj{sub},subj{sub},hemfile{j},surftype{i})]);
            
            clear g
        end
        
    %% Combine L & R Hemisphere surfaces
    cd(surfdir)
   
    lh=sprintf('ro_lh.%s.gii',surftype{i});
    rh=sprintf('ro_rh.%s.gii',surftype{i});
    combined=sprintf('ro_%s.gii',surftype{i});
    combine_surfaces(lh, rh, combined);
    
    clear surfdir
    
    end

    
    %% Combine white & pial surfaces to create grey matter surface
    white = 'ro_white.gii';
    sl_w = gifti(white);
    pial = 'ro_pial.gii';
    sl_p = gifti(pial);
    sl_combine = sl_p;
    sl_combine.vertices = (sl_w.vertices+sl_p.vertices)/2;
    save(sl_combine, 'ro_white-pial.gii');
    
    
    
end

clearvars -except i sub

% %% combine white & pial surfaces
% 
% white=fullfile('ro_white.gii');
% pial=fullfile('ro_pial.gii');
% combined=fullfile('ro_white-pial.gii');
% combine_surfaces(white, pial, combined);
