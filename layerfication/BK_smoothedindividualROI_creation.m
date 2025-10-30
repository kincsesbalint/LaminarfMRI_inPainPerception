function BK_smoothedindividualROI_creation(subid,region)
% This function is called by its own. This is a modification of BK_ROI_creation_GROUPLVLactivationconj function.
% This function aim to convert to surface space the individual smoothed
% ROIs and then load it in the matlab space. It should return for each ROI
% an index number that reflect the columns which need to be selected from
% the already sampled matrices. Therefore, it can be later loaded by the
% fitting procedure to and restrict the averaging across those columns. 
% 
% TODO the rest from here 
% The following steps are implemented here:
%   1. load boundaries to matlab which coming from freesurfer to the matrix layer_boundaries [vertex,[x,y,z],[white,pial]] : 3D matrix of the surfaces(VPF_load_layer_boundaries)
%   (checked: does what I thought it does)
%   2. load ROIs coming from the group average activation (BK_convert_load_ROIs)
%   (checked: does what I thought it does)
%   3. transform boundaries to matlab space (from freesurfer space), the
%        resulting space will be in "voxel coordinate".
% 
% 
% .
% 
% 
% The ROI
%
%INPUTs:
% subid [str]                : subject id
% subpath [str]              : path to the derivatives folder
% fspath [str]               : path to the freesurfer folder
% T1path [str]               : full path and filename of the T1 image
% visualize
% region
% sampledimg
%OUTPUTs: a saved file in the indiviudal local folder (hardcoded, see
%below)
% EXAMPLE CALL:
% for subj={7349}; BK_smoothedindividualROI_creation(subj,'S1'); end
%
% Balint Kincses
% 12.2024

fspath = 'E:\pain_layers\main_project\derivatives\freesurfer\';

if subid{1} <= 7405
    subpath = 'E:/pain_layers/main_project/derivatives/pipeline/';
else
    subpath = 'D:/main_project/derivatives/pipeline/';
end
if subid{1}==7356
    T1path = [subpath num2str(subid{1}) '' ...
        '\ses-01\anat\presurf_MPRAGEise\presurf_UNI\UNI_MPRAGEised_biascorrected.nii'];                  
else
    T1path = [subpath num2str(subid{1}) '' ...
          '\ses-01\anat\presurf_MPRAGEise\presurf_UNI\UNI_MoCo_MPRAGEised_biascorrected.nii'];
end

if ~ischar(subid)
    subid = char(string(subid));
end
%hardcode the number of layers what we use for interpolation
N_layers = 20;
    
%local folder where individual subject folders are located
roipath='C:\Users\lenov\Documents\layerfMRI_DATA\groupavg_correctBET\';
% the output of this function can be saved.
outputpath=[roipath subid '\smootheddata\'];
%     sampledfilenm=[outputpath 'interimdata_rwls_' region sampledimgtype '.mat'];
    %check if the output is already exist. If so, and no visualization is
    %set, tehn it will exit:
%     if ~exist(sampledfilenm ,'file')
%         BK_displaytxt('The interim file does NOT exist, so we create one!')
%     end
        %load boundaries coming from freesurfer
layer_boundaries = VPF_load_layer_boundaries(subid,fspath);
        
%load ROIs 
try
    ROIs = BK_convert_load_ROIs(subid,roipath,fspath,size(layer_boundaries,1),region);
catch ME            
    BK_displaytxt('Set correctly the mask! maybe only unzipping is necessary?')
%             throw(ME)
    return
end
        fprintf('Size of different ROIs:%i\n', [sum(ROIs,1)])
        % Transform boundaries to matlab space (from freesurfer space), the resulting space will be in "voxel coordinate".
        %this is crucial step and need to be checked thouroughly. the
        %sampling depends on this:!!
        [layer_boundaries,T1_mat,old_layer_boudnaries] = VPF_transform_layers_to_matlab_space(layer_boundaries,T1path);

        %this would be the visualization of the 
%         if strcmp(visualizationtype,'no') || strcmp(visualizationtype,'ROI')
%             if strcmp(visualizationtype,'no')
%                 BK_displaytxt('No visualization of the ROI was selected but we DO the sampling.')
%             elseif strcmp(visualizationtype,'ROI')
                 BK_displaytxt('Visualization of the ROI + sampling')
                 BK_plotROIverticesinmatlabspace(T1path,ROIs,layer_boundaries,region,outputpath)
%             end
            fprintf(sprintf('Starting layer sampling...\n'));
%             %sample layers
%             tic
% %             fprintf('time for layer sampling:')
% %             if strcmp(sampledimgtype,'derived')
% %                 disp("We sample the derived images!?")
% %                 %todo: this was not refactored as I do not use eventually
% %                 %this approach:
% % 
% %                 [columnspecificts,layeractivation,allroicolumnsize,sampled_img_list]= BK_select_active_columns_basedontmap(subid,subpath,layer_boundaries,ROIs,N_layers, T1_mat);
% %                 interimdata_columns ={columnspecificts,layeractivation,allroicolumnsize,sampled_img_list};
% %             elseif strcmp(sampledimgtype,'raw')
% %                 disp("PLUGIN the external HD!! we go with raw data sampling!!")
% %                 [columnspecificts,layeractivation,allroicolumnsize]= BK_select_active_columns_basedonts(subid,subpath,layer_boundaries,ROIs,N_layers, T1_mat,old_layer_boudnaries);
% %                 interimdata_columns ={columnspecificts,layeractivation,allroicolumnsize};
% %             end
% %             
% %             save(sampledfilenm,'interimdata_columns','-v7.3');
% %             fprintf('Total time with sampling and saving:')
% %             toc
% %         elseif strcmp(visualizationtype,'onlyROI')
% %             BK_displaytxt('Only visualization of the ROI was selected !!NO sampling!!')
% %             BK_plotROIverticesinmatlabspace(T1path,ROIs,layer_boundaries,region,outputpath)
% %         elseif strcmp(visualizationtype,'sampledROI')
% %             BK_displaytxt('This is invalid with the argument sample! Try with the glmestimate argument.')
%         end
end

%% load layers function (checked -keep it)
function layer_boundaries = VPF_load_layer_boundaries(subid,fspath)
    %loads the layer boundaries (white/pial) generated by Freesurfer
    %INPUT:
    %subid [str]                : subject id
    %fspath [str]               : path to the freesurfer folder
    %OUTPUT:
    % layer_boundaries [vertex,[x,y,z],[white,pial]] : 3D matrix of the surfaces
    %this part loads in the vertex_coordinates of the surface file
    %(https://github.com/fieldtrip/fieldtrip/blob/master/external/freesurfer/read_surf.m)
    %and create a huge matrix which each row represent the cooridate of a
    %vertex (starting from left hemisphere and continou with right hemishpere)
    %and teh 3rd dimension contains white matter and pial information.
    fn = [fspath '\' subid '\surf\lh.white'];
    if strcmp(subid,'7356')
        asciiload=1; %it does not really change as the header is not readable
    else
        asciiload=0;
    end
    layer_boundaries = read_surf(fn,asciiload);
    fn = [fspath '\' subid '\surf\rh.white']; 
    tmp = read_surf(fn,asciiload);
 
    layer_boundaries  = cat(1,layer_boundaries ,tmp);
    
    
    fn = [fspath '\' subid '\surf\lh.pial'];
    tmp = read_surf(fn,asciiload);

    fn = [fspath '\' subid '\surf\rh.pial'];
    tmp = cat(1,tmp,read_surf(fn,asciiload));
    layer_boundaries = cat(3,layer_boundaries,tmp);

end

%% convert the ROI to surface (wb_command) - main function (checked -keep it)
function ROIs = BK_convert_load_ROIs(subid,roipath,fspath,N_vertex,region)
    % This is an adaptation of VPF_convert_load_ROIs function. It converts the
    % ROIs to surfaces and loads them. As luckily it was already run by Viktor,
    % the freesurfer part can be commented out (but keep them to see the whole
    % picture), bc the output files are available. The checkfile function test
    % if the file is available, and output a message in case of its missing.
    % The 
    %INPUT:
    %   subid [str]                : subject id
    %   subpath [str]              : path to the derivatives folder
    %   fspath [str]               : path to the freesurfer folder
    %   N_vertex                   : the number of vertices comming from freesurfer
    %   region [str]               : the ROI derived from smoothed group average and transformed back to the indiviudal space after intersecting with some anatomical mask. Usually left and right are differentiated.t
    %
    %OUTPUT:
    %ROIs [N_vertex,N_ROI] : The ROI mask in loaded surface space. 
    %                        N_roi is the number of regions provided
    %                        as an arfument(left and right are
    %                        concatenated, but one can provide separeate
    %                        ROI masks, so only one hemispheres are
    %                        processed)
    
    % fs_base = ['export FREESURFER_HOME=/usr/local/freesurfer/6.0.0; ' ...
    %     'export SUBJECTS_DIR=' fspath '; '...
    %     'source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
    

    
    localizer = dir(fullfile([roipath subid '\smootheddata\'], '**',['roi_' region '*thr0.nii']));
    %I assume here that the first localizer is always the left
    %one(alphabetic order), double check and trow an error if not
    idxtostart=strfind(localizer(1).name,region)+length(region);
    if ~strcmp(localizer(1).name(idxtostart:idxtostart+3),'left')
        BK_displaytxt('The order of ROI was not left and right! Double check what is the issue')
    end
    if isempty(localizer)        
        return
    end
    N_total_ROIs = length(localizer);
    %first convert freesurfer's xh.white to xh.white.gii
    %theoretically this step is needed so wb_command then capable of using the
    %surface data (othe solruion would be to use mri_vol2surf??) but most
    %probably that result in different result--> as the file exist, it is good
    %as it is.
    %(https://www.mail-archive.com/freesurfer@nmr.mgh.harvard.edu/msg66816.html)
    
    %check if these files exist
    % cmd = ['mris_convert ' '"' fspath '/' subid '/surf/lh.white" lh.white.gii'];
    % system([fs_base cmd]);
    BK_checkfile([fspath '\' subid '\surf\lh.white.gii'])
    % cmd = ['mris_convert ' '"' fspath '/' subid '/surf/rh.white" rh.white.gii'];
    % system([fs_base cmd]);
    BK_checkfile([fspath '\' subid '\surf\rh.white.gii'])
    
    ROIs = zeros([N_vertex,N_total_ROIs],'logical');
    for ROI = 1:N_total_ROIs
        data = [localizer(ROI).folder '\' localizer(ROI).name];
        wls_data=strrep(data, '\', '/');
        if data(1)=='E'
            wls_data=strrep(wls_data,'E:','/mnt/e');
        elseif data(1)=='D'
            wls_data=strrep(wls_data,'D:','/mnt/d');
        elseif data(1)=='C'
            wls_data=strrep(wls_data,'C:','/mnt/c');
        end
        wls_fspath=strrep(fspath, '\', '/');
        wls_fspath=strrep(wls_fspath,'E:','/mnt/e');
        
        try %try to load to ROI surfaces
            ROIs(:,ROI) = VPF_average_mask_and_thres(data(1:end-4));
        catch %apparently, not created yet so create them first
            hemis = {'lh','rh'};
            for hemi = 1:2
                cmd = ['wsl wb_command -volume-to-surface-mapping ' '"' wls_data '"'...
                    ' "' wls_fspath '/' subid '/surf/' hemis{hemi} '.white.gii' '"'...
                    ' "' wls_data(1:end-4) '_' hemis{hemi} '.shape.gii' '"'...
                    ' -cubic'];
                system(cmd);
            end
            ROIs(:,ROI) = VPF_average_mask_and_thres(data(1:end-4));
        end
    end
end

% convert the ROI to surface (wb_command) - supporting function
function cdata = VPF_average_mask_and_thres(in)
    %loads the ROI surfaces. As they are transformed from voxel to surface space, 
    %some interpolation errors occur which are fixed here. Finally, the left 
    %and right hemispheres are concatenated
    %INPUT:
    %in [str]                : full path and name of ROI WITHOUT extension (i.e.
    %                          without .nii)      
    %
    %
    %OUTPUT:
    %cdata [N_vertex]        : logical array of the ROI in surface space
    
    ROI = gifti([in '_lh.shape.gii']);
    cdata_lh = ROI.cdata;
    cdata_lh(abs(cdata_lh)>=1e-2) = 1;
    cdata_lh(abs(cdata_lh)<1e-2) = 0;
    
    ROI = gifti([in '_rh.shape.gii']);
    cdata_rh = ROI.cdata;
    cdata_rh(abs(cdata_rh)>=1e-2) = 1;
    cdata_rh(abs(cdata_rh)<1e-2) = 0;
    
    cdata = logical(cat(1,cdata_lh,cdata_rh));
end

%% transform between fs and matlab (voxel) spaces. - checked the output of the visualization and it is good,keep it
%this is important as the height of the columns can be calculate from here
function [layer_boundaries, T1_mat,old_layer_boudnaries ] = VPF_transform_layers_to_matlab_space(layer_boundaries,T1path)
    %Transforms the layers from freesurfer space to matlab space. After that, 
    %signal can be sampled using spm_sample_vol
    %INPUT:
    %layer_boundaries [vertex,[x,y,z],[white,pial]] : 3D matrix of the surfaces
    %T1path [str]                                   : full path and filename of 
    %                                                 the T1 image
    %
    %OUTPUT:
    %layer_boundaries        : 3D matrix of surfaces in matlab space
    %
    
    sz = size(layer_boundaries);
    
    %bring surfaces into matlab space
    T1_mat = spm_get_space(T1path); %Get/set the voxel-to-world mapping of an image
    old_layer_boudnaries=layer_boundaries ; %just to compare them for transformation

    boundaries_2_cat = ones([sz(1) 1 sz(end)]);
    layer_boundaries = cat(2,layer_boundaries,boundaries_2_cat);
    

    for ii = 1:sz(end)
        tmp = T1_mat\squeeze(layer_boundaries(:,:,ii)).'; 
        %  \ - left matrix division (also known as the backslash operator):
        %  solves the following equstion: T1_mat X tmp =
        %  squeeze(layer_boundaries(:,:,ii)).'--> T1mat^(−1)×laybound
        % .' - non-conjugate transpose (simple transpose)
        layer_boundaries(:,:,ii) = tmp.';
    end
    
    layer_boundaries= layer_boundaries(:,1:3,:);
        %for subject 7356,rotation issue - solved with an additional
        %translation when creating the .gii file. see detailed description
        %of the issue and its solution in my documentation.
% %        img_T1 = spm_read_vols(spm_vol(T1path));
% % %        img_func=spm_read_vols(spm_vol('E:\pain_layers\main_project\derivatives\pipeline\7356\ses-02\func\layers\run1\func\mag_POCS_r1_1000_Warped-to-Anat.nii.gz'));
%     for slice = 180 %205
%         %this is the y direction/AP
%         idx = find(layer_boundaries(:,1,1)>slice & layer_boundaries(:,1,1)<slice+2);
%         idx3 = find(layer_boundaries(:,1,2)>slice & layer_boundaries(:,1,2)<slice+2);
%     
%         figure;imagesc(squeeze(img_T1(slice,:,:)), [0 1200]), colormap(gray(256)), title(['Slice ' num2str(slice)])
% %         figure;imagesc(squeeze(img_func(:,:,slice)), [0 1200]), colormap(gray(256)), title(['Slice ' num2str(slice)])
%         hold on;
%     
%         plot(layer_boundaries(idx,3,1),layer_boundaries(idx,2,1),'r.');
%         plot(layer_boundaries(idx3,3,2),layer_boundaries(idx3,2,2),'y.');
%     
%         legend('white','pial');
%         hold off
%         %this is the z direction
%         slice=250;
%         idx = find(layer_boundaries(:,2,1)>slice & layer_boundaries(:,2,1)<slice+2);
%         idx3 = find(layer_boundaries(:,2,2)>slice & layer_boundaries(:,2,2)<slice+2);
%     
%         figure;imagesc(squeeze(img_T1(:,slice,:)), [0 1200]), colormap(gray(256)), title(['Slice ' num2str(slice)])
% %         figure;imagesc(squeeze(img_func(:,:,slice)), [0 1200]), colormap(gray(256)), title(['Slice ' num2str(slice)])
%         hold on;
%     
%         plot(layer_boundaries(idx,3,1),layer_boundaries(idx,1,1),'r.');
%         plot(layer_boundaries(idx3,3,2),layer_boundaries(idx3,1,2),'y.');
%         hold off
% 
%         %this is the x direction
%     end
end

%% sample columnswise information based on ROI in the raw fmri data (output is ts) - it is checked,keep it
function [columns_out,layers_out,allroicolumnsize] = BK_select_active_columns_basedonts(subid,subpath,layer_boundaries,ROIs,N_layers, T1_mat,old_layer_boudnaries)
% This is a modificaiton of VPF core sampling function (VPF_sample_layers).
% It does three things:
%   1. get the timeseries signal from the deep 75% of the column (unbias (less effect of venous bias) estimation of
%   columnar activation) 
%   2. get the timeseries signal from the whole column (layerwise)
%   3. calculate the "column" thickness from all those vertices where we
%   sampled.
%
% Is uses the previously defined ROIs.
% Searches for folders with 'run' in their name and
% samples the containing *Warped-to-Anat.nii.gz images within the predefined 
% ROIs.
%
% INPUT:
%subid [str]                : subject id
%subpath [str]              : path to the derivatives folder
%layer_boundaries           : 3D matrix of surfaces in matlab space
%ROIs [N_vertex,N_ROI] :    : 2D logical matrix containing the ROIs in
%                             surface space
%N_layers [int]             : number of layers (default: 20)
%T1_mat : affine matrix of T1img
%OUTPUT:
% columns_out {N_ROIx1} [N_columns  N_totalvol]
% layer_out {N_ROI x1} [N_columns N_layers N_totalvol]
% allroicolumnsize {N_ROI} the distribution of the columnar thickness

    currpath = pwd;
    SAMPLING_ORDER = 3;
    runs = dir([subpath '\' subid '\ses-02\func\layers\run*']);
    N_run = numel(runs);
    N_ROI = size(ROIs,2);
    allroicolumnsize = cell(N_ROI);
    for run = 1:N_run
            
        cd([runs(run).folder '\' runs(run).name '\func']);
        sampled_img_list{run} = dir([runs(run).folder '\' runs(run).name '\func\*Warped-to-Anat.nii.gz']);
        N_vols(run) = length(sampled_img_list{run});
    end
    tic
    layers_out=cell(N_ROI,1);
    columns_out=cell(N_ROI,1);
    fullvol=1;
    for run=1:N_run
        fprintf(sprintf(['run ' num2str(run) '...']));
        for vol = 1:N_vols(run)
            cd([runs(run).folder '\' runs(run).name '\func']);
    %         sampled_img = load_nifti([sampled_img_list(vol).name]).vol;
            sampled_img = niftiread([sampled_img_list{run}(vol).name]); %I ahve chekcd on a random subject random image and it is the same as the previous line(which was originally implemented by VPF)
            for ROI = 1:N_ROI
                ind = find(ROIs(:,ROI)); %columns under the group lvl mask
                %% randomly pick 200 voxelpairs
                %TODO remove this, this is only for the auditory cortex, to
                %have a "null" ROI.
%                 randomIndices = randperm(length(ind), 200);
%                 ind = ind(randomIndices);
                %%
                tmp = zeros(size(ind,1),N_layers);
                boundaries = layer_boundaries;
                if run==1 && vol==1
                    columns_out{ROI,1} = zeros([length(ind),sum(N_vols)]);
                    layers_out{ROI,1} = zeros([length(ind),N_layers,sum(N_vols)]);
                    columnsize=nan(1,size(ind,1));
                end
               
                %I am not sure what is the vertex coordinate here, but I can
                %         %believe that is in mm (mean+/-sd=3.55+/-1.1)
                %         %double check if that is the case. 
                %         % think about to use only columns which are within a range, which
                %         % suggest that the definition of that column makes sense(too big or
                %         % too small is OK??)
                for k = 1:size(ind,1) 
                    X  = linspace(boundaries(ind(k),1,1),boundaries(ind(k),1,2),N_layers);
                    Y  = linspace(boundaries(ind(k),2,1),boundaries(ind(k),2,2),N_layers);
                    Z  = linspace(boundaries(ind(k),3,1),boundaries(ind(k),3,2),N_layers);
                    if vol==1 && run==1 %calculate the column height. This need to be calculated only once.
                        columnsize(k)=sqrt((X(end)-X(1))^2+(Y(end)-Y(1))^2+(Z(end)-Z(1))^2) * T1_mat(1,1);
%                         columnsize_orig(k)=sqrt((old_layer_boudnaries(ind(k),1,2)-old_layer_boudnaries(ind(k),1,1))^2+...
%                         (old_layer_boudnaries(ind(k),2,2)-old_layer_boudnaries(ind(k),2,1))^2+...
%                         (old_layer_boudnaries(ind(k),3,2)-old_layer_boudnaries(ind(k),3,1))^2);

                    end
                    tmp(k,:) = spm_sample_vol(sampled_img,X,Y,Z,SAMPLING_ORDER); %this seems to work fine if we use the already loaded image, but one can give the 
                end
                if vol==1 && run==1
                    allroicolumnsize{ROI}=columnsize;
                end
                layers_out{ROI}(:,:,fullvol) = tmp;
                columns_out{ROI}(:,fullvol) = mean(tmp(:,1:15),2); % ,'omitnan' - 15 is hardcoded here as Peter recommended to skip the top 0.75mm, which would be in case of 3mm avg thickness the same as 15/20.It can be adjusted later by the thickness of the column as we can easily calculate that.
                
            end
            fullvol=fullvol+1;
        end
        toc
        fprintf(sprintf(' done \n'));
    end
cd(currpath)
end

%% sample columnswise information based on ROI in the ß image fmri data - have not checked, but keep it
function [columns_out,layers_out,allroicolumnsize,sampled_img_list] = BK_select_active_columns_basedontmap(subid,subpath,layer_boundaries,ROIs,N_layers, T1_mat)
% This is a modificaiton of VPF core sampling function (VPF_sample_layers).
% It does two things:
%   1. get the t values from the deep 75% of the column (unbias (less effect of venous bias) estimation of
%   columnar activation) 
%   2. get the t values from the whole column (layerwise)
% Is uses the previously defined ROIs.
% Searches for local folders in which the corresponding tmaps are included.
%
% INPUT:
%subid [str]                : subject id
%subpath [str]              : path to the derivatives folder
%layer_boundaries           : 3D matrix of surfaces in matlab space
%ROIs [N_vertex,N_ROI] :    : 2D logical matrix containing the ROIs in
%                             surface space
%N_layers [int]             : number of layers (default: 20)
%T1_mat : affine matrix of T1img
%OUTPUT:
% columns_out [N_columns N_ROI N_totalvol]
% layer_out [N_columns N_layers N_ROI N_totalvol]
% allroicolumnsize {N_ROI} the distribution of the 

    currpath = pwd;
    SAMPLING_ORDER = 3;
%     runs = dir([subpath '\' subid '\ses-02\func\layers\run*']);
%     N_run = numel(runs);
    N_ROI = size(ROIs,2);
    allroicolumnsize = cell(N_ROI);

    layers_out=cell(N_ROI,1);
    columns_out=cell(N_ROI,1);
%      runs = dir([subpath '\' subid '\ses-02\func\layers\rwls_stats_compcor_UNsmoothed_hpf180']);    
    sampled_img_list_con4 = dir([subpath '\' subid '\ses-02\func\layers\rwls_stats_compcor_UNsmoothed_hpf180\con_0004*']); %interaction is not necessary, 4 would be the main effect of cognition
    sampled_img_list_con8 = dir([subpath '\' subid '\ses-02\func\layers\rwls_stats_compcor_UNsmoothed_hpf180\con_0008*']);% and 8 would be the main effect of pain.
    sampled_img_list_tstat4 = dir([subpath '\' subid '\ses-02\func\layers\rwls_stats_compcor_UNsmoothed_hpf180\spmT_0004*']); %interaction is not necessary, 4 would be the main effect of cognition and 8 would be the main effect of pain.
    sampled_img_list_tstat8 = dir([subpath '\' subid '\ses-02\func\layers\rwls_stats_compcor_UNsmoothed_hpf180\spmT_0008*']); %interaction is not necessary, 4 would be the main effect of cognition and 8 would be the main effect of pain.
    sampled_img_list=[sampled_img_list_con4;...
                      sampled_img_list_con8;...
                      sampled_img_list_tstat4;...
                      sampled_img_list_tstat8];
    
    nproc=length(sampled_img_list);
    for process=1:nproc % we have a bottom-up(8) and a top-down(4) process. Take care that the topdown should be reversed
        sampled_img = niftiread([sampled_img_list(process).folder '\' sampled_img_list(process).name]);
        for ROI = 1:N_ROI
                ind = find(ROIs(:,ROI)); %columns under the group lvl mask
                tmp = zeros(size(ind,1),N_layers);
                boundaries = layer_boundaries;

                if process==1 
                    columns_out{ROI,1} = zeros([length(ind),nproc]);
                    layers_out{ROI,1} = zeros([length(ind),N_layers,nproc]);
                    columnsize=nan(1,size(ind,1));
                end
             for k = 1:size(ind,1) 
                    %sanity check of the cooridnates 
                    X  = linspace(boundaries(ind(k),1,1),boundaries(ind(k),1,2),N_layers);
                    Y  = linspace(boundaries(ind(k),2,1),boundaries(ind(k),2,2),N_layers);
                    Z  = linspace(boundaries(ind(k),3,1),boundaries(ind(k),3,2),N_layers);
%                     if vol==1 && run==1 %calculate the column height. This need to be calculated only once.
%                         columnsize(k)=sqrt((X(end)-X(1))^2+(Y(end)-Y(1))^2+(Z(end)-Z(1))^2) * T1_mat(1,1); 
%                         %it seems that as we are in matlab space, we need
%                         %to transform the distance between the two
%                         %points (column bottom and top) back to world
%                         %space. This can be done simply applyint the affine
%                         %matrix but as it is the same for all subjects and
%                         %isovoxel resolution with a 0.75mm, we can simple
%                         %multiply the distance with this constant. I
%                         %hardcoded it!!!
%                     end
                    tmp(k,:) = spm_sample_vol(sampled_img,X,Y,Z,SAMPLING_ORDER);
              end
                if process==1
                    allroicolumnsize{ROI}=columnsize;
                end
    
                layers_out{ROI}(:,:,process) = tmp;
                columns_out{ROI}(:,process) = mean(tmp(:,1:15),2); % ,'omitnan' - 15 is hardcoded here as Peter recommended to skip the top 0.75mm, which would be in case of 3mm avg thickness the same as 15/20.It can be adjusted later by the thickness of the column as we can easily calculate that.
                
        end

    end
    
cd(currpath)
end

%% visualize ROI vertices in matlab space - checked ,keep it
function BK_plotROIverticesinmatlabspace(imgpath,ROIs,layer_boundaries,region,outputpath)
    img_T1 = spm_read_vols(spm_vol(imgpath));
    fig = figure;%('Visible', 'off');
        leftsidemask=find(ROIs(:,1));
        rightsidemask=find(ROIs(:,2));
        
        mostcolumns=mode(round(layer_boundaries([leftsidemask],3,1)));
        slice_idx=0;
%         zoom=30;
        for slice = mostcolumns-2:mostcolumns+3
            idx = find(layer_boundaries(:,3,1)>slice & layer_boundaries(:,3,1)<slice+1);
            idx3 = find(layer_boundaries(:,3,2)>slice & layer_boundaries(:,3,2)<slice+1);
            
%             idx_roimask = intersect(idx,leftsidemask);
            idx_roimask = intersect(idx,[leftsidemask; rightsidemask]);
%             idx3_roimask = intersect(idx3,leftsidemask);
            idx3_roimask = intersect(idx3,[leftsidemask; rightsidemask]);

            slice_idx = slice_idx + 1;
    
            % Create a subplot for each slice
            ax = subplot(2, 3, slice_idx); % 2 rows, ceil() for columns
            pos = get(ax, 'Position'); % Get current position
            pos(3) = pos(3) * 1.2;     % Increase width by 20%
            pos(4) = pos(4) * 1.2;     % Increase height by 20%
            set(ax, 'Position', pos);
            
            imagesc(squeeze(img_T1(:,:,slice)), [0 1200]), colormap(gray(256)), title(['Slice ' num2str(slice)])
%             zoomedincoord=round([min(layer_boundaries(idx_roimask,1,1))-zoom  ...
%                 max(layer_boundaries(idx_roimask,1,1))+zoom  ...
%                 min(layer_boundaries(idx_roimask,2,1))-zoom ...
%                 max(layer_boundaries(idx_roimask,2,1))+zoom]);
%             axis(zoomedincoord);
            
            
            hold on;
            
            plot(layer_boundaries(idx,2,1),layer_boundaries(idx,1,1),'r.');
            plot(layer_boundaries(idx3,2,2),layer_boundaries(idx3,1,2),'y.');
            
            


            plot(layer_boundaries(idx_roimask,2,1),layer_boundaries(idx_roimask,1,1),'g.');
            plot(layer_boundaries(idx3_roimask,2,2),layer_boundaries(idx3_roimask,1,2),'b.');
%               
            columntogether=intersect(idx_roimask,idx3_roimask);
            x1=layer_boundaries(columntogether,2,1);
            x2=layer_boundaries(columntogether,2,2);
            y1=layer_boundaries(columntogether,1,1);
            y2=layer_boundaries(columntogether,1,2);
            for togethervertices=1:length(columntogether)
                plot([x1(togethervertices), x2(togethervertices)], [y1(togethervertices), y2(togethervertices)], '-', 'LineWidth', 2, 'MarkerSize', 8,'Color','m');
            end

            legend('white','pial');
            axis off;
        end
    
        txtROIcolumns=['The number of vertices in the ROI:', num2str(sum(ROIs))];
        annotation('textbox', [0.5, 0.01, 0.1, 0.05], ...
                   'String', txtROIcolumns, ...
                   'HorizontalAlignment', 'center', ...
                   'VerticalAlignment', 'bottom', ...
                   'FontSize', 14, 'LineStyle', 'none');
        set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); 
        set(fig, 'PaperUnits', 'inches', 'PaperSize', [24, 32],'PaperPosition', [0 0 24 32]); % Match the figure size
        print(fig, [outputpath region '_ROI_forallposvertex.pdf'], '-dpdf', '-bestfit');
        savefig(fig, [outputpath region '_ROI_forallposvertex.fig']);  % Save as .fig
        close(fig);
end