function [columnwisestat,columndistribution,layerwisestat] = BK_firstlvlanalysis(subid,subpath,fspath,T1path,visualize,region,N_layers)
%
% ToDo: THe function ToDo description
% 
% This is a adaptation of VPF's similarly called function (VPF__layer_sampling_pain_study)., but the functional ROIs(based on smoothed group level result)
% were used on the unsmoothed data, to find columns which activate.
% The ROI
% Function to sample layers within predefined regions (see VPF_ROI_creation.m). 
% It loads freesurfer surfaces, transforms ROI volumetric .nii to surfaces 
% and samples layers of images saved in folders with 'run' in their name
%INPUT:
%subid [str]                : subject id
%subpath [str]              : path to the derivatives folder
%fspath [str]               : path to the freesurfer folder
%T1path [str]               : full path and filename of the T1 image
%N_layers [int] optional    :number of layers to be sampled (default: 20)
%OUTPUT:
% activatedcolumns [column_index N_ROI contrastofinterest(cogn,pain)
% N_estimation(reml,rwsl)] : the fitted 
%
% 
% 
% Balint Kincses
% 12.2024
    if ~ischar(subid)
        subid = char(string(subid));
    end
    if nargin < 7
        N_layers = 20;
    end
    roipath='C:\Users\lenov\Documents\layerfMRI_DATA\groupavg_correctBET\';
    outputpath=[roipath subid '\functionalmasks\'];

% %     if ~exist([outputpath 'interimdata_rwls.mat'],'file')
% %         tic
%         %load boundaries coming from freesurfer
%         layer_boundaries = VPF_load_layer_boundaries(subid,fspath);
% %         fprintf('time for loading:')
%         %~0sec
% %         toc
%         % load ROIs coming from the group average activation
%         % time for converting:Elapsed time is 92.710407 seconds.
% %         tic
%         ROIs = BK_convert_load_ROIs(subid,roipath,fspath,size(layer_boundaries,1));
%         fprintf('Size of different ROIs:%i\n', [sum(ROIs,1)])
% %         fprintf('time for converting:')
%         %~4sec
% %         toc
%         %transform boundaries to matlab space (from freesurfer space), the
%         %resulting space will be in "voxel coordinate". It also outputs the 
% %         tic
%         [layer_boundaries,T1_mat] = VPF_transform_layers_to_matlab_space(layer_boundaries,T1path);
% %         fprintf('time for transforming:')
%         %~0sec
% %         toc
%         %sample layers
%         %todo check the code from here:
%         fprintf(sprintf('Starting layer sampling...\n'));
% %         tic
%         [columnspecificts,layeractivation,allroicolumnsize]= BK_select_active_columns_basedonts(subid,subpath,layer_boundaries,ROIs,N_layers, T1_mat);
% %         fprintf('time for layer sampling:')
%         %~500sec
% %         toc
%         interimdata_columns=struct('columnspecificts',columnspecificts,'layeractivation',layeractivation,'allroicolumnsize',allroicolumnsize);
%         save([outputpath 'interimdata_rwls.mat'],'interimdata_columns');
%         tic
        %define inputs
        %% S1 stuff
        if strcmp(region,'S1')
            load([outputpath 'interimdata_rwls.mat'],'interimdata_columns')
            columnspecificts=interimdata_columns.columnspecificts;
            sz=size(columnspecificts);
            if sz(3)~=1221
                warning('The ts is shorter!!! why?')
            end
            layeractivation=interimdata_columns.layeractivation;
            columnsize=interimdata_columns.allroicolumnsize;
            %-----------------------------------------------------------------------------------------------%
            %check if the profile depends on the way we select the voxels:
%             wholecolumnactivation=mean(layeractivation(:,10:20,:,:),2);
%             [ncol,~,~,tslengt]=size(layeractivation);
%             wholecolumnactivationts=reshape(wholecolumnactivation,ncol,1,tslengt);%size(muki)
%             [T,Tcrit,beta,pmax] = BK_column_analysis_stats(wholecolumnactivationts,subid,subpath);
            %-----------------------------------------------------------------------------------------------%
            [T,Tcrit,beta,pmax] = BK_column_analysis_stats(columnspecificts,subid,subpath);
            T(:,:,3)=min(T(:,:,:),[],3); %conjuction analysis columnwise. In the third dimension, 1-cognition(lc-hc), 2-
            columnwisestat = struct('beta',beta,'T',T,'T_crit',Tcrit,'p_max',pmax);
    
    %         interimdata_columns=struct('columnwisestat',columnwisestat,'columnspecificts',columnspecificts,'layeractivation',layeractivation,'allroicolumnsize',allroicolumnsize);
    %         save([outputpath 'interimdata_rwls.mat'],'interimdata_columns');
            
    
            numberofcolumnsinmask=length(columnsize);
            tthrsss=[0 1];
            thrsidx=0;
            for tthrs=tthrsss
                thrsidx=thrsidx+1;
                numberofcognactivcolumn(thrsidx)=sum(columnwisestat.T(:,1,1)>tthrs);
                
                numberofpainactivecolumns(thrsidx)=sum(columnwisestat.T(:,1,2)>tthrs);
                
                numberofactcolumnsbasedonconj(thrsidx)=sum(columnwisestat.T(:,1,3)>tthrs);
                activevoxelsid=find((columnwisestat.T(:,1,3)>tthrs));
    %             disp(['Subject:' subid 'The number of columns within the funcitonal mask: ',num2str(numberofcolumnsinmask)])
    %             disp(['Subject:' subid 'The number of columns acitvating in each run(t>', num2str(tthrs), ') for cognition: ',num2str(numberofcognactivcolumn(thrsidx))])
    %             disp(['Subject:' subid 'The number of columns acitvating in each run(t>', num2str(tthrs), ') for pain: ',num2str(numberofpainactivecolumns(thrsidx))])
    %             disp(['Subject:' subid 'The number of columns acitvating in each run(t>', num2str(tthrs), ') for cognition AND pain: ',num2str(numberofactcolumnsbasedonconj(thrsidx))])
            end
            
            columndistribution = struct('numberofcognactivcolumn',numberofcognactivcolumn, ...
                                        'numberofpainactivecolumns',numberofpainactivecolumns, ...
                                        'numberofactcolumnsbasedonconj',numberofactcolumnsbasedonconj, ...
                                        'numberofcolumnsinmask',numberofcolumnsinmask);
    
    %         estimate model params of active columns
            topcolumnnumber=200; %200
           
            [top200Values, top200Indices] =maxk(columnwisestat.T(:,1,3),topcolumnnumber);
    
            % Filter the indices where the values are positive
            positiveIndices = top200Indices(top200Values > 0);
            % Check if any negative values were removed
            if length(positiveIndices) < topcolumnnumber
                warning('The top 200 values contained negative numbers. Only %d positive values are retained.', length(positiveIndices));
            end
        
    %         tic
            sz=size(layeractivation);
            layerts_significant=mean(layeractivation(positiveIndices,:,:,:),1,'omitnan'); 
            layerts_significant_forfunc=reshape(layerts_significant,[sz(2), sz(3),sz(4)]);
            [T,Tcrit,beta,pmax] = BK_layer_analysis_stats(layerts_significant_forfunc,subid,subpath);
            layerwisestat = struct('beta',beta,'T',T,'T_crit',Tcrit,'p_max',pmax);
%         layerwisestat = [];
        %% pIns stuff
        elseif strcmp(region,'pIns')
            load([outputpath 'interimdata_rwls_pIns.mat'],'interimdata_columns')
            columnspecificts=interimdata_columns{1};
            sz=size(columnspecificts{1,1});
            if sz(2)~=1221
                warning('The ts is shorter!!! why?')
            end
            layeractivation=interimdata_columns{2};
            columnsize=interimdata_columns{3};
            [T,Tcrit,beta,pmax] = BK_column_analysis_stats_roicell(columnspecificts,subid,subpath);
            T(:,:,3)=min(T(:,:,:),[],3); %conjuction analysis columnwise. In the third dimension, 1-cognition(lc-hc), 2-
            columnwisestat = struct('beta',beta,'T',T,'T_crit',Tcrit,'p_max',pmax);
    
    %         interimdata_columns=struct('columnwisestat',columnwisestat,'columnspecificts',columnspecificts,'layeractivation',layeractivation,'allroicolumnsize',allroicolumnsize);
    %         save([outputpath 'interimdata_rwls.mat'],'interimdata_columns');
            
            for ROI=1:length(layeractivation)
                numberofcolumnsinmask=length(columnsize{ROI,1});
                tthrsss=[0 1];
                thrsidx=0;
                for tthrs=tthrsss
                    thrsidx=thrsidx+1;
                    numberofcognactivcolumn(thrsidx)=sum(columnwisestat.T(:,ROI,1)>tthrs);
                    
                    numberofpainactivecolumns(thrsidx)=sum(columnwisestat.T(:,ROI,2)>tthrs);
                    
                    numberofactcolumnsbasedonconj(thrsidx)=sum(columnwisestat.T(:,ROI,3)>tthrs);
                    activevoxelsid=find((columnwisestat.T(:,1,3)>tthrs));
                    disp(['Subject:' subid 'The number of columns within the funcitonal mask: ',num2str(numberofcolumnsinmask)])
                    disp(['Subject:' subid 'The number of columns acitvating in each run(t>', num2str(tthrs), ') for cognition: ',num2str(numberofcognactivcolumn(thrsidx))])
                    disp(['Subject:' subid 'The number of columns acitvating in each run(t>', num2str(tthrs), ') for pain: ',num2str(numberofpainactivecolumns(thrsidx))])
                    disp(['Subject:' subid 'The number of columns acitvating in each run(t>', num2str(tthrs), ') for cognition AND pain: ',num2str(numberofactcolumnsbasedonconj(thrsidx))])
                end
                
                columndistribution{ROI} = struct('numberofcognactivcolumn',numberofcognactivcolumn, ...
                                            'numberofpainactivecolumns',numberofpainactivecolumns, ...
                                            'numberofactcolumnsbasedonconj',numberofactcolumnsbasedonconj, ...
                                            'numberofcolumnsinmask',numberofcolumnsinmask);
            
        %         estimate model params of active columns
                topcolumnnumber=200; %200 or 100?
               
                [top200Values, top200Indices] =maxk(columnwisestat.T(:,ROI,3),topcolumnnumber);
        
                % Filter the indices where the values are positive
                positiveIndices = top200Indices(top200Values > 0);
                % Check if any negative values were removed
                if length(positiveIndices) < topcolumnnumber
                    warning('The top 200 values contained negative numbers. Only %d positive values are retained.', length(positiveIndices));
                end
                layerts=layeractivation{ROI,1};
        %         tic
                sz=size(layerts);
                layerts_significant=mean(layerts(positiveIndices,:,:),1); %,'omitnan'
                layerts_significant=squeeze(layerts_significant);
        %         layerts_significant_forfunc=reshape(layerts_significant,[sz(2), sz(3),sz(4)]);
                [T,Tcrit,beta,pmax] = BK_layer_analysis_stats_roicell(layerts_significant,subid,subpath);
                layerwisestat{ROI} = struct('beta',beta,'T',T,'T_crit',Tcrit,'p_max',pmax);
        %         layerwisestat = [];
            end
    
        
            fprintf('time for laminar activation:')
            toc
        elseif strcmp(region,'S2')
            load([outputpath 'interimdata_rwls_S2.mat'],'interimdata_columns')
            columnspecificts=interimdata_columns{1};
            sz=size(columnspecificts{1,1});
            if sz(2)~=1221
                warning('The ts is shorter!!! why?')
            end
            layeractivation=interimdata_columns{2};
            columnsize=interimdata_columns{3};
            [T,Tcrit,beta,pmax] = BK_column_analysis_stats_roicell(columnspecificts,subid,subpath);
            T(:,:,3)=min(T(:,:,:),[],3); %conjuction analysis columnwise. In the third dimension, 1-cognition(lc-hc), 2-
            columnwisestat = struct('beta',beta,'T',T,'T_crit',Tcrit,'p_max',pmax);
    
    %         interimdata_columns=struct('columnwisestat',columnwisestat,'columnspecificts',columnspecificts,'layeractivation',layeractivation,'allroicolumnsize',allroicolumnsize);
    %         save([outputpath 'interimdata_rwls.mat'],'interimdata_columns');
            
            for ROI=1:length(layeractivation)
                numberofcolumnsinmask=length(columnsize{ROI,1});
                tthrsss=[0 1];
                thrsidx=0;
                for tthrs=tthrsss
                    thrsidx=thrsidx+1;
                    numberofcognactivcolumn(thrsidx)=sum(columnwisestat.T(:,ROI,1)>tthrs);
                    
                    numberofpainactivecolumns(thrsidx)=sum(columnwisestat.T(:,ROI,2)>tthrs);
                    
                    numberofactcolumnsbasedonconj(thrsidx)=sum(columnwisestat.T(:,ROI,3)>tthrs);
                    activevoxelsid=find((columnwisestat.T(:,1,3)>tthrs));
                    disp(['Subject:' subid 'The number of columns within the funcitonal mask: ',num2str(numberofcolumnsinmask)])
                    disp(['Subject:' subid 'The number of columns acitvating in each run(t>', num2str(tthrs), ') for cognition: ',num2str(numberofcognactivcolumn(thrsidx))])
                    disp(['Subject:' subid 'The number of columns acitvating in each run(t>', num2str(tthrs), ') for pain: ',num2str(numberofpainactivecolumns(thrsidx))])
                    disp(['Subject:' subid 'The number of columns acitvating in each run(t>', num2str(tthrs), ') for cognition AND pain: ',num2str(numberofactcolumnsbasedonconj(thrsidx))])
                end
                
                columndistribution{ROI} = struct('numberofcognactivcolumn',numberofcognactivcolumn, ...
                                            'numberofpainactivecolumns',numberofpainactivecolumns, ...
                                            'numberofactcolumnsbasedonconj',numberofactcolumnsbasedonconj, ...
                                            'numberofcolumnsinmask',numberofcolumnsinmask);
            
        %         estimate model params of active columns
                topcolumnnumber=200; %200 or 100?
               
                [top200Values, top200Indices] =maxk(columnwisestat.T(:,ROI,3),topcolumnnumber);
        
                % Filter the indices where the values are positive
                positiveIndices = top200Indices(top200Values > 0);
                % Check if any negative values were removed
                if length(positiveIndices) < topcolumnnumber
                    warning('The top 200 values contained negative numbers. Only %d positive values are retained.', length(positiveIndices));
                end
                layerts=layeractivation{ROI,1};
        %         tic
                sz=size(layerts);
                layerts_significant=mean(layerts(positiveIndices,:,:),1); %,'omitnan'
                layerts_significant=squeeze(layerts_significant);
        %         layerts_significant_forfunc=reshape(layerts_significant,[sz(2), sz(3),sz(4)]);
                [T,Tcrit,beta,pmax] = BK_layer_analysis_stats_roicell(layerts_significant,subid,subpath);
                layerwisestat{ROI} = struct('beta',beta,'T',T,'T_crit',Tcrit,'p_max',pmax);
        %         layerwisestat = [];
            end
    
        
            fprintf('time for laminar activation:')
            toc
        end
        
%         interimdata=struct('columnwisestat',columnwisestat,'layerwisestat',layerwisestat,'allroicolumnsize',allroicolumnsize,'activevoxelids',activevoxelsid);
        
%         save([outputpath 'interimdata_rwls.mat'],'interimdata');
        if visualize
            layer_boundaries = VPF_load_layer_boundaries(subid,fspath);
            [layer_boundaries,~] = VPF_transform_layers_to_matlab_space(layer_boundaries,T1path);
            ROIs = BK_convert_load_ROIs(subid,roipath,fspath,size(layer_boundaries,1));
%             load([outputpath 'interimdata_rwls.mat']);
%             activevoxelsid=interimdata.activevoxelids;
            BK_plotactiveclusterscolumn(T1path,ROIs,activevoxelsid,layer_boundaries,outputpath)
        end

%     elseif exist([outputpath 'interimdata_rwls.mat'],'file') && visualize
%         %% visualize activated columns in matlab:
%         layer_boundaries = VPF_load_layer_boundaries(subid,fspath);
%         [layer_boundaries,~] = VPF_transform_layers_to_matlab_space(layer_boundaries,T1path);
%         ROIs = BK_convert_load_ROIs(subid,roipath,fspath,size(layer_boundaries,1));
%         load([outputpath 'interimdata_rwls.mat']);
%         activevoxelsid=interimdata.activevoxelids;
%         BK_plotactiveclusterscolumn(T1path,ROIs,activevoxelsid,layer_boundaries,outputpath)
%     else exist([outputpath 'interimdata_rwls.mat'],'file')
%         warning('The interimfile does already exist, and no visualization was selected. Load and play!')
%     end

end
%% load layers function 
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
    layer_boundaries = read_surf(fn);
    fn = [fspath '\' subid '\surf\rh.white'];
    tmp = read_surf(fn);
    layer_boundaries  = cat(1,layer_boundaries ,tmp);
    
    
    fn = [fspath '\' subid '\surf\lh.pial'];
    tmp = read_surf(fn);
    fn = [fspath '\' subid '\surf\rh.pial'];
    tmp = cat(1,tmp,read_surf(fn));
    layer_boundaries = cat(3,layer_boundaries,tmp);

end
%% convert the ROI to surface (wb_command) - main function
function ROIs = BK_convert_load_ROIs(subid,roipath,fspath,N_vertex)
    % This is a modification of VPF_convert_load_ROIs function. It converts the
    % ROIs to surfaces and loads them. As luckily it was already run by Viktor,
    % the freesurfer part can be commented out (but keep them to see the whole
    % picture), bc the output files are available. The checkfile function test
    % if the file is available, and output a message in case of its missing.
    %INPUT:
    %   subid [str]                : subject id
    %   subpath [str]              : path to the derivatives folder
    %   fspath [str]               : path to the freesurfer folder
    %
    %
    %OUTPUT:
    %ROIs [N_vertex,N_ROI] : 2D logical matrix containing the ROIs in surface
    %                        space. ROIs are sorted the following way... TODO
    
    % fs_base = ['export FREESURFER_HOME=/usr/local/freesurfer/6.0.0; ' ...
    %     'export SUBJECTS_DIR=' fspath '; '...
    %     'source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
    
    %localizer derived from smoothed group average and transformed back to the indiviudal space after intersecting with some anatomical mask.
    
    localizer = dir(fullfile([roipath subid '\functionalmasks\'], '**','roi_*.nii'));
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
%% transform between fs and matlab (voxel) spaces.
%this is important as the height of the columns can be calculate from here
function [layer_boundaries, T1_mat] = VPF_transform_layers_to_matlab_space(layer_boundaries,T1path)
    %Transforms the layers from freesurfer space to matlab space. After that, 
    %signal can be sampled using spm_sample_vol
    %INPUT:
    %layer_boundaries [vertex,[x,y,z],[white,pial]] : 3D matrix of the surfaces
    %T1path [str]                                   : full path and filename of 
    %                                                 the T1 image
    
    %OUTPUT:
    %layer_boundaries        : 3D matrix of surfaces in matlab space
    %
    
    sz = size(layer_boundaries);
    
    %bring surfaces into matlab space
    T1_mat = spm_get_space(T1path); %Get/set the voxel-to-world mapping of an image
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
    
    % img_T1 = spm_read_vols(spm_vol(T1path));
    % for slice = 205
    %     idx = find(layer_boundaries(:,3,1)>slice & layer_boundaries(:,3,1)<slice+1);
    %     idx3 = find(layer_boundaries(:,3,2)>slice & layer_boundaries(:,3,2)<slice+1);
    % 
    %     figure;imagesc(squeeze(img_T1(:,:,slice)), [0 1200]), colormap(gray(256)), title(['Slice ' num2str(slice)])
    %     hold on;
    % 
    %     plot(layer_boundaries(idx,2,1),layer_boundaries(idx,1,1),'r.');
    %     plot(layer_boundaries(idx3,2,2),layer_boundaries(idx3,1,2),'y.');
    % 
    %     legend('white','pial');
    % end

end

%% sample columnswise information based on ROI
function [columns_out,layers_out,allroicolumnsize] = BK_select_active_columns_basedonts(subid,subpath,layer_boundaries,ROIs,N_layers, T1_mat)
% This is a modificaiton of VPF core sampling function (VPF_sample_layers).
% It does two things:
%   1. get the timeseries signal from the deep 75% of the column (unbias (less effect of venous bias) estimation of
%   columnar activation) 
%   2. get the timeseries signal from the whole column
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
% columns_out [N_columns N_ROI N_totalvol]
% layer_out [N_columns N_layers N_ROI N_totalvol]
% allroicolumnsize {N_ROI} the distribution of the 

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
    fullvol=1;
    for run=1:N_run
        fprintf(sprintf(['run ' num2str(run) '...']));
        for vol = 1:N_vols(run)
            cd([runs(run).folder '\' runs(run).name '\func']);
    %         sampled_img = load_nifti([sampled_img_list(vol).name]).vol;
            sampled_img = niftiread([sampled_img_list{run}(vol).name]);
    %         "E:\pain_layers\main_project\derivatives\pipeline\7349\ses-02\func\layers\rwls_stats_compcor_UNsmoothed_hpf180\spmT_0005.nii.gz"
    %         parfor ROI = 1:N_ROI
            for ROI = 1:N_ROI
                ind = find(ROIs(:,ROI)); %columns under the group lvl mask
                tmp = zeros(size(ind,1),N_layers);
                boundaries = layer_boundaries;
                if run==1 && ROI==1 && vol==1
                    columns_out = zeros([length(ind),N_ROI,sum(N_vols)]);
                    layers_out = zeros([length(ind),N_layers,N_ROI,sum(N_vols)]);
                    columnsize=nan(1,size(ind,1));
                end
               
                %I am not sure what is the vertex coordinate here, but I can
%         %believe that is in mm (mean+/-sd=3.55+/-1.1)
%         %double check if that is the case. 
%         % think about to use only columns which are within a range, which
%         % suggest that the definition of that column makes sense(too big or
%         % too small is OK??)
    %we go here from vertex to vertex(==column to column), and sample the funcitonal images. Here we should include an additional function, to select only columns active mostly in the deep layer. However, this need to be happen on the unsmoothed main effects (smoothed actually is available) conjuction map. That is calculate indiviudal conjuction maps, and apply the mask on those and estimate the layer profile. However, we would need to have a decent t value at the deeper part of the cortex. (my gut feeleing that interpolate statistical map is not that correct.-think about it a bit more)
                for k = 1:size(ind,1) 
                    X  = linspace(boundaries(ind(k),1,1),boundaries(ind(k),1,2),N_layers);
                    Y  = linspace(boundaries(ind(k),2,1),boundaries(ind(k),2,2),N_layers);
                    Z  = linspace(boundaries(ind(k),3,1),boundaries(ind(k),3,2),N_layers);
                    if vol==1 && run==1 %calculate the column height. This need to be calculated only once.
                        columnsize(k)=sqrt((X(end)-X(1))^2+(Y(end)-Y(1))^2+(Z(end)-Z(1))^2) * T1_mat(1,1); 
                        %it seems that as we are in matlab space, we need
                        %to transform the distance between the two
                        %points (column bottom and top) back to world
                        %space. This can be done simply applyint the affine
                        %matrix but as it is the same for all subjects and
                        %isovoxel resolution with a 0.75mm, we can simple
                        %multiply the distance with this constant. I
                        %hardcoded it!!!
                    end
                    tmp(k,:) = spm_sample_vol(sampled_img,X,Y,Z,SAMPLING_ORDER);
                end
                if vol==1 && run==1
                    allroicolumnsize{ROI}=columnsize;
                end
    %             layers_out(:,ROI,vol,run) = squeeze(mean(tmp,1,'omitnan'));
                layers_out(:,:,ROI,fullvol) = tmp;
                columns_out(:,ROI,fullvol) = mean(tmp(:,1:15),2,'omitnan'); %15 is hardcoded here as Peter recommended to skip the top 0.75mm, which would be in case of 3mm avg thickness the same as 15/20.It can be adjusted later by the thickness of the column as we can easily calculate that.
                fullvol=fullvol+1;
            end
        end
        toc
        fprintf(sprintf(' done \n'));
    end
cd(currpath)
end

%% estimate 1st lvl GLM
% column wise
function [T,T_crit,beta,p_max]=BK_column_analysis_stats(columnspecificts,subid,subpath,ZTRANS)%subpath,ZTRANS)

    if nargin < 4
        ZTRANS = false;
    end
    % layerpath = [subpath '/ses-02/func/layers'];
    % load([layerpath '/layers.mat'],'layers');
    
    [N_columns,N_ROIS,~] = size(columnspecificts);
    contrasts = [4 8]; %main effects --> calculate later some conjunction stuff?, I would go now with the average effect, and not with the individual run effects.
    N_contrasts=length(contrasts);
    % T = zeros(N_columns,N_ROIS,N_runs,N_contrasts,2); % 2 for no_compcor vs compcor
    % beta = zeros(N_columns,N_ROIS,N_runs,N_contrasts,2);
    % T_crit = zeros(N_ROIS,N_runs,N_contrasts,2);
    % p_max = zeros(N_ROIS,N_runs,N_contrasts,2);
    
    %maybe implement 
    T = zeros(N_columns,N_ROIS,N_contrasts); 
    beta = zeros(N_columns,N_ROIS,N_contrasts);
    T_crit = zeros(N_ROIS,N_contrasts);
    p_max = zeros(N_ROIS,N_contrasts);
    % for run = 1:N_runs
%     statspath = dir(['E:\pain_layers\main_project\derivatives\pipeline\' num2str(subid) '/ses-02/func/layers/*_stats_compcor_UNsmoothed_hpf180']);
    statspath = dir([subpath num2str(subid) '/ses-02/func/layers/rwls_stats_compcor_UNsmoothed_hpf180']);
%     statspath = dir([subpath num2str(subid) '/ses-02/func/layers/reml_stats_compcor_UNsmoothed_hpf180']);
    
    % "E:\pain_layers\main_project\derivatives\pipeline\7349\ses-02\func\layers\rwls_stats_compcor_UNsmoothed_hpf180"

%     for folder = 1:length(statspath)
%         load([statspath(folder).folder '/' statspath(folder).name '/SPM.mat'],'SPM');
        load([statspath(1).folder '/SPM.mat'],'SPM');
        W = SPM.xX.W;
        GLM = SPM.xX.pKX;
%         con = zeros(size(GLM,1),1);
        for ROI = 1:N_ROIS
            Y = squeeze(columnspecificts(:,ROI,:)); %in the original pipeline, the size of Y is layernnum X tslength.
            Y = Y(:,any(Y ~= 0, 1)); %selects only the columns where at least one element is non-zero, effectively removing all-zero columns from Y
            if ZTRANS
                %baseline z-transform. I take the volumes corresponding
                % to 0 in the sum of all pain trials as baseline
                idx = find(sum(SPM.xX.X(:,3:22),2)==0);
                m = mean(Y(:,idx),2);
                s = std(Y(:,idx),[],2);
                Y = (Y-m)./s;
            end
    
            KWY = spm_filter(SPM.xX.K,W*Y.');
            b   = GLM*KWY;
    
            res      = spm_sp('r',SPM.xX.xKXs,KWY);        %-Residuals
            ResSS    = sum(res.^2);                    %-Residual SSQ
            ResMS = ResSS / SPM.xX.trRV;
    
            for idxcontastofinterest = 1:numel(contrasts)
                %We always assume a one-sided effect, i.e. for ROIs where
                %the localizer showed a negative effect, we assume that mu < 0.
                
%                 con(con==1) = 0;
%                 con(ii+2) = 1;
                contastofinterest=contrasts(idxcontastofinterest);
                if contrasts(idxcontastofinterest)<5
                    reversefact=-1; %this would mean that the activation is lower in the high cognition condition.
                else
                    reversefact=1;
                end
                %     {'anticipation_high_cognition'}
                %     {'anticipation_low_cognition' }
                %     {'pain_high_cogn_high_pain'   }
                %     {'pain_high_cogn_low_pain'    }
                %     {'pain_low_cogn_high_pain'    }
                %     {'pain_low_cogn_low_pain'     }
                %     {'rating'                     }
                con=SPM.xCon(contastofinterest).c*reversefact; 
                [T(:,ROI,idxcontastofinterest),...
                 T_crit(ROI,idxcontastofinterest),...
                 beta(:,ROI,idxcontastofinterest),... %this should be a cell array instead of a matrix as the number of columns iwhtin a ROI most porbably change.
                 p_max(ROI,idxcontastofinterest)] = BK_Tmap_from_SPM_columns(SPM,b,ResMS,con,0.05,'FDR');
            end
        end
%     end
    % end
    
%     rwls_results = struct('beta',beta,'T',T,'T_crit',T_crit,'p_max',p_max);
%     save([layerpath '/rwls_results.mat'],'rwls_results');
end

% column wise - multiple ROI as cell
function [T,T_crit,beta,p_max]=BK_column_analysis_stats_roicell(columnspecificts,subid,subpath,ZTRANS)%subpath,ZTRANS)

    if nargin < 4
        ZTRANS = false;
    end
    % layerpath = [subpath '/ses-02/func/layers'];
    % load([layerpath '/layers.mat'],'layers');
    
%     [N_columns,N_ROIS,~] = size(columnspecificts);
    N_ROIS=size(columnspecificts,1);
    contrasts = [4 8]; %main effects --> calculate later come conjunction stuff?, I would go now with the average effect, and not with the individual run effects.
    N_contrasts=length(contrasts);
    for roi=1:N_ROIS
        columnsss(roi)=size(columnspecificts{roi,1},1);
    end
    N_columns=max(columnsss);
    % T = zeros(N_columns,N_ROIS,N_runs,N_contrasts,2); % 2 for no_compcor vs compcor
    % beta = zeros(N_columns,N_ROIS,N_runs,N_contrasts,2);
    % T_crit = zeros(N_ROIS,N_runs,N_contrasts,2);
    % p_max = zeros(N_ROIS,N_runs,N_contrasts,2);
    
    %maybe implement 
    
    % for run = 1:N_runs
%     statspath = dir(['E:\pain_layers\main_project\derivatives\pipeline\' num2str(subid) '/ses-02/func/layers/*_stats_compcor_UNsmoothed_hpf180']);
    statspath = dir([subpath num2str(subid) '/ses-02/func/layers/rwls_stats_compcor_UNsmoothed_hpf180']);
    % "E:\pain_layers\main_project\derivatives\pipeline\7349\ses-02\func\layers\rwls_stats_compcor_UNsmoothed_hpf180"

%     for folder = 1:length(statspath)
%         load([statspath(folder).folder '/' statspath(folder).name '/SPM.mat'],'SPM');
        load([statspath(1).folder '/SPM.mat'],'SPM');
        W = SPM.xX.W;
        GLM = SPM.xX.pKX;
        T = nan(N_columns,N_ROIS,N_contrasts); 
        beta = nan(N_columns,N_ROIS,N_contrasts);
        T_crit = nan(N_ROIS,N_contrasts);
        p_max = nan(N_ROIS,N_contrasts);
%         con = zeros(size(GLM,1),1);
        for ROI = 1:N_ROIS
            
            Y = squeeze(columnspecificts{ROI,1}); %in the original pipeline, the size of Y is layernnum X tslength.
            Y = Y(:,any(Y ~= 0, 1)); %selects only the columns(rows in this matrix) where at least one element is non-zero, effectively removing all-zero columns from Y
            if ZTRANS
                %baseline z-transform. I take the volumes corresponding
                % to 0 in the sum of all pain trials as baseline
                idx = find(sum(SPM.xX.X(:,3:22),2)==0);
                m = mean(Y(:,idx),2);
                s = std(Y(:,idx),[],2);
                Y = (Y-m)./s;
            end
    
            KWY = spm_filter(SPM.xX.K,W*Y.');
            b   = GLM*KWY;
    
            res      = spm_sp('r',SPM.xX.xKXs,KWY);        %-Residuals
            ResSS    = sum(res.^2);                    %-Residual SSQ
            ResMS = ResSS / SPM.xX.trRV;
    
            for idxcontastofinterest = 1:numel(contrasts)
                %We always assume a one-sided effect, i.e. for ROIs where
                %the localizer showed a negative effect, we assume that mu < 0.
                
%                 con(con==1) = 0;
%                 con(ii+2) = 1;
                contastofinterest=contrasts(idxcontastofinterest);
                if contrasts(idxcontastofinterest)<5
                    reversefact=-1; %this would mean that the activation is lower in the high cognition condition.
                else
                    reversefact=1;
                end
                %     {'anticipation_high_cognition'}
                %     {'anticipation_low_cognition' }
                %     {'pain_high_cogn_high_pain'   }
                %     {'pain_high_cogn_low_pain'    }
                %     {'pain_low_cogn_high_pain'    }
                %     {'pain_low_cogn_low_pain'     }
                %     {'rating'                     }
                con=SPM.xCon(contastofinterest).c*reversefact; 
                [T(1:columnsss(ROI),ROI,idxcontastofinterest),...
                 T_crit(ROI,idxcontastofinterest),...
                 beta(1:columnsss(ROI),ROI,idxcontastofinterest),... %this should be a cell array instead of a matrix as the number of columns iwhtin a ROI most porbably change.
                 p_max(ROI,idxcontastofinterest)] = BK_Tmap_from_SPM_columns(SPM,b,ResMS,con,0.05,'FDR');
            end
        end
%     end
    % end
    
%     rwls_results = struct('beta',beta,'T',T,'T_crit',T_crit,'p_max',p_max);
%     save([layerpath '/rwls_results.mat'],'rwls_results');
end
% layerwise
function [T,T_crit,beta,p_max]=BK_layer_analysis_stats(laminarts,subid,subpath,ZTRANS)%subpath,ZTRANS)
    if nargin < 4
        ZTRANS = false;
    end
    % layerpath = [subpath '/ses-02/func/layers'];
    % load([layerpath '/layers.mat'],'layers');
    
    [N_layer,N_ROIS,~] = size(laminarts);
    %we hardcode the used contrast here. this is something which we set,
    %and related to the design:
    %     {'anticipation_high_cognition'}
    %     {'anticipation_low_cognition' }
    %     {'pain_high_cogn_high_pain'   }
    %     {'pain_high_cogn_low_pain'    }
    %     {'pain_low_cogn_high_pain'    }
    %     {'pain_low_cogn_low_pain'     }
    %     {'rating'                     }
%     contrasts = [4 8]; %main effects --> calculate later some conjunction stuff?, I would go now with the average effect, and not with the individual run effects.
%     N_contrasts=length(contrasts);
    contrasts = [3 4 5 6]; %condition effects. to model every condition separately
%     contrasts = [3 5]; %condition effects.
    numberofregressors=18; %this is the total number of regression which used in the 1st lvl GLM. We have 7 regressor to model (see above) the design, and an additional 11 regressor of noise (6 motion + 5 WM compcor).
    N_contrasts=length(contrasts);
    % T = zeros(N_columns,N_ROIS,N_runs,N_contrasts,2); % 2 for no_compcor vs compcor
    % beta = zeros(N_columns,N_ROIS,N_runs,N_contrasts,2);
    % T_crit = zeros(N_ROIS,N_runs,N_contrasts,2);
    % p_max = zeros(N_ROIS,N_runs,N_contrasts,2);
    
    %maybe implement 
    T = zeros(N_layer,N_ROIS,N_contrasts); 
    beta = zeros(N_layer,N_ROIS,N_contrasts);
    T_crit = zeros(N_ROIS,N_contrasts);
    p_max = zeros(N_ROIS,N_contrasts);
    % for run = 1:N_runs
%     statspath = dir(['E:\pain_layers\main_project\derivatives\pipeline\' num2str(subid) '/ses-02/func/layers/*_stats_compcor_UNsmoothed_hpf180']);
%     statspath = dir(['E:\pain_layers\main_project\derivatives\pipeline\' num2str(subid) '/ses-02/func/layers/rwls_stats_compcor_UNsmoothed_hpf180']);
    statspath = dir([subpath num2str(subid) '/ses-02/func/layers/rwls_stats_compcor_UNsmoothed_hpf180']);
%     statspath = dir([subpath num2str(subid) '/ses-02/func/layers/reml_stats_compcor_UNsmoothed_hpf180']);
%     statspath = dir(['E:\pain_layers\main_project\derivatives\pipeline\'
%     num2str(subid) '/ses-02/func/layers/*_stats_compcor_UNsmoothed_hpf180']); c
    % "E:\pain_layers\main_project\derivatives\pipeline\7349\ses-02\func\layers\rwls_stats_compcor_UNsmoothed_hpf180"

%     for folder = 1:length(statspath)
%         disp([statspath(folder).folder '/' statspath(folder).name]);
%         load([statspath(folder).folder '/' statspath(folder).name '/SPM.mat'],'SPM');
        load([statspath(1).folder '/SPM.mat'],'SPM');
    
        W = SPM.xX.W;
        GLM = SPM.xX.pKX;
%         con = zeros(size(GLM,1),1);
        for ROI = 1:N_ROIS
            Y = squeeze(laminarts(:,ROI,:)); %in the original pipeline, the size of Y is layernnum X tslength.
            Y = Y(:,any(Y ~= 0, 1)); %selects only the columns where at least one element is non-zero, effectively removing all-zero columns from Y
            if ZTRANS
                %baseline z-transform. I take the volumes corresponding
                % to 0 in the sum of all pain trials as baseline
                idx = find(sum(SPM.xX.X(:,3:22),2)==0);
                m = mean(Y(:,idx),2);
                s = std(Y(:,idx),[],2);
                Y = (Y-m)./s;
            end
    
            KWY = spm_filter(SPM.xX.K,W*Y.');
            b   = GLM*KWY;
    
            res      = spm_sp('r',SPM.xX.xKXs,KWY);        %-Residuals
            ResSS    = sum(res.^2);                    %-Residual SSQ
            ResMS = ResSS / SPM.xX.trRV;
    

%             for idxcontastofinterest = 1:numel(contrasts)
%                 %We always assume a one-sided effect, i.e. for ROIs where
%                 %the localizer showed a negative effect, we assume that mu < 0.
%                 
% %                 con(con==1) = 0;
% %                 con(ii+2) = 1;
%                 contastofinterest=contrasts(idxcontastofinterest);
%                 if contrasts(idxcontastofinterest)<5
%                     reversefact=-1; %this would mean that the activation is lower in the high cognition condition.
%                 else
%                     reversefact=1;
%                 end
%                 %     {'anticipation_high_cognition'}
%                 %     {'anticipation_low_cognition' }
%                 %     {'pain_high_cogn_high_pain'   }
%                 %     {'pain_high_cogn_low_pain'    }
%                 %     {'pain_low_cogn_high_pain'    }
%                 %     {'pain_low_cogn_low_pain'     }
%                 %     {'rating'                     }
%                 con=SPM.xCon(contastofinterest).c*reversefact;
            idxcontastofinterest=0;
            for k=contrasts
                idxcontastofinterest=idxcontastofinterest+1;
                 
                con=zeros(57,1); %hardcoded here, but can be changed to the length(SPM.xCon(1).c)
                con([k, k+numberofregressors,k+(numberofregressors*2)],:)=1; %each condition separately
%                 con([k, k+numberofregressors, k+(numberofregressors*2)],:)=1;
%                 con([k+1, k+1+numberofregressors, k+1+(numberofregressors*2)],:)=-1;
%                 disp(con)
    
                
                [T(:,ROI,idxcontastofinterest),...
                 T_crit(ROI,idxcontastofinterest),...
                 beta(:,ROI,idxcontastofinterest),... %this should be a cell array instead of a matrix as the number of columns iwhtin a ROI most porbably change.
                 p_max(ROI,idxcontastofinterest)] = BK_Tmap_from_SPM_columns(SPM,b,ResMS,con,0.05,'FDR');
            end
        end
%     end
    % end
    
%     rwls_results = struct('beta',beta,'T',T,'T_crit',T_crit,'p_max',p_max);
%     save([layerpath '/rwls_results.mat'],'rwls_results');
end
%layerwise cell roi
function [T,T_crit,beta,p_max]=BK_layer_analysis_stats_roicell(laminarts,subid,subpath,ZTRANS)%subpath,ZTRANS)
    if nargin < 4
        ZTRANS = false;
    end
    % layerpath = [subpath '/ses-02/func/layers'];
    % load([layerpath '/layers.mat'],'layers');
    
    [N_layer,~] = size(laminarts);
    %we hardcode the used contrast here. this is something which we set,
    %and related to the design:
    %     {'anticipation_high_cognition'}
    %     {'anticipation_low_cognition' }
    %     {'pain_high_cogn_high_pain'   }
    %     {'pain_high_cogn_low_pain'    }
    %     {'pain_low_cogn_high_pain'    }
    %     {'pain_low_cogn_low_pain'     }
    %     {'rating'                     }
%     contrasts = [3 4 5 6]; %condition effects. to model every condition separately
%     contrasts = [3 5]; %condition effects.
%     numberofregressors=18; %this is the total number of regression which used in the 1st lvl GLM. We have 7 regressor to model (see above) the design, and an additional 11 regressor of noise (6 motion + 5 WM compcor).
    contrasts = [4 8];
    N_contrasts=length(contrasts);
    % T = zeros(N_columns,N_ROIS,N_runs,N_contrasts,2); % 2 for no_compcor vs compcor
    % beta = zeros(N_columns,N_ROIS,N_runs,N_contrasts,2);
    % T_crit = zeros(N_ROIS,N_runs,N_contrasts,2);
    % p_max = zeros(N_ROIS,N_runs,N_contrasts,2);
    
    %maybe implement 
    T = zeros(N_layer,N_contrasts); 
    beta = zeros(N_layer,N_contrasts);
    T_crit = zeros(N_contrasts,1);
    p_max = zeros(N_contrasts,1);
    % for run = 1:N_runs
%     statspath = dir(['E:\pain_layers\main_project\derivatives\pipeline\' num2str(subid) '/ses-02/func/layers/*_stats_compcor_UNsmoothed_hpf180']);
%     statspath = dir(['E:\pain_layers\main_project\derivatives\pipeline\' num2str(subid) '/ses-02/func/layers/rwls_stats_compcor_UNsmoothed_hpf180']);
    statspath = dir([subpath num2str(subid) '/ses-02/func/layers/rwls_stats_compcor_UNsmoothed_hpf180']);
%     statspath = dir(['E:\pain_layers\main_project\derivatives\pipeline\'
%     num2str(subid) '/ses-02/func/layers/*_stats_compcor_UNsmoothed_hpf180']); c
    % "E:\pain_layers\main_project\derivatives\pipeline\7349\ses-02\func\layers\rwls_stats_compcor_UNsmoothed_hpf180"

%     for folder = 1:length(statspath)
%         disp([statspath(folder).folder '/' statspath(folder).name]);
%         load([statspath(folder).folder '/' statspath(folder).name '/SPM.mat'],'SPM');
        load([statspath(1).folder '/SPM.mat'],'SPM');
    
        W = SPM.xX.W;
        GLM = SPM.xX.pKX;
%         con = zeros(size(GLM,1),1);
%         for ROI = 1:N_ROIS
            Y = squeeze(laminarts(:,:)); %in the original pipeline, the size of Y is layernnum X tslength.
            Y = Y(:,any(Y ~= 0, 1)); %selects only the columns where at least one element is non-zero, effectively removing all-zero columns from Y
            if ZTRANS
                %baseline z-transform. I take the volumes corresponding
                % to 0 in the sum of all pain trials as baseline
                idx = find(sum(SPM.xX.X(:,3:22),2)==0);
                m = mean(Y(:,idx),2);
                s = std(Y(:,idx),[],2);
                Y = (Y-m)./s;
            end
    
            KWY = spm_filter(SPM.xX.K,W*Y.');
            b   = GLM*KWY;
    
            res      = spm_sp('r',SPM.xX.xKXs,KWY);        %-Residuals
            ResSS    = sum(res.^2);                    %-Residual SSQ
            ResMS = ResSS / SPM.xX.trRV;
    
%             idxcontastofinterest=0;
            for idxcontastofinterest = 1:numel(contrasts)
                contastofinterest=contrasts(idxcontastofinterest);
                if contrasts(idxcontastofinterest)<5
                    reversefact=-1; %this would mean that the activation is lower in the high cognition condition.
                else
                    reversefact=1;
                end
                %     {'anticipation_high_cognition'}
                %     {'anticipation_low_cognition' }
                %     {'pain_high_cogn_high_pain'   }
                %     {'pain_high_cogn_low_pain'    }
                %     {'pain_low_cogn_high_pain'    }
                %     {'pain_low_cogn_low_pain'     }
                %     {'rating'                     }
                con=SPM.xCon(contastofinterest).c*reversefact;
%             for k=contrasts
%                 idxcontastofinterest=idxcontastofinterest+1;
%                  
%                 con=zeros(57,1); %hardcoded here, but can be chenged to the length(SPM.xCon(1).c)
%                 con([k, k+numberofregressors,k+(numberofregressors*2)],:)=1; %each condition separately
%                 con([k, k+numberofregressors, k+(numberofregressors*2)],:)=1;
%                 con([k+1, k+1+numberofregressors, k+1+(numberofregressors*2)],:)=-1;
%                 disp(con)
    
                
                [T(:,idxcontastofinterest),...
                 T_crit(idxcontastofinterest),...
                 beta(:,idxcontastofinterest),... %this should be a cell array instead of a matrix as the number of columns iwhtin a ROI most porbably change.
                 p_max(idxcontastofinterest)] = BK_Tmap_from_SPM_columns(SPM,b,ResMS,con,0.05,'FDR');
            end
%         end
%     end
    % end
    
%     rwls_results = struct('beta',beta,'T',T,'T_crit',T_crit,'p_max',p_max);
%     save([layerpath '/rwls_results.mat'],'rwls_results');
end

%% calculate columnwise t-values
function [T,Tcrit,con_array,pmax] = BK_Tmap_from_SPM_columns(SPM,beta,ResMS,con,alpha,flag)

    if any(isnan(beta))
        T = nan(1,size(beta,2));
        Tcrit = nan;
        con_array = T;
        pmax = nan;
        return
    end
    Vc  = con'*SPM.xX.Bcov*con;
    SE  = sqrt(ResMS*Vc);    
    
    beta_index = find(abs(con) > 0);
    beta_use   = beta(beta_index,:);    
    
    con_array     = zeros(1,size(beta,2));
    for j=1:size(beta_use,1)
        con_array = con_array + con(beta_index(j)) * beta_use(j,:);
    end
    
    
    T   = con_array./SE;
    
    switch flag
        case 'FWE'
            Tcrit = spm_uc(alpha,[1 SPM.xX.erdf],'T',SPM.xVol.R,1,numel(beta));
        case 'FDR'
            p = 2 * (1 - spm_Tcdf(abs(T), SPM.xX.erdf));
            p = spm_P_FDR(p);
            Tcrit = min(abs(T(p<alpha)));
    
            if isempty(Tcrit)
                Tcrit = nan;
            end
            pmax = max(p(p<alpha));
            if isempty(pmax)
                pmax = 1;
            end
        case 'none'
            Tcrit = spm_u(alpha,[1 SPM.xX.erdf],'T');
            pmax = alpha;
    end
    
    if sum(T<0) > sum(T>0)
        Tcrit = -Tcrit;
    end


end

%% visualize active columns in matlab space
function BK_plotactiveclusterscolumn(T1path,ROIs,activevoxelsid,layer_boundaries,outputpath)
    img_T1 = spm_read_vols(spm_vol(T1path));
    figure;
    % set(gcf, 'Position', [100, 100, 800, 600]); % [x, y, width, height]
    slice_idx = 0; % To keep track of subplot index
    roimaskbycolumn=find(ROIs(:,1));
    
    visualizeddimension=2;
%     numberofimgs=range(layer_boundaries(roimaskbycolumn(activevoxelsid),visualizeddimension,1))/2;
%     initialslice=round(min(layer_boundaries(roimaskbycolumn(activevoxelsid),visualizeddimension,1)));
%     endslice=round(max(layer_boundaries(roimaskbycolumn(activevoxelsid),visualizeddimension,1)));
    slicethickness=0.5;
    zoom=50;
    [lowerinterval,~]=BK_define_intervalsforvisualization(layer_boundaries(roimaskbycolumn(activevoxelsid),visualizeddimension,1),slicethickness,4);
    
        for slice = round(lowerinterval)%round(linspace(initialslice, endslice, 4)); %initialslice:2:endslice
                % Update subplot index
            slice_idx = slice_idx + 1;
        
            % Create a subplot for each slice
            ax = subplot(2, 2, slice_idx); % 2 rows, ceil() for columns
            pos = get(ax, 'Position'); % Get current position
            pos(3) = pos(3) * 1.2;     % Increase width by 20%
            pos(4) = pos(4) * 1.2;     % Increase height by 20%
            set(ax, 'Position', pos);
            %define image parameters>:
        %     figure
        %     subplot(1,2,1)
            rotatedimg=rot90(squeeze(img_T1(:,slice,:)));
            %     layer_boundaries(numberofactcolumnsbasedonconj(4),:,:);
        %     layer_boundaries(roimaskbycolumn(activevoxels),:,1)
            
            %surface
            idx_surf = find(layer_boundaries(:,visualizeddimension,1)>=slice & layer_boundaries(:,visualizeddimension,1)<slice+slicethickness);
            idx3_surf = find(layer_boundaries(:,visualizeddimension,2)>=slice & layer_boundaries(:,visualizeddimension,2)<slice+slicethickness);
            %functional mask:
            idx_roimask = intersect(idx_surf,roimaskbycolumn);%find(layer_boundaries(roimaskbycolumn,visualizeddimension,1)>slice & layer_boundaries(roimaskbycolumn,visualizeddimension,1)<slice+1);
%             idx3_roismask = intersect(idx3_surf,roimaskbycolumn);%find(layer_boundaries(roimaskbycolumn,visualizeddimension,2)>slice & layer_boundaries(roimaskbycolumn,visualizeddimension,2)<slice+1);
            % significant voxels in the functional mask:
            idx_signfuncmask = intersect(idx_surf,roimaskbycolumn(activevoxelsid));%find(layer_boundaries(roimaskbycolumn(activevoxels),visualizeddimension,1)>slice & layer_boundaries(roimaskbycolumn(activevoxels),visualizeddimension,1)<slice+1);
%             idx3_signfuncmask = intersect(idx3_surf,roimaskbycolumn(activevoxelsid));%find(layer_boundaries(roimaskbycolumn(activevoxels),visualizeddimension,2)>slice & layer_boundaries(roimaskbycolumn(activevoxels),visualizeddimension,2)<slice+1);
            fprintf('the number of columns we visualize:%i\n',length(idx_signfuncmask))
            % do the plotting:
            %highres anat
            imagesc(rotatedimg, [0 1200]), colormap(gray(256)), title(['Slice ' num2str(slice)])
        %     axis([100 200 150 250]);
            axis(round([min(layer_boundaries(idx_roimask,1,1))-zoom  ...
                max(layer_boundaries(idx_roimask,1,1))+zoom  ...
                min(size(rotatedimg, 1)-layer_boundaries(idx_roimask,3,1))-30 ...
                max(size(rotatedimg, 1)-layer_boundaries(idx_roimask,3,1))+30]));
            hold on;
            % surface
            plot(layer_boundaries(idx_surf,1,1),size(rotatedimg, 1)-layer_boundaries(idx_surf,3,1),'r.');
            plot(layer_boundaries(idx3_surf,1,2),size(rotatedimg, 1)-layer_boundaries(idx3_surf,3,2),'y.');
        %     subplot(1,2,2)
        %     imagesc(squeeze(img_T1(:,slice,:)), [0 1200]), colormap(gray(256)), title(['Slice ' num2str(slice)])
        %     hold on
        
            %funcitonal ROImask
            plot(layer_boundaries(idx_roimask,1,1),size(rotatedimg, 1)-layer_boundaries(idx_roimask,3,1),'g.');
        %     plot(layer_boundaries(idx3_roismask,1,2),size(rotatedimg, 1)-layer_boundaries(idx3_roismask,3,2),'b.');
            %significant columns
            plot(layer_boundaries(idx_signfuncmask,1,1),size(rotatedimg, 1)-layer_boundaries(idx_signfuncmask,3,1),'b.');
        %     plot(layer_boundaries(idx3_signfuncmask,1,2),size(rotatedimg, 1)-layer_boundaries(idx3_signfuncmask,3,2),'b.');
            xlabel(''); % Clear x-axis label
            ylabel(''); % Clear y-axis label
            axis off;
            
        %     legend('white','pial');
        end
        % t = tiledlayout(2, ceil((endslice-initialslice)/4), 'Padding', 'compact', 'TileSpacing', 'compact');
        % tiledlayout(2, ceil((endslice-initialslice)/4), 'Padding', 'compact', 'TileSpacing', 'compact');
        % sgtitle(['The number of signifiantly activating columns:', num2str(length(activevoxels))])
        txtsignifcolumnnum=['The number of signifiantly activating columns:', num2str(length(activevoxelsid))];
        annotation('textbox', [0, 0.5, 0.1, 0.05], 'String', txtsignifcolumnnum, ...
                   'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                   'FontSize', 14, 'LineStyle', 'none');
        set(gcf, 'PaperUnits', 'inches', 'PaperSize', [12, 9]); % Match the figure size
        print(gcf, [outputpath 'S1_activecluster.pdf'], '-dpdf', '-bestfit');
        
        
        
    % layers = VPF_sample_layers(subid,subpath,layer_boundaries,ROIs,N_layers);
end
%% deprecated stuff
function layers_out = BK_select_active_columns(subid,subpath,layer_boundaries,ROIs,N_layers)

% This is a modificaiton of VPF core sampling function (VPF_sample_layers).
% The function aims to select columns which active within the functional ROI(smoothed group avg) based on the unsmoothed . 
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
%
%OUTPUT:
%layer_out [N_layers N_ROI N_vol N_run]
%                           : sampled layer array

currpath = pwd;
SAMPLING_ORDER = 3;
runs = dir([subpath '\' subid '\ses-02\func\layers\run*']);
N_run = numel(runs);
N_ROI = size(ROIs,2);
maki="E:\pain_layers\main_project\derivatives\pipeline\7349\ses-02\func\layers\rwls_stats_compcor_UNsmoothed_hpf180\con_0005.nii.gz";
sampled_img = niftiread(maki);
layers_out = zeros([N_layers,N_ROI,N_run]);
for ROI = 1:N_ROI
    ind = find(ROIs(:,ROI));
    tmp = zeros(size(ind,1),N_layers);
    boundaries = layer_boundaries;
%we go here from vertex to vertex(==column to column), and sample the funcitonal images. Here we should include an additional function, to select only columns active mostly in the deep layer. However, this need to be happen on the unsmoothed main effects (smoothed actually is available) conjuction map. That is calculate indiviudal conjuction maps, and apply the mask on those and estimate the layer profile. However, we would need to have a decent t value at the deeper part of the cortex. (my gut feeleing that interpolate statistical map is not that correct.-think about it a bit more)
    for k = 1:size(ind,1) 
        X  = linspace(boundaries(ind(k),1,1),boundaries(ind(k),1,2),N_layers);
        Y  = linspace(boundaries(ind(k),2,1),boundaries(ind(k),2,2),N_layers);
        Z  = linspace(boundaries(ind(k),3,1),boundaries(ind(k),3,2),N_layers);
        %I am not sure what is the vertex coordinate here, but I can
        %believe that is in mm (mean+/-sd=3.55+/-1.1)
        %double check if that is the case. 
        % think about to use only columns which are within a range, which
        % suggest that the definition of that column makes sense(too big or
        % too small is OK??)
        columnsize(k)=sqrt((X(end)-X(1))^2+(Y(end)-Y(1))^2+(Z(end)-Z(1))^2);
        tmp(k,:) = spm_sample_vol(sampled_img,X,Y,Z,SAMPLING_ORDER);
    end
    layers_out(:,ROI,vol,run) = squeeze(mean(tmp,1,'omitnan'));
end
% for run = 1:N_run
%     fprintf(sprintf(['run ' num2str(run) '...']));
% 
%     cd([runs(run).folder '\' runs(run).name '\func']);
%     sampled_img_list = dir([runs(run).folder '\' runs(run).name '\func\*Warped-to-Anat.nii.gz']);
%     N_vols = length(sampled_img_list);
% 
%     if run==1
%         layers_out = zeros([N_layers,N_ROI,N_vols,N_run]);
%     end
%     for vol = 1:N_vols
% %         sampled_img = load_nifti([sampled_img_list(vol).name]).vol;
%         sampled_img = niftiread([sampled_img_list(vol).name]);
% %         "E:\pain_layers\main_project\derivatives\pipeline\7349\ses-02\func\layers\rwls_stats_compcor_UNsmoothed_hpf180\spmT_0005.nii.gz"
% %         parfor ROI = 1:N_ROI
%         for ROI = 1:N_ROI
%             ind = find(ROIs(:,ROI));
%             tmp = zeros(size(ind,1),N_layers);
%             boundaries = layer_boundaries;
% %we go here from vertex to vertex(==column to column), and sample the funcitonal images. Here we should include an additional function, to select only columns active mostly in the deep layer. However, this need to be happen on the unsmoothed main effects (smoothed actually is available) conjuction map. That is calculate indiviudal conjuction maps, and apply the mask on those and estimate the layer profile. However, we would need to have a decent t value at the deeper part of the cortex. (my gut feeleing that interpolate statistical map is not that correct.-think about it a bit more)
%             for k = 1:size(ind,1) 
%                 X  = linspace(boundaries(ind(k),1,1),boundaries(ind(k),1,2),N_layers);
%                 Y  = linspace(boundaries(ind(k),2,1),boundaries(ind(k),2,2),N_layers);
%                 Z  = linspace(boundaries(ind(k),3,1),boundaries(ind(k),3,2),N_layers);
%                 
%                 tmp(k,:) = spm_sample_vol(sampled_img,X,Y,Z,SAMPLING_ORDER);
%             end
%             layers_out(:,ROI,vol,run) = squeeze(mean(tmp,1,'omitnan'));
%         end
%     end
%     fprintf(sprintf(' done \n'));
% end
% cd(currpath)
end