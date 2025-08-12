function [] = BK_ROI_creation(subid,pipe_path,fs_path)
% This is a modifed version of VPF's similar function. As the original
% function saved many interim files in the permanent locatoin (external
% HDs), I just check if those files exist or not.
% This function aims to handle the smoothed main data. The localitAlso, the 
% WM_path = [pipe_path '/' subid '/ses-02/func/WM_localizer/rwls_stats'];

%if the post calibration does not show anything reasonable, try using pain_calib
% pain_path = [pipe_path '/' subid '/ses-02/func/post_calib/rwls_stats'];
path_maintask= [pipe_path '/' subid '/ses-02/func/layers/rwls_stats_compcor_smoothed/'];
cognition='spmT_0004.nii.gz';
path_WM=[path_maintask cognition];
pain='spmT_0008.nii.gz';
path_pain=[path_maintask pain];
subj_fs_path_mri = [fs_path '/' subid '/mri'];

%used for working memory localizer
% aparc     = [subj_fs_path_mri '/aparc+aseg.mgz'];
% [~,~,aparc_ext] = fileparts(aparc);
% 
% %used for pain localizer
% aparc2009 =  [subj_fs_path_mri '/aparc.a2009s+aseg.mgz'];

% thresholding (t>1)
t_thres = 1;
%% transform atlases into native subject space and from .mgz to .nii.gz
% fs_base = ['export FREESURFER_HOME=/usr/local/freesurfer/6.0.0; ' ...
%     'export SUBJECTS_DIR=' subj_fs_path_mri '; '...
%     'source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];

%aparc
%aprac.mgz is created
%not sure if it is an issue whether in the mri_label2vol funciton the same
%file (aseg.mgz) was added to regheader for converting aparc and aparc2009 to .mgz
%extension. But I would say that is OK.
% cmd = ['mri_label2vol --seg ' '"' aparc '"' ' --temp ' '"' subj_fs_path_mri '/rawavg.mgz" --o '...
%     '"' subj_fs_path_mri '/aparc' aparc_ext '"' ' --regheader ' '"' subj_fs_path_mri '/aseg.mgz"'];
% system([fs_base cmd]);
% path_aseg=[subj_fs_path_mri '/aseg.mgz'];
% BK_checkfile(path_aseg)
% cmd = ['mri_convert ' '"' subj_fs_path_mri '/aparc.mgz" ' '"' subj_fs_path_mri '/aparc.nii.gz"'];
% system([fs_base cmd]);
% aparc = [subj_fs_path_mri  '/aparc.nii.gz'];
% aparc = [subj_fs_path_mri  '/aparc.nii.gz'];
path_aparc = [subj_fs_path_mri  '/aparc.nii.gz'];
BK_checkfile(path_aparc)

%aparc2009
%aparc2009.mgz is created
% cmd = ['mri_label2vol --seg ' '"' aparc2009 '"' ' --temp ' '"' subj_fs_path_mri '/rawavg.mgz" --o '...
%     '"' subj_fs_path_mri '/aparc2009' aparc_ext '"' ' --regheader ' '"' subj_fs_path_mri '/aseg.mgz"'];
% system([fs_base cmd]);

% cmd = ['mri_convert ' '"' subj_fs_path_mri '/aparc2009.mgz" ' '"' subj_fs_path_mri '/aparc2009.nii.gz"'];
% system([fs_base cmd]);
% 
% aparc2009 = [subj_fs_path_mri '/aparc2009.nii.gz'];
path_aparc2009 = [subj_fs_path_mri '/aparc2009.nii.gz'];
BK_checkfile(path_aparc2009)
%% working memory localizer
%load smoothed t-map
% if ~isfile([WM_path '/sspmT_0001.nii'])
%     spm_smooth([WM_path '/spmT_0001.nii'],[WM_path '/sspmT_0001.nii'],[4.5 4.5 4.5]);
% end
% WM_localizer_t_hdr = load_nifti([WM_path '/sspmT_0001.nii']);
% WM_localizer_t = WM_localizer_t_hdr.vol;

BK_checkfile(path_WM)
WM_localizer_t = niftiread(path_WM);
WM_localizer_t_hdr=niftiinfo(path_WM);

% intersect with anatomical masks coming from Freesurfer's aparc atlas (fist number = left hemisphere, 2nd = right)
ROI_idx = {1003,2003;1008,2008;1024,2024;1027,2027;1028,2028;1028,2028;1029,2029};
% atlas = load_nifti(path_aparc).vol;
atlas = niftiread(path_aparc);
ROIs = size(ROI_idx,1);

for ROI = 1:ROIs
    tmp = zeros(size(WM_localizer_t));
    tmp(atlas==ROI_idx{ROI,1}) = WM_localizer_t(atlas==ROI_idx{ROI,1});
    tmp(atlas==ROI_idx{ROI,2}) = WM_localizer_t(atlas==ROI_idx{ROI,2});

    for ii = [1,-1] %threshold positive and negative t-values seperately
        mask = zeros(size(WM_localizer_t),'logical');
        if ii == 1
            mask(tmp>t_thres) = true;
        else
            mask(tmp<-t_thres) = true;
        end

        WM_localizer_t_forsave = mask;
        tmp_name = num2str(ROI_idx{ROI,1});
        if ii == 1
            outname = [path_maintask '/WM_localizer_ROI' tmp_name(2:end) '_pos.nii'];
        else
            outname = [path_maintask '/WM_localizer_ROI' tmp_name(2:end) '_neg.nii'];
        end
        WM_localizer_t_hdr.Filename=outname;
        niftiwrite(single(WM_localizer_t_forsave),outname,WM_localizer_t_hdr,'Compressed',1);
%         gzip(outname);
%         delete(outname);
    end
end


%% pain localizer
%load smoothed t-map
% if ~isfile([pain_path '/sspmT_0001.nii'])
%     spm_smooth([pain_path '/spmT_0001.nii'],[pain_path '/sspmT_0001.nii'],[4.5 4.5 4.5]);
% end
% pain_localizer_t_hdr = load_nifti(path_pain);
% pain_localizer_t = pain_localizer_t_hdr.vol;

BK_checkfile(path_pain);
pain_localizer_t = niftiread(path_pain);
pain_localizer_t_hdr = niftiinfo(path_pain);


% intersect with anatomical masks coming from Freesurfer's aparc2009 atlas
%the last 4 ROIs will be merged into 2 ROIs
ROI_idx = {11104,12104;11106,12106;11107,12107;11128,12128;11168,12168;...
    11117,12117;11149,12149;...
    11118,12118;11150,12150};
% atlas = load_nifti(path_aparc2009).vol;
atlas = niftiread(path_aparc2009);

ROIs = size(ROI_idx,1);
for ROI = 1:ROIs-4
    tmp = zeros(size(pain_localizer_t));
    tmp(atlas==ROI_idx{ROI,1}) = pain_localizer_t(atlas==ROI_idx{ROI,1});
    tmp(atlas==ROI_idx{ROI,2}) = pain_localizer_t(atlas==ROI_idx{ROI,2});

    for ii = [1,-1]%threshold positive and negative t-values seperately
        mask = zeros(size(pain_localizer_t),'logical');
        if ii == 1
            mask(tmp>t_thres) = true;
        else
            mask(tmp<-t_thres) = true;
        end

        pain_localizer_t_forsave = mask;
        tmp_name = num2str(ROI_idx{ROI,1});
        if ii == 1
            outname = [path_maintask '/pain_localizer_ROI' tmp_name(2:end) '_pos.nii'];
        else
            outname = [path_maintask '/pain_localizer_ROI' tmp_name(2:end) '_neg.nii'];
        end
        pain_localizer_t_hdr.Filename=outname;
        niftiwrite(single(pain_localizer_t_forsave),outname,pain_localizer_t_hdr,'Compressed',1);
%         gzip(outname);
%         delete(outname);
    end
end

for ROI = ROIs-3:2:ROIs
    tmp = zeros(size(pain_localizer_t));
    tmp(atlas==ROI_idx{ROI,1}) = pain_localizer_t(atlas==ROI_idx{ROI,1});
    tmp(atlas==ROI_idx{ROI,2}) = pain_localizer_t(atlas==ROI_idx{ROI,2});

    tmp(atlas==ROI_idx{ROI+1,1}) = pain_localizer_t(atlas==ROI_idx{ROI+1,1});
    tmp(atlas==ROI_idx{ROI+1,2}) = pain_localizer_t(atlas==ROI_idx{ROI+1,2});

    for ii = [1,-1]%threshold positive and negative t-values seperately
        mask = zeros(size(pain_localizer_t),'logical');
        if ii == 1
            mask(tmp>t_thres) = true;
        else
            mask(tmp<-t_thres) = true;
        end

        pain_localizer_t_forsave = mask;
        tmp_name = num2str(ROI_idx{ROI,1});
        if ii == 1
            outname = [path_maintask '/pain_localizer_ROI' tmp_name(2:end) '_pos.nii'];
        else
            outname = [path_maintask '/pain_localizer_ROI' tmp_name(2:end) '_neg.nii'];
        end
        pain_localizer_t_hdr.Filename=outname;
        niftiwrite(single(pain_localizer_t_forsave),outname,pain_localizer_t_hdr,'Compressed',1);
%         gzip(outname);
%         delete(outname);
    end
end