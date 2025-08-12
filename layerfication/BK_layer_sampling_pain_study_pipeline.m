function [columnwisestat,layerwisestat,columndistribution]=BK_layer_sampling_pain_study_pipeline(subid,aim,region,samplingtype,visualizationtype,typeofcon,selectedvoxels)
% This is a modified script of VPF's same function. 
% I removed many part of it. Thos parts which has saved the output in a
% permanent location (in the external HDs), were removed. Additionally,
% output warning messages when these files are somehow not available. 
% 
% The function calls
% 
% Balint Kincses
% 08.2024
% 
% List of all subjects:
%7349, 7356,7361, 7375, 7376,7383, 7402, 7403, 7404,7405, 7408, 7414, 7415, 7425, 7426, 7433, 7434, 7435, 7443, 7444, 7445, 7448, 7449, 7452, 7453, 7454, 7455, 7456, 7457, 7468, 7469, 7482, 7484, 7485
% 
%
% S1 sampling was done for all subjects (interimdata_rwls.mat locally)
% S1derived 
% pIns sampling was done for all subjects (interimdata_rwls_pIns.mat locally)
% pIns derived
% S2 sampling was done for all (interimdata_rwls_S2.mat locally)
% S2 sampling for processed images (con and t) was done for all (interimdata_rwls_S2stat.mat locally)
% DLPFC 
% DLPFCderived is done for 7485
% 
% "Null" region - to test whether the seen layer profile is region
% specific: Aud:16 and 117 with no
% 
% 
% INPUTS:
%   subid [cell]:               the subjectID based on the ELHID
%   aim [char]:                 sample OR glmestimate - specify if layer level sampling OR the estimation of 1st lvl GLM should be done. The sampling must be done before estimation. Sampling saves an interim mat file, which is loaded by tge estimator.
%   region [char]:              specify the ROI. e.g.: S1,S2,pIns, DLPFC. The indiviudal ROI must be in the input folder of the subject as an .nii file (mind unzipping!)
%   samplingtype [char]:        raw/derived - specify the images which should be sampled. The raw ts data OR the already derived images ß/t. 
%   visualizationtype [char]:   no, onlyROI, ROI, sampledROI, . - specify whether a visualization of the surface points should be outputted. The onlyROI/ROI should be used with sampling and sampledROI should be used with glmestimation.
%                                   no - NO visualization is done
%                                   onlyROI - (should be used with aim==sample) no sampling is done, only the visualization of the ROI.
%                                   ROI - (should be used with aim==sample) whole ROI mask is visualized
%                                   sampledROI - (should be used with aim==glmestimate) is a subset of the indiviudal functional mask from which the sampling was done
% 
% OUTPUTS:
% Depends on the specified arguments:
% sample - interim data saved in the subject folder, containing timeseries
%   of the "columns" (incldued in the individual functional mask) deeper part, and layer level timeseries.
% glmestimate - layer level 1st estimation of contrasts. This could be later used as input in statistcal modelling and visualization. 
% visualization -  a pdf file which shows the vertices of interest
% 
% Example call of the function:
% to sample the raw fMRI data, and fit a 1st lvl GLM to the averaged ts
% signal (across column and across layers). The benefit of this approach is
% that the sampling happens before any modelling. 
% for subj={7485};BK_layer_sampling_pain_study_pipeline(subj,'sample','pIns','raw','no'); end
% 
% For sampling the derived image(t/ß) one can call the following ways
% for subj={7485};BK_layer_sampling_pain_study_pipeline(subj,'sample','S1','derived','no'); end
% 
% 
% 
% s
% diary subj8subj_layerif; for subj={7408,  7415,  7426,  7434,  7443,  7445,  7449,  7453}; BK_layer_sampling_pain_study_pipeline(subj,'glmestimate','S1','raw','no'); end; diary off
% 
% 
% 
% Balint Kincses
% balint.kincses@uk-essen.de
% 2025
% 

%hardcode freesurfer output path so the T1 images can be located
fspath = 'E:\pain_layers\main_project\derivatives\freesurfer\';

tic
    if subid{1} <= 7405
        subpath = 'E:/pain_layers/main_project/derivatives/pipeline/';
    else
        subpath = 'D:/main_project/derivatives/pipeline/';
    end
    if subid{1}==7356
        T1path = [subpath num2str(subid{1}) '' ...
            '\ses-01\anat\presurf_MPRAGEise\presurf_UNI\UNI_MPRAGEised_biascorrected.nii'];
%             '\ses-02\func\layers\run1\func\mag_POCS_r1_1000_Warped-to-Anat.nii.gz'];
                      
    else
        T1path = [subpath num2str(subid{1}) '' ...
              '\ses-01\anat\presurf_MPRAGEise\presurf_UNI\UNI_MoCo_MPRAGEised_biascorrected.nii'];
    end
        BK_checkfile(T1path)
    
        if strcmp(aim,'sample')
            if strcmp(samplingtype,'raw')
                BK_displaytxt('We sample the raw fmri image!')
            elseif strcmp(samplingtype,'derived')
                BK_displaytxt('We sample the t stat/ß images!')
            end
%             if strcmp(region,'Aud')
%                 BK_displaytxt('We only select a radnom 200 vertex pairs from a "null ROI"')
%                 BK_ROI_creation_GROUPLVLactivationconj_nullroi(subid{1},subpath,fspath,T1path,visualizationtype,region,samplingtype); % visualization!!!the 5th argument.
%             else
                BK_ROI_creation_GROUPLVLactivationconj(subid{1},subpath,fspath,T1path,visualizationtype,region,samplingtype); % visualization!!!the 5th argument.
%             end

        elseif strcmp(aim,'glmestimate')
            [columnwisestat,layerwisestat,columndistribution]=BK_firstlvlanalysis(subid{1},subpath,fspath,T1path,visualizationtype,region,typeofcon,selectedvoxels); %no visualization            
        else
            warning('Nothing was selected!')
        end
    fprintf('Time for one subject to process:')
    toc
end

