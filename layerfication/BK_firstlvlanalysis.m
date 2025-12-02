function [columnwisestat,layerwisestat,columndistribution] = BK_firstlvlanalysis(subid,subpath,fspath,T1path,visualizationtype,region,typeofcon,selectedvoxels)
% This funtion is called by the BK_layer_sampling_pain_study_pipeline
% function. The aim of this function is to do parameter estimation of columnar and layer level
% information per subject level. The input of this data is saved in the indiviudal files as
% interimdata_rwls_'ROINAME'(raw/derived).mat. This contains many information (see deatils in the corresponding fucntion) 
% regarding vertex wise (only taking the deeper 75% of column) timeseries, laminar timeseries from each vertex from the functinoal ROI, and the size of the
% mask (as the number of vertices). The script process ROI wise information for left and right side separately.
% The output is collection of structures, containing information about the
% raw ß parameters column and layer-wise. Also additional statistical
% information (pmax,T value,Tcrit). The information is provided per side
% separately. Additionally, the number of vertices active in the funcitonal
% mask is outputted. 
%
%
%INPUT:
% subid [str]                : subject id
% subpath [str]              : path to the derivatives folder
% fspath [str]               : path to the freesurfer folder
% T1path [str]               : full path and filename of the T1 image
% visualize [integer]        : 0-do not create, 1-create image of the
% selected vertices          
% region [string]            :the region of interest
% typeofcon [string]         : the type of contrast we base outr selection
% of voxels. Originally, I used the maineff (possible values:
% 'maineff','conditions','painlvls','cognlvls','main+int')
% selectedvoxels [string]    : whether if we slecet all voxels 'all', or only 200 (see topcolumnnumber)  
% e.g.:S1,S2,pIns,DLPFC,...
% typeofcon[str]             :the contrast to define functional ROIs. Originally, the
% group level mask of conjucntion of main effect of pain and cognition was
% used. Therefore, this is the default. But one can specify different type
% of conjucntion (painlvl,cognlvl- the conjucntion of the effect of
% pain/cognition in each level of cognition/pain). Also the contrasts of the laminea can be speccified: it can be the main effect, but each condition separatly is also possible.
%OUTPUT:
%   columnwisestat: [2x1 struct]: ß,T,T_crit,Pmax. it contains the fitted
% first lvl ß parameter estimates for all the veritces(e.g.: ß - # of columns X # of contrast) +optional statistical
% test result (e.g.: oncjunction analysis for selection of columns)
% 
%   layerwisestat: [2x1 cell]:  ß,T,T_crit,p_max. it contains the fitted
%   first lvl ß parameter estimated for all layers (sampled from the most
%   active columns, and averaged timeseries). e.g.: ß- # of layers X # of
%   contrasts
%
% Example call:
% [columnwisestat,columndistribution,layerwisestat]=BK_firstlvlanalysis(subid{:},subpath,fspath,T1path,0,region)
% 
% Balint Kincses
% balint.kincses@uk-essen.de
% 12.2024
    if ~ischar(subid)
        subid = char(string(subid));
    end
   
    N_layers = 20;
   %hardcode here the location of the imaging files.
    roipath='C:\Users\lenov\Documents\layerfMRI_DATA\groupavg_correctBET\';
    outputpath=[roipath subid '\functionalmasks\'];
    % load the data
    load([outputpath 'interimdata_rwls_' region 'raw.mat'],'interimdata_columns')
    %%this works with the averaged deep layer activation. it can be changed
    %%to prvoide information about the whole column (averaged within the
    %%column)
    %%original
    columnspecificts=interimdata_columns{1};
    %%modifed to work with the whole columns:
%     for ikj=1:2 %side
%         columnspecificts{ikj,1}=squeeze(mean(interimdata_columns{2}{ikj},2));
%     end
    sz=size(columnspecificts{1,1});
    if sz(2)<1221
        warning('The ts is SHORTER!!! why?')
    elseif sz(2)>1221
        warning('The ts is LONGER!!! why?')
    end
    
    layeractivation=interimdata_columns{2};
    columnsize=interimdata_columns{3};
    if strcmp(selectedvoxels,'all')
        for ROI=1:length(layeractivation)
            layerts=layeractivation{ROI,1};
            layerts_significant=mean(layerts(:,:,:),1); %,'omitnan'
            layerts_significant=squeeze(layerts_significant);
            [T,Tcrit,beta,pmax,percentsignalchange] = BK_layer_analysis_stats(layerts_significant,subid,subpath,region,'maineff+conditions'); %it always spits out the conditions in themself.
            columndistribution{ROI}=[];
            columnwisestat=[];
            layerwisestat{ROI} = struct('beta',beta,'T',T,'T_crit',Tcrit,'p_max',pmax);
            
        end
    elseif strcmp(selectedvoxels,'smoothed1thr')
        
    else
        %fit a 1st lvl GLM on each vertices.
        [T,Tcrit,beta,pmax] = BK_column_analysis_stats(columnspecificts,subid,subpath,region,typeofcon);
        if contains(typeofcon,'maineff') %maineff/maineff+conditions
            for side=1:length(T)
                T{side,1}(:,3)=min(T{side,1}(:,:),[],2); %conjuction analysis columnwise. In the third dimension, 1-cognition(lc-hc) OR (hc-lc), 2- pain (hp-lp)
            end
        elseif strcmp(typeofcon,'conditions')
            %TODO
        elseif strcmp(typeofcon,'painlvls')
            for side=1:length(T)
                T{side,1}(:,3)=min(T{side,1}(:,:),[],2); %conjuction analysis columnwise. In the third dimension, 1-effect of pain in the high cognition, 2- effect of pain on the low cognition
            end
        elseif strcmp(typeofcon,'cognlvls')
            for side=1:length(T)
                T{side,1}(:,3)=min(T{side,1}(:,:),[],2); %conjuction analysis columnwise. In the third dimension, 1-effect of pain in the high cognition, 2- effect of pain on the low cognition
            end
        elseif contains(typeofcon,'nociception+int') %maineff/maineff+conditions
            for side=1:length(T)
                T{side,1}(:,4)=min(T{side,1}(:,2:3),[],2); %conjuction analysis columnwise. In the third dimension,  2- pain (hp-lp) and interaction (positive like in beavhior),but it also includes the interaction effect tvalues

            end
        end
        
       
        columnwisestat = struct('beta',beta,'T',T,'T_crit',Tcrit,'p_max',pmax);
        %go through all the "ROI" (this would mena left and right in the
        %current implementation, and each anatomical ROI should be called
        %separately)
        for ROI=1:length(layeractivation)
            numberofcolumnsinmask=length(columnsize{ROI,1});
            tthrsss=[0 1];
            thrsidx=0;
            for tthrs=tthrsss
                thrsidx=thrsidx+1;
    %             numberofcognactivcolumn(thrsidx)=sum(columnwisestat.T(:,ROI,1)>tthrs);
                numberofcognactivcolumn(thrsidx)=sum(columnwisestat(ROI).T(:,1)>tthrs);
    %             numberofpainactivecolumns(thrsidx)=sum(columnwisestat.T(:,ROI,2)>tthrs);
                numberofpainactivecolumns(thrsidx)=sum(columnwisestat(ROI).T(:,2)>tthrs);
    %             numberofactcolumnsbasedonconj(thrsidx)=sum(columnwisestat.T(:,ROI,3)>tthrs);
                numberofactcolumnsbasedonconj(thrsidx)=sum(columnwisestat(ROI).T(:,3)>tthrs);
                activevoxelsid{ROI}=find((columnwisestat(ROI).T(:,3)>tthrs));
                disp(['Subject:' subid 'The number of columns within the funcitonal mask: ',num2str(numberofcolumnsinmask)])
                disp(['Subject:' subid 'The percent/number of columns acitvating in each run(t>', num2str(tthrs), ') for cognition: ',num2str(round(numberofcognactivcolumn(thrsidx)/numberofcolumnsinmask*100)),'%/', num2str(numberofcognactivcolumn(thrsidx))])
                disp(['Subject:' subid 'The percent/number of columns acitvating in each run(t>', num2str(tthrs), ') for pain: ',num2str(round(numberofpainactivecolumns(thrsidx)/numberofcolumnsinmask*100)),'%/',num2str(numberofpainactivecolumns(thrsidx))])
                disp(['Subject:' subid 'The percent/number of columns acitvating in each run(t>', num2str(tthrs), ') for cognition AND pain: ',num2str(round(numberofactcolumnsbasedonconj(thrsidx)/numberofcolumnsinmask*100)),'%/',num2str(numberofactcolumnsbasedonconj(thrsidx))])
                disp('-------------------------------------------------------------')
            end
            
            columndistribution{ROI} = struct('numberofcognactivcolumn',numberofcognactivcolumn, ...
                                        'numberofpainactivecolumns',numberofpainactivecolumns, ...
                                        'numberofactcolumnsbasedonconj',numberofactcolumnsbasedonconj, ...
                                        'numberofcolumnsinmask',numberofcolumnsinmask);
            %original, select the top 200 columns,so the most signficant voxel
            %pairs are included        
            topcolumnnumber=200; %200 or 100? I think the most important would be that, in most of the participant this 200 vertices are sampled (see below for checking if any of the conjunction is below 0, if so only the positives are estimated)
           %we can add the number here to be the total mask, so only the
           %positive values will be included
%            topcolumnnumber=numberofcolumnsinmask;
    
    
            %original, select top XX positive values
            if contains(typeofcon,'nociception+int') %maineff/maineff+conditions
                %we prescreen the vertices for those which do not have "strong"
                %interaction term.
%                 verticeswithnointeraction=columnwisestat(ROI).T(:,4)>1;
                %weonly focus on those vertices which response to
                %nocicpetion and also modified by attention.
                [top200Values, top200Indices] =maxk(columnwisestat(ROI).T(:,4),topcolumnnumber);
            else
                [top200Values, top200Indices] =maxk(columnwisestat(ROI).T(:,3),topcolumnnumber);
            end
            
            %Filter the indices where the values are positive
            positiveIndices = top200Indices(top200Values > 0);
            %Check if any negative values were removed
            if length(positiveIndices) < topcolumnnumber
                warning('The top 200 values contained negative numbers. Only %d positive values are retained in the %d ROI.', length(positiveIndices),ROI);
%                 disp(['The top 200 values contained negative numbers. Only ' num2str(length(positiveIndices)) ...
%                     ' positive values are retained in the ' num2str(ROI) ' ROI.']);
            end
    %%
            %select the top negative values
    %         [top200Values, top200Indices] =mink(columnwisestat(ROI).T(:,3),topcolumnnumber);
    %         positiveIndices = top200Indices(top200Values < 0);
    %         % Check if any negative values were removed
    %         if length(positiveIndices) < topcolumnnumber
    %             warning('The top 200 values contained positive numbers. Only %d negative values are retained.', length(positiveIndices));
    %         end
    %%
            %select all positive voxelpairs based on the t values from the
            %conjunction analysis
    %         positiveIndices=find(columnwisestat(ROI).T(:,3)>0);
            %attention values
    %         positiveIndices=find(columnwisestat(ROI).T(:,1)>0);
            % pain values
    %         positiveIndices=find(columnwisestat(ROI).T(:,2)>0);
    
    %%
            %select 200 random voxelpairs from around 0 t value of the
            %conjunction map
    %         nonsignmask=columnwisestat(ROI).T(:,3)<0.5 & columnwisestat(ROI).T(:,3)>-0.5;
    %         nonsignidx=find(nonsignmask);
    %         randomIndices = randperm(length(nonsignidx), topcolumnnumber);
    %         positiveIndices = nonsignidx(randomIndices);
            %cognition map
    %         nonsignmask=columnwisestat(ROI).T(:,1)<0.5 & columnwisestat(ROI).T(:,1)>-0.5;
    %         nonsignidx=find(nonsignmask);
    %         randomIndices = randperm(length(nonsignidx), topcolumnnumber);
    %         positiveIndices = nonsignidx(randomIndices);        
    
            layerts=layeractivation{ROI,1};
            activecluster{ROI,1}=positiveIndices; %this is the indices within the active cluster, so they do not represent an index in the whole surface, only in this subset
    %         sz=size(layerts);
            %calculate the mean timeseries from that 200 (or less) vertices.
    %         if any(isnan(layerts(positiveIndices,:,:)),'all')
    %             warining('The ts from the top vertices contain NA values')
    %         end
            %average the layer level information derived from the top 200
            %column:
            %original KEEP this:
            layerts_significant=mean(layerts(positiveIndices,:,:),1); %,'omitnan'
            layerts_significant=squeeze(layerts_significant);
    %%
            %average across all vertexpairs
    %         layerts_significant=mean(layerts(:,:,:),1); %,'omitnan'
    %%

    %%do some magic here and only inlcude a certain "bin" of this 200active
    %%columns
%     mybins=round(linspace(1,length(positiveIndices), 5));
%     for mybin=1:4
%         binstart=mybins(mybin);
%         binend=mybins(mybin+1);
%         layerts_significant=mean(layerts(positiveIndices(binstart:binend),:,:),1); %,'omitnan'
%         layerts_significant=squeeze(layerts_significant);
%         [T,Tcrit,beta,pmax] = BK_layer_analysis_stats(layerts_significant,subid,subpath,region,'maineff+conditions');
%             
%         layerwisestat{ROI,mybin} = struct('beta',beta,'T',T,'T_crit',Tcrit,'p_max',pmax);
%     end
    %         layerts_significant_forfunc=reshape(layerts_significant,[sz(2), sz(3),sz(4)]);
%             [T,Tcrit,beta,pmax] =BK_layer_analysis_stats(layerts_significant,subid,subpath,region,typeofcon);%this
%             was previously the call(until 10.08.2025),and it always
%             followed the contrast definition in the column selection.
%             However, it make more sense to calculate the effect for all
%             conditons separately, as a simple aritmetic can turn the
%             stuff into main effect,...

%this is the main stuf which we need to KEEP THIS
            [T,Tcrit,beta,pmax,percentsignalchange] = BK_layer_analysis_stats(layerts_significant,subid,subpath,region,'maineff+conditions');
            
            layerwisestat{ROI} = struct('beta',beta,'T',T,'T_crit',Tcrit,'p_max',pmax,'percentsignalchange',percentsignalchange);
        end
    end
    
    fprintf('time for laminar activation:')
    toc

    if strcmp(visualizationtype,'sampledROI')
        fprintf('Now create an image and save it...')
        layer_boundaries = VPF_load_layer_boundaries(subid,fspath);
        [layer_boundaries,~] = VPF_transform_layers_to_matlab_space(layer_boundaries,T1path);
        ROIs = BK_convert_load_ROIs(subid,roipath,fspath,size(layer_boundaries,1),region);
        BK_plotactiveclusterscolumn(T1path,ROIs,activecluster,layer_boundaries,outputpath,region)
    end
end

%% Column wise 1st lvl GLM estimate -  ROIs (sides) as different cells - checked, works just fine
function [T,T_crit,beta,p_max]=BK_column_analysis_stats(columnspecificts,subid,subpath,region,typeofcon)%ZTRANS)
% This is a function which loads the previosuly runned whole brain rwls
% estimation 
    %this is not implemented now, consider for visualizing raw ts data for
    %each trial
%     if nargin < 4
%         ZTRANS = false;
%     end

    N_ROIS=size(columnspecificts,1);    
    for roi=1:N_ROIS
        columnsss(roi)=size(columnspecificts{roi,1},1);
    end
%     N_columns=max(columnsss);
    T = cell(N_ROIS,1);
    beta = cell(N_ROIS,1);
    T_crit = cell(N_ROIS,1);
    p_max = cell(N_ROIS,1);

    statspath = dir([subpath num2str(subid) '/ses-02/func/layers/rwls_stats_compcor_UNsmoothed_hpf180']);
    load([statspath(1).folder '/SPM.mat'],'SPM');

    %define contrasts:
    % a reminder about the first lvl GLM
    %     {'anticipation_high_cognition'}
    %     {'anticipation_low_cognition' }
    %     {'pain_high_cogn_high_pain'   }
    %     {'pain_high_cogn_low_pain'    }
    %     {'pain_low_cogn_high_pain'    }
    %     {'pain_low_cogn_low_pain'     }
    %     {'rating'                     }
    if contains(typeofcon,'maineff') %maineff/maineff+conditions
        contrast = [4 8]; %main effects --> calculate later come conjunction stuff?, I would go now with the average effect, and not with the individual run effects.
        ncontrast = length(contrast);
        for idxcon=1:ncontrast
                if contrast(idxcon)==4 && ~strcmp(region,'DLPFC')
                    reversefact=-1; %this would mean that the activation is lower in the high cognition condition, that is higher in the low cognitive load. In other words, some supression occured in the high cognitive load condition. This should be reversed only for the praimary pain areas, but not the DLPFC.
                else
                    reversefact=1;
                end
                contastofinterest(:,idxcon)=SPM.xCon(contrast(idxcon)).c*reversefact;
        end
                
            
    elseif strcmp(typeofcon,'conditions')
        contrast = [3, 4,5 ,6]; %condition effects. to model every condition separately
        numberofregressors=18;
        contrasts=[contrast;contrast+numberofregressors;contrast+(2*numberofregressors)];
        ncontrast = length(contrast);
        nregr=length(SPM.xCon(1).c);
        contastofinterest=zeros(nregr,ncontrast);
        for idxcon=1:ncontrast
            contastofinterest(contrasts(:,idxcon),idxcon)=1;
        end
    elseif strcmp(typeofcon,'painlvls')
        %these contrast are for effect of pain in the high cognition condition and
        %effect of pain in the low cognition condition:
        contrast_pos = [3, 5 ];
        contrast_neg = [4, 6 ];
        numberofregressors=18;
        contrasts_pos=[contrast_pos;contrast_pos+numberofregressors;contrast_pos+(2*numberofregressors)];
        contrasts_neg=[contrast_neg;contrast_neg+numberofregressors;contrast_neg+(2*numberofregressors)];
        ncontrast = 2;
        nregr=length(SPM.xCon(1).c);
        contastofinterest=zeros(nregr,ncontrast);
        for idxcon=1:ncontrast
            contastofinterest(contrasts_pos(:,idxcon),idxcon)=1;
            contastofinterest(contrasts_neg(:,idxcon),idxcon)=-1;
        end
    elseif strcmp(typeofcon,'cognlvls')
        %these contrast are for effect of cognition in the high pain condition and
        %effect of cogn in the low pain condition:
        contrast_pos = [5 ,6];
        contrast_neg = [3, 4];
        numberofregressors=18;
        contrasts_pos=[contrast_pos;contrast_pos+numberofregressors;contrast_pos+(2*numberofregressors)];
        contrasts_neg=[contrast_neg;contrast_neg+numberofregressors;contrast_neg+(2*numberofregressors)];
        ncontrast = 2;
        nregr=length(SPM.xCon(1).c);
        contastofinterest=zeros(nregr,ncontrast);
        for idxcon=1:ncontrast
            contastofinterest(contrasts_pos(:,idxcon),idxcon)=1;
            contastofinterest(contrasts_neg(:,idxcon),idxcon)=-1;
        end
        if strcmp(region,'DLPFC')
            contastofinterest=contastofinterest*(-1);
        end
    elseif strcmp(typeofcon,'nociception+int') %maineff/maineff+conditions
        contrast = [4 8 12]; %main effects --> calculate later come conjunction stuff?, I would go now with the average effect, and not with the individual run effects.
        ncontrast = length(contrast);
        for idxcon=1:ncontrast
                if contrast(idxcon)==4 && ~strcmp(region,'DLPFC')
                    reversefact=-1; %this would mean that the activation is lower in the high cognition condition, that is higher in the low cognitive load. In other words, some supression occured in the high cognitive load condition. This should be reversed only for the praimary pain areas, but not the DLPFC.
                else
                    reversefact=1;
                end
                contastofinterest(:,idxcon)=SPM.xCon(contrast(idxcon)).c*reversefact;
        end
    end

    %load SPM matrix information for estimation
    W = SPM.xX.W;
    GLM = SPM.xX.pKX;

    for ROI = 1:N_ROIS
                
        Y = squeeze(columnspecificts{ROI,1}); %in the original pipeline, the size of Y is layernnum X tslength as here.
        if ~all(any(Y ~= 0, 1))
            BK_displaytxt("There is (at least) one column which has 0 value")
        end
        if any(isnan(Y),'all')
            BK_displaytxt("There is (at least) one column NA value")
        end
        Y = Y(:,any(Y ~= 0, 1)); %selects only the vertices/columns(rows in this matrix) where at least one element of the ts is non-zero, effectively removing all-zero columns from Y
%         if ZTRANS
%             %baseline z-transform. I take the volumes corresponding
%             % to 0 in the sum of all pain trials as baseline
%             idx = find(sum(SPM.xX.X(:,3:22),2)==0);
%             m = mean(Y(:,idx),2);
%             s = std(Y(:,idx),[],2);
%             Y = (Y-m)./s;
%         end

        KWY = spm_filter(SPM.xX.K,W*Y.');
        b   = GLM*KWY;

        res      = spm_sp('r',SPM.xX.xKXs,KWY);        %-Residuals
        ResSS    = sum(res.^2);                    %-Residual SSQ
        ResMS = ResSS / SPM.xX.trRV;

        for idxcontastofinterest = 1:ncontrast
            %We always assume a one-sided effect.            
%             contastofinterest=contrasts(:,idxcontastofinterest);
          
            con=contastofinterest(:,idxcontastofinterest);
            [T{ROI,1}(1:columnsss(ROI),idxcontastofinterest),...
             T_crit{ROI,1}(idxcontastofinterest),...
             beta{ROI,1}(1:columnsss(ROI),idxcontastofinterest),... %this should be a cell array instead of a matrix as the number of columns iwhtin a ROI most porbably change.
             p_max{ROI,1}(idxcontastofinterest),...
             ~] = BK_Tmap_from_SPM(SPM,b,ResMS,con,0.05,'none','column');
       end
    end

end

%% Layerwise 1st lvl GLM estimation - 
function [T,T_crit,beta,p_max,percentsignalchange]=BK_layer_analysis_stats(laminarts,subid,subpath,region,typeofcon)%subpath,ZTRANS)
%     if nargin < 4
%         ZTRANS = false;
%     end
    
    [N_layer,~] = size(laminarts);
    %     {'anticipation_high_cognition'}
    %     {'anticipation_low_cognition' }
    %     {'pain_high_cogn_high_pain'   }
    %     {'pain_high_cogn_low_pain'    }
    %     {'pain_low_cogn_high_pain'    }
    %     {'pain_low_cogn_low_pain'     }
    %     {'rating'                     }

    statspath = dir([subpath num2str(subid) '/ses-02/func/layers/rwls_stats_compcor_UNsmoothed_hpf180']);

    load([statspath(1).folder '/SPM.mat'],'SPM');


    if strcmp(typeofcon,'maineff')
        contrast = [4 8]; %main effects --> calculate later come conjunction stuff?, I would go now with the average effect, and not with the individual run effects.
        ncontrast = length(contrast);
        for idxcon=1:ncontrast
                if contrast(idxcon)==4 && ~strcmp(region,'DLPFC')
                    reversefact=-1; %this would mean that the activation is lower in the high cognition condition, that is higher in the low cognitive load. In other words, some supression occured in the high cognitive load condition. This should be reversed only for the praimary pain areas, but not the DLPFC.
                else
                    reversefact=1;
                end
                contastofinterest(:,idxcon)=SPM.xCon(contrast(idxcon)).c*reversefact;
        end                        
    elseif strcmp(typeofcon,'maineff+conditions')
        contrast = [3,4,5,6]; %condition effects. to model every condition separately
        numberofregressors=18;
        contrasts=[contrast;contrast+numberofregressors;contrast+(2*numberofregressors)];
        ncontrast = length(contrast);
        nregr=length(SPM.xCon(1).c);
        contastofinterest=zeros(nregr,ncontrast);
        for idxcon=1:ncontrast
            contastofinterest(contrasts(:,idxcon),idxcon)=1;
        end
    elseif strcmp(typeofcon,'main+int') %maineff/maineff+conditions
%         contrast = [4 8 12]; %main effects --> calculate later come conjunction stuff?, I would go now with the average effect, and not with the individual run effects.
%         ncontrast = length(contrast);
%         for idxcon=1:ncontrast
%                 if contrast(idxcon)==4 && ~strcmp(region,'DLPFC')
%                     reversefact=-1; %this would mean that the activation is lower in the high cognition condition, that is higher in the low cognitive load. In other words, some supression occured in the high cognitive load condition. This should be reversed only for the praimary pain areas, but not the DLPFC.
%                 else
%                     reversefact=1;
%                 end
%                 contastofinterest(:,idxcon)=SPM.xCon(contrast(idxcon)).c*reversefact;
%         end

        contrast = [3, 4,5 ,6]; %condition effects. to model every condition separately
        numberofregressors=18;
        contrasts=[contrast;contrast+numberofregressors;contrast+(2*numberofregressors)];
        ncontrast = length(contrast);
        nregr=length(SPM.xCon(1).c);
        contastofinterest=zeros(nregr,ncontrast);
        for idxcon=1:ncontrast
            contastofinterest(contrasts(:,idxcon),idxcon)=1;
        end
    end
%     N_contrasts=length(contrasts);

    T = zeros(N_layer,ncontrast); 
    beta = zeros(N_layer,ncontrast);
    T_crit = zeros(ncontrast,1);
    p_max = zeros(ncontrast,1);
    %added later
    percentsignalchange= zeros(N_layer,ncontrast);


    W = SPM.xX.W;
    GLM = SPM.xX.pKX;

    Y = squeeze(laminarts(:,:)); %in the original pipeline, the size of Y is layernnum X tslength.
    Y = Y(:,any(Y ~= 0, 1)); %selects only the columns where at least one element is non-zero, effectively removing all-zero columns from Y
%             if ZTRANS
%                 %baseline z-transform. I take the volumes corresponding
%                 % to 0 in the sum of all pain trials as baseline
%                 idx = find(sum(SPM.xX.X(:,3:22),2)==0);
%                 m = mean(Y(:,idx),2);
%                 s = std(Y(:,idx),[],2);
%                 Y = (Y-m)./s;
%             end
    %data prewhitening. align with spm, but we need to transpose the Y as
    %its columns represent time
    KWY = spm_filter(SPM.xX.K,W*Y.');
    %this is literally the same as in SPM
    b   = GLM*KWY;
    %SPM.xX.xKXs is the whitened and filtered design matrix (ensuring htat
    %residuals account for temporal autocorr)
    %literally the same as in SPM
    res      = spm_sp('r',SPM.xX.xKXs,KWY);        %-Residuals
    ResSS    = sum(res.^2);                    %-Residual SSQ
    ResMS = ResSS / SPM.xX.trRV;

%     for idxcontastofinterest = 1:numel(contrasts)
%         contastofinterest=contrasts(idxcontastofinterest);
%         if contrasts(idxcontastofinterest)==4 && ~strcmp(region,'DLPFC')
%             reversefact=-1; %this would mean that the activation is lower in the high cognition condition.
%         else
%             reversefact=1;
%         end
%         %     {'anticipation_high_cognition'}
%         %     {'anticipation_low_cognition' }
%         %     {'pain_high_cogn_high_pain'   }
%         %     {'pain_high_cogn_low_pain'    }
%         %     {'pain_low_cogn_high_pain'    }
%         %     {'pain_low_cogn_low_pain'     }
%         %     {'rating'                     }
%         con=SPM.xCon(contastofinterest).c*reversefact;
    for idxcontastofinterest = 1:ncontrast
        con=contastofinterest(:,idxcontastofinterest);
        
        [T(:,idxcontastofinterest),...
         T_crit(idxcontastofinterest),...
         beta(:,idxcontastofinterest),... %this should be a cell array instead of a matrix as the number of columns iwhtin a ROI most porbably change.
         p_max(idxcontastofinterest),...
         percentsignalchange(:,idxcontastofinterest)] = BK_Tmap_from_SPM(SPM,b,ResMS,con,0.05,'none','layer');
    end

end
%% calculate columnwise t-values - not checked,but trust in VPF that the implementation is correct.
function [T,Tcrit,con_array,pmax,percentsignal_array] = BK_Tmap_from_SPM(SPM,beta,ResMS,con,alpha,flag,typeofanalysis)

    if any(isnan(beta))
        warning('beta values contain NaNs! this is an issue and has to be futher investigated!')
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

    percentsignal_array     = zeros(1,size(beta,2));

    if strcmp(typeofanalysis,'column')
        disp('Column wise stats are calculated')
    elseif strcmp(typeofanalysis,'layer')
        runspecificconstant=beta(end-2:end,:);
        
        for j=1:size(beta_use,1)
            runspecificpercentsignal=(con(beta_index(j)) * beta_use(j,:))./runspecificconstant(j,:);
            percentsignal_array = percentsignal_array + runspecificpercentsignal;
        end
        percentsignal_array=percentsignal_array/3;
    end
%     con_array_spm=con'*beta; %might change to this as it is more stable
%     then the implementation above
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
%% functions for visualization
% load layers function 
%load layer boundaries information from .gii files
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

% transform between fs and matlab (voxel) spaces.
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

% Load the ROI surfaces and convert before if it does not exist (wb_command)
% convert the ROI to surface (wb_command) - main function (checked -keep it)
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
    

    
    localizer = dir(fullfile([roipath subid '\functionalmasks\'], '**',['roi_' region '*.nii']));
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

% Main visualization function - active columns in matlab space
function BK_plotactiveclusterscolumn(imgpath,ROIs,activevoxelsid,layer_boundaries,outputpath,region)
% function BK_plotROIverticesinmatlabspace(imgpath,ROIs,layer_boundaries,region,outputpath)
    img_T1 = spm_read_vols(spm_vol(imgpath));
    fig = figure;%('Visible', 'off');
        leftsidemask=find(ROIs(:,1));
%         roimaskbycolumn=find(ROIs(:,1));
        rightsidemask=find(ROIs(:,2));
%         leftsideclust=layer_boundaries(leftsidemask(activevoxelsid{1,1}),:,1);
        leftsideclust=leftsidemask(activevoxelsid{1,1});
        rightsideclust=rightsidemask(activevoxelsid{2,1});
        mostcolumns=mode(round(layer_boundaries([leftsideclust],3,1)));
%         leftsideclust=leftsidemask(104);
        if mostcolumns<5
            mostcolumns=5;
        end
        slice_idx=0;
%         zoom=30;
%force to go with a minimun value for the column
        sliceslab=2;
        for slice = mostcolumns-4:2:mostcolumns+6
            idx = find(layer_boundaries(:,3,1)>slice & layer_boundaries(:,3,1)<slice+sliceslab);
            idx3 = find(layer_boundaries(:,3,2)>slice & layer_boundaries(:,3,2)<slice+sliceslab);
            
%             idx_roimask = intersect(idx,leftsidemask);
            idx_roimask = intersect(idx,[leftsideclust; rightsideclust]);
%             idx3_roimask = intersect(idx3,leftsidemask);
            idx3_roimask = intersect(idx3,[leftsideclust; rightsideclust]);

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
        print(fig, [outputpath region '_ROIclust_forallposvertex.pdf'], '-dpdf', '-bestfit');
        savefig(fig, [outputpath region '_ROIclust_forallposvertex.fig']);  % Save as .fig
        close(fig);

end

%% this is some separate visualization thingy:
% sliceslab=2;
%         for slice = 78;%mostcolumns%-4:2:mostcolumns+6
%             idx = find(layer_boundaries(:,1,1)>slice & layer_boundaries(:,1,1)<slice+sliceslab);
%             idx3 = find(layer_boundaries(:,1,2)>slice & layer_boundaries(:,1,2)<slice+sliceslab);
%             
%             idx_roimask = intersect(idx,leftsidemask);
% %             idx_roimask = intersect(idx,[leftsideclust; rightsideclust]);
%             idx3_roimask = intersect(idx3,leftsidemask);
% %             idx3_roimask = intersect(idx3,[leftsideclust; rightsideclust]);
% 
% %             slice_idx = slice_idx + 1;
%     
%             % Create a subplot for each slice
% %             ax = subplot(2, 3, slice_idx); % 2 rows, ceil() for columns
% %             pos = get(ax, 'Position'); % Get current position
% %             pos(3) = pos(3) * 1.2;     % Increase width by 20%
% %             pos(4) = pos(4) * 1.2;     % Increase height by 20%
% %             set(ax, 'Position', pos);
%             
%             imagesc(squeeze(img_T1(slice,:,:))', [0 1200]), colormap(gray(256)), title(['Slice ' num2str(slice)])
% %             zoomedincoord=round([min(layer_boundaries(idx_roimask,1,1))-zoom  ...
% %                 max(layer_boundaries(idx_roimask,1,1))+zoom  ...
% %                 min(layer_boundaries(idx_roimask,2,1))-zoom ...
% %                 max(layer_boundaries(idx_roimask,2,1))+zoom]);
% %             axis(zoomedincoord);
%             
%             
%             hold on;
%             
%             plot(layer_boundaries(idx,2,1),layer_boundaries(idx,3,1),'r.');
%             plot(layer_boundaries(idx3,2,2),layer_boundaries(idx3,3,2),'y.');
%             
%             plot(layer_boundaries(idx_roimask,2,1),layer_boundaries(idx_roimask,3,1),'g.');
%             plot(layer_boundaries(idx3_roimask,2,2),layer_boundaries(idx3_roimask,3,2),'b.');
% %               
%             columntogether=[103846,104390]';%intersect(idx_roimask,idx3_roimask);
%             x1=layer_boundaries(columntogether,2,1);
%             x2=layer_boundaries(columntogether,2,2);
%             y1=layer_boundaries(columntogether,3,1);
%             y2=layer_boundaries(columntogether,3,2);
%             for togethervertices=1:length(columntogether)
%                 plot([x1(togethervertices), x2(togethervertices)], [y1(togethervertices), y2(togethervertices)], '-', 'LineWidth', 2, 'MarkerSize', 8,'Color','m');
%             end
% 
%             legend('white','pial');
%             axis off;
%         end
