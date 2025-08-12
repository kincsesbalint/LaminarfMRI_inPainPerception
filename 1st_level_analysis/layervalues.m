function [allsubjlayer,allprocsubjid,allsubjcoldistr] = layervalues(region,typeofcon,varargin)
% Run the first lvl GLM on each subject data in a specific ROI and save the ß/contrast estimates for each layer.
%
% ToDo: Details of the script
% example call:
% [S1_nociceptivecolumns,procesubjid,columndistr]=layervalues('S1','painlvls','visualization','no');
% by default the function use all subjects' data, do not do visualization
% of the indiviudal funcitonal ROI, and select only the top 200 vertex
% pairs("columns")
% sampledROI
% Balint Kincses
% balint.kincses@uk-essen.de
% 2025

    p = inputParser;

    % Required arguments
    addRequired(p, 'region', @ischar);
    addRequired(p, 'typeofcon', @ischar);

    % Optional arguments with default values
    %these are all the 34 subjects:
    allsubjid={7349, 7356, 7361, 7375, 7376,...
         7383, 7402, 7403, 7404,...
         7405, 7408, 7414, 7415,...
         7425, 7426, 7433, 7434,...
         7435, 7443, 7444, 7445,...
         7448, 7449, 7452, 7453,...
         7454, 7455, 7456, 7457,...
         7468, 7469, 7482, 7484, 7485};
    defpathtosubjdata='C:\Users\lenov\Documents\layerfMRI_DATA\groupavg_correctBET\';
    selectedcolumnvalue='200';
    addOptional(p, 'subjectids', allsubjid, @iscell);
    addOptional(p, 'visualization', 'no', @ischar);
    addOptional(p, 'pathtosubjdata', defpathtosubjdata, @ischar);
    addOptional(p, 'SELECTEDCOLUMNS', selectedcolumnvalue, @ischar);
    

    % Parse inputs
    parse(p, region,typeofcon, varargin{:});
    
    % Extract values
    subjid = p.Results.subjectids;
    visualize = p.Results.visualization;
    pathtosubj = p.Results.pathtosubjdata;
    selectedcolumns = p.Results.SELECTEDCOLUMNS;

    % Display values
    fprintf('We process in region: %s, number of subjects: %d \n', region,length(subjid));
    if length(subjid)<7
        fprintf('the following subjects: %d\n', subjid{:});
    else
        warning('The number of subjects >7, only the first seven is written out')
        fprintf('the following subjects: %d\n', subjid{1:7});
    end

    plotidx=0;
    allsubjlayer=cell(length(subjid),1);
    allsubjcoldistr=cell(length(subjid),1);
    allsubjstat=cell(length(subjid),1);
    for subj=subjid       
        columnartsfile=[pathtosubj num2str(subj{:}) '\functionalmasks\interimdata_rwls_' region 'raw.mat'];
        if exist(columnartsfile,"file")

            disp([num2str(subj{:}) ' subject columnar info r'])
            plotidx=plotidx+1;
            [stat,layerwisestat,columndistribution]=BK_layer_sampling_pain_study_pipeline(subj,'glmestimate',region,'raw',visualize,typeofcon,selectedcolumns);
            
            allsubjlayer{plotidx}=layerwisestat;
            allsubjcoldistr{plotidx}=columndistribution;
            allsubjstat{plotidx}=stat;
            allprocsubjid(plotidx)=subj{:};
        else
            fprintf('%d subject has missing %s file\n',subj{:},columnartsfile)
        end
%     
    end

end