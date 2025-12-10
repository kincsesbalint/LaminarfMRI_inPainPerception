%create some cluster statistics just for understanding the surfaces waht we
%work with:
% allcluster - is the output of the BK_smoothedsurfaceROI_creation with the
% following call:
allsubj={7349, 7356, 7361, 7375, 7376,...
         7383, 7402, 7403, 7404,...
         7405, 7408, 7414, 7415,...
         7425, 7426, 7433, 7434,...
         7435, 7443, 7444, 7445,...
         7448, 7449, 7452, 7453,...
         7454, 7455, 7456, 7457,...
         7468, 7469, 7482, 7484, 7485};
% maki=1;for subj=allsubj; [allcluster{maki},laminarresults{maki}]=BK_smoothedsurfaceROI_creation(subj,{'S2'},'grouplvlmask','expansion',false);maki=maki+1; end
% allcluster=cell(34,1);maki=1;for subj=allsubj; allcluster{maki}=BK_smoothedsurfaceROI_creation(subj,{'S2'},'thrs1','expansion',false);maki=maki+1; end
% allcluster=cell(34,1);maki=1;for subj=allsubj; allcluster{maki}=BK_smoothedsurfaceROI_creation(subj,{'S2'},'thrs0','expansion',false);maki=maki+1; end
% allcluster=cell(34,1);laminarresults=cell(34,1);maki=1;for subj=allsubj; allcluster{maki}=BK_smoothedsurfaceROI_creation(subj,{'S2'},'top200','expansion',false);maki=maki+1; end
% 
% Additional idea is to expand the mask! now that is the default call of
% the function above. The clusters are very small for the top200 and thrs1
% option, that means they need to be expanded additionally. The idea here
% is that if we capture only veins, which reflect also something related to
% neuronal activity, we can assume that the neuronal activity is close in
% proximity, that is we include the neighbours of those vein containing
% vertices.
% allcluster=cell(34,1);maki=1;for subj=allsubj; allcluster{maki}=BK_smoothedsurfaceROI_creation(subj,{'S2'},'top200');maki=maki+1; end
% allcluster=cell(34,1);maki=1;for subj=allsubj; allcluster{maki}=BK_smoothedsurfaceROI_creation(subj,{'S2'},'thrs1');maki=maki+1; end
%
% The above calls ended up in too big masks per individual. Most probably
% we cover too much then. Workaround:
% Focus at least clusters with two vertices, and include those neighbours:
% 
% 
% allcluster=cell(34,1);maki=1;for subj=allsubj; allcluster{maki}=BK_smoothedsurfaceROI_creation(subj,{'S2'},'thrs1','expansion',false,'twoverticesexpansion',true);maki=maki+1; end
% 
% implemented also the layer estimation from the selected vertices:
% 
% allcluster=cell(34,1);laminarresults=cell(34,1);maki=1;for subj=allsubj; [allcluster{maki},laminarresults{maki}]=BK_smoothedsurfaceROI_creation(subj,{'S2'},'thrs1','expansion',false,'twoverticesexpansion',true,'laminarest',true);maki=maki+1; end
% 
% 
% s
nameofregion='pIns';
pathtofigures=['C:\Users\lenov\Documents\layerfMRI\E_stats\images\taskfmri\surfaceclusterinformation\' nameofregion '_thrs1_twoverticesexpansion\'];

if ~isfolder(pathtofigures) % use exist(folderPath,'dir') for older MATLAB versions 
    mkdir(pathtofigures);
    fprintf('Created folder: %s\n', pathtofigures); 
else 
    warning('Folder "%s" already exists. Stopping script.', pathtofigures);
    return; % stops the current script/function after issuing the warning 
end
% Notes:
% If you need a hard stop regardless of context, replace the warning/return with: error('Folder "%s" already exists. Stopping.', folderPath);

thresholdvalue=[0.8,1,200,20];
for myside={'left','right'}
    interestingdata=struct([]);
    myidx=1;
    visualizedsubjdata=[];
    for i=1:34
        lens = cellfun(@numel, allcluster{i,1}.(myside{1}));
        [maxLen, idx] = max(lens);
        interestingdata(i).mainclusterratio=maxLen/sum(lens);
        interestingdata(i).numberofclusters=length(lens);
        
        interestingdata(i).totalmasksize=sum(lens);
        interestingdata(i).biggestclustersize=maxLen;
        if ~isempty(lens)
            visualizedsubjdata(myidx)=allsubj{i};
            myidx=myidx+1;
        end

    end
    whatweareinterestedin=fieldnames(interestingdata);
    for maki=1:length(whatweareinterestedin)
        if strcmp(whatweareinterestedin{maki},'numberofclusters') || strcmp(whatweareinterestedin{maki},'totalmasksize')
            tofigureindices=allsubj;
        else
            tofigureindices=visualizedsubjdata;
        end
        fig=figure;
        bar(vertcat(interestingdata(:).(whatweareinterestedin{maki})))
        hold on
        yline(thresholdvalue(maki),'r')
        ax = gca;
        ax.XTick = 1:length(tofigureindices);
        ax.XTickLabel = tofigureindices;
        ax.XTickLabelRotation = 90;
        title([myside{1} ' ' nameofregion ' ' whatweareinterestedin{maki} ' for each indiviudal'])
        myfilename=[pathtofigures whatweareinterestedin{maki} nameofregion myside{1} '.png'];
        saveas(fig,myfilename)
        close(fig)
    end
end
%leftside

% figure
% bar(vertcat(interestingdata(:).mainclusterratio))
% hold on
% yline(0.8,'r')
% title('Left S2 main cluster ratio in the grouplevel maskfor each indiviudal')
% 
% 
% interestingdata=struct([]);
% %rightside
% for i=1:34
%     lens = cellfun(@numel, allcluster{1,i}.right);
%     interestingdata(i).mainclusterratio=max(lens)/sum(lens);
%     interestingdata(i).numberofclusters=length(lens);
%     [maxLen, idx] = max(lens);
%     interestingdata(i).totalmasksize=sum(lens);
% end
% figure
% bar(vertcat(interestingdata(:).mainclusterratio))
% hold on
% yline(0.8,'r')
% title('Right S2 main cluster ratio in the grouplevel maskfor each indiviudal')

%it is still decent number of vertices belong hte biggest cluster(around
%300??). Anyway, this is a bit weird that they are different in the two
%sides...