function [clusterindices,layerwisestat]=BK_smoothedsurfaceROI_creation(subid,regions,masktype,varargin)
% This function is called by its own. This is a modification of BK_ROI_creation_GROUPLVLactivationconj function.
%
% This function aims to create continoues surfaces, which form a cluster.
% The initial data is the surface information from the individual, and we
% select only those vertices which form a continous surfcae. That means, we
% based our selection on vertex level stats, and we apply not a simple
% threshold (as previously to pick the top vertices) but also some cluster
% forming condition has to be met.
% 
% 
%  E
% a
% TODO the rest from here 
% The following steps are implemented here:
%   1. load boundaries to matlab which coming from freesurfer to the matrix layer_boundaries [vertex,[x,y,z],[white,pial]] : 3D matrix of the surfaces(VPF_load_layer_boundaries)
%   2. load ROIs coming from the group average activation (BK_convert_load_ROIs)
%   3. transform boundaries to matlab space (from freesurfer space), the
%        resulting space will be in "voxel coordinate".
% 
% 
% .
% 
% 
% EXAMPLE CALL:
% see
% C:\Users\lenov\Documents\GitHub\LaminarfMRI_inPainPerception\supplementary\clusterstat.m
% for examples how to call.
%
% Balint Kincses
% 12.2025
    p = inputParser;

    % Required arguments
    addRequired(p, 'subid', @iscell);
    addRequired(p, 'regions', @iscell); 
    addRequired(p, 'masktype', @ischar);
    
    expansiondef=true;
    addOptional(p, 'expansion', expansiondef, @islogical);
    twovertices=false;
    addOptional(p, 'twoverticesexpansion', twovertices, @islogical);
    laminarestimationdef=false;
    addOptional(p, 'laminarest', laminarestimationdef, @islogical);
    % Parse inputs
    parse(p, subid,regions,masktype, varargin{:});
    
    % Extract values
    expansion = p.Results.expansion;
    twoverticesexpansion = p.Results.twoverticesexpansion;
    laminarestimation = p.Results.laminarest;

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
fprintf('the subject IDwhich we work with:%s\n',subid)
%local folder where individual subject folders are located
roipath='C:\Users\lenov\Documents\layerfMRI_DATA\groupavg_correctBET\';

% the output of this function can be saved.
outputpath=[roipath subid '\smootheddata\'];
interimdatalocation=[roipath subid '\functionalmasks\'];

% smoothedmaskfilenm=[outputpath 'interimdata_smoothedmask_' region '.mat'];
% smoothedmaskfilenm=[outputpath 'interimdata_smoothedmask.mat'];%this is with 0 threshold 
smoothedmaskfilenm=[outputpath 'interimdata_smoothedmask_thr1.mat'];%this is with 1 threshold 

% modified helper function to load faces and vertex coordinates in a side
% specific manner(in the original it was concatenated  along the first
% dimension, which was OK for that purpose, but here the exact number and
% the indices are needed for defining the clusteres)
[layer_boundaries,faces,vertexcoordinates] = BK_load_layer_boundaries_sidespecific(subid,fspath);
fprintf('The number of vertices in the left hemisphere: %i\n',size(vertexcoordinates.left,1))
% Nv = size(layer_boundaries,1);
%the indices start with 0, which cannot be handled by MATLAB so handle the situation:
if min(faces.rh_w(:)) == 0 
    faces.rh_w = faces.rh_w + 1;
end
if min(faces.lh_w(:)) == 0 
    faces.lh_w = faces.lh_w + 1;
end

clusterindices=struct('left', [], 'right', []);
%go through all regions
for myregion=regions
    region=myregion{1};
    display(region)
    vertexwiseparamfilenm=[outputpath 'vertexwiseparamestimation_' region '.mat'];
    %create OR load the columnwosestat from the originally sampled mask
    if ~exist(vertexwiseparamfilenm ,'file')
        BK_displaytxt('The vertexwiseparamestimate file does NOT exist, so we create one and save it!')
        %hardcoded arguments of this function
        [columnwisestat,~,~]=BK_firstlvlanalysis(subid,subpath,fspath,T1path,'no',region,'maineff+int','200'); %no visualization            
        save(vertexwiseparamfilenm,'columnwisestat','-v7.3');
    else
        load(vertexwiseparamfilenm,'columnwisestat')
    end 
    
    
    maki=10;
    %load ROIs 
    localizer_original = dir(fullfile([roipath subid '\functionalmasks\'], '**',['roi_' region '*.nii.gz']));
        %load the group smoothed mask images' surface reconstruction (this
        %is the original, we sampled all the vertex pairs in this mask and
        %can use to individually select vertex pairs). The statistics of
        %each vertex pairs/columns can be found in the columnwisestat(each
        %vertex pair has a value of T which can be used for selection)

    ROIs_original = BK_convert_load_ROIs(subid,roipath,fspath,size(layer_boundaries,1),region,localizer_original); %the definition of hte function is differ a bit from the one which is in the BK_firstlvlanalysis as we use the zipped mask files and not the unzipped ones.
    for myside={'left','right'}
        if strcmp(myside,'left')
           vertexIDs=find(ROIs_original{1,1}.lh); %left side
           facestolookinto='lh_w';
           defineside=1;
        elseif strcmp(myside,'right')
           vertexIDs=find(ROIs_original{2,1}.rh); %right side
           facestolookinto='rh_w';
           defineside=2;
        end
        
        if strcmp(masktype,'grouplvlmask') % check the continouty of the indiviudally back- and surface-projected masks:
            vrtindicesforclustertesting=vertexIDs;
        elseif strcmp(masktype,'top200') % check the continouty of the top200 vertices. 
            [hami,laci]=maxk(columnwisestat(defineside).T(:,3),200);            
            vrtindicesforclustertesting=vertexIDs(laci);
        elseif strcmp(masktype,'thrs1') %check the continouty of the vertices with a value larger then 1
            laci=columnwisestat(defineside).T(:,3)>1;
            vrtindicesforclustertesting=vertexIDs(laci);
        elseif strcmp(masktype,'thrs0') %check the continouty of the vertices with a value larger then 1
            laci=columnwisestat(defineside).T(:,3)>0;
            vrtindicesforclustertesting=vertexIDs(laci);
        end
        if expansion
            vrtindicesforclustertesting=neigbourexpansion(faces.(facestolookinto),vrtindicesforclustertesting);
        end
        maki=clusteridentification(faces.(facestolookinto),vrtindicesforclustertesting);
%         [clustersMarkedOnly, clustersWithBridges] =clustersFromTripScanning_TwoHop(faces.(facestolookinto),vrtindicesforclustertesting);
        if twoverticesexpansion && ~expansion
            lens = cellfun(@numel, maki);
            hopi=maki(~(lens<2));
            twovertexclusters=vertcat(hopi{:});
            vrtindicesforclustertesting=neigbourexpansion(faces.(facestolookinto),twovertexclusters);
            maki=clusteridentification(faces.(facestolookinto),vrtindicesforclustertesting);
        end
        
        
        clusterindices.(myside{1})=maki;   
        maskindices.(myside{1})=ismember(vertexIDs,vertcat(maki{:}));
% find(hammmi);

    end

    if laminarestimation
        load([interimdatalocation 'interimdata_rwls_' region 'raw.mat'],'interimdata_columns')
        layeractivation=interimdata_columns{2};
        
%         positiveIndices=%this has to be something which codes within the groupmask(so it uses the group mask indices and not the whole mesh)
        for myside={'left','right'}
            if strcmp(myside,'left')
               ROI=1;
            elseif strcmp(myside,'right')
               ROI=2;
            end
            layerts=layeractivation{ROI,1};
            vertexindices=maskindices.(myside{1});
            if all(vertexindices==0) %if no clusters are found. This can be happen if we go with strict threhold(t>1) and with at least 2 neighbouring vertices)
                layerwisestat.(myside{1}) = struct('beta',[],'T',[],'T_crit',[],'p_max',[],'percentsignalchange',[]);
            else
                layerts_significant=mean(layerts(vertexindices,:,:),1); %,'omitnan'
                layerts_significant=squeeze(layerts_significant);
                [T,Tcrit,beta,pmax,percentsignalchange] = BK_layer_analysis_stats(layerts_significant,subid,subpath,region,'maineff+conditions');
                layerwisestat.(myside{1}) = struct('beta',beta,'T',T,'T_crit',Tcrit,'p_max',pmax,'percentsignalchange',percentsignalchange);
            end
        end
    else
       layerwisestat=struct(); 

    end
    
    
 
    
% %         ROIs_original = BK_convert_load_ROIs_originalmask(subid,roipath,fspath,size(layer_boundaries,1),region,localizer_original);
%         %select the top 200 voxels (irrespective of the exact value, only
%         %focus on the most active (it allows that all the values e.g. below
%         %t=1, which is a bit troublesome)
%         defineside=1; %1 is for left and 2 is for right
%         [hami,laci]=maxk(columnwisestat(defineside).T(:,3),200);
%         %only keep values which are above t=1, this is a bit stricter and
%         %makes more sense. However, this ends up only a few vertex per
%         %subject.
%         
%         locationoffullmask=find(ROIs_original{defineside,1}.lh);
%         maskof200=locationoffullmask(laci);
%         clusterindices=clusteridentification(faces.lh_w,maskof200);
%         [clustersMarkedOnly, clustersWithBridges] =clustersFromTripScanning_TwoHop(faces.lh_w,maskof200);


        

        
        % Transform boundaries to matlab space (from freesurfer space), the resulting space will be in "voxel coordinate".
        %this is crucial step and need to be checked thouroughly. the
        %sampling depends on this:!!
        %%for visualization:
%         [layer_boundaries,T1_mat,old_layer_boudnaries] = VPF_transform_layers_to_matlab_space(layer_boundaries,T1path);
%         BK_displaytxt('Visualization of the ROI + sampling')
%         BK_plotROIverticesinmatlabspace(T1path,ROIs_smoothed,layer_boundaries,region,outputpath)
        %this would code the original ROIs from which we made the sampling
        %back then:
%         for ROI=1:2
%             idxtostart=strfind(localizer_smoothed(ROI).name,region)+length(region);
%             myroiname=localizer_smoothed(ROI).name(idxtostart:idxtostart+3);
% %             ind = find(ROIs_smoothed(:,ROI));
%             intersectedwithgroupavgmask.(region).(myroiname)=ROIs_smoothed(ROIs_original(:,ROI),ROI);
%             %the original mask size can be calculate with length(intersectedwithgroupavgmask.(region).(myroiname))
%             originalmasksize=length(intersectedwithgroupavgmask.(region).(myroiname));
%             %the new mask size can be calcaulted with sum(ROIs_smoothed(:,ROI))
%             newmasksize=sum(ROIs_smoothed(:,ROI));
%             %the intersection can be calculated with sum(intersectedwithgroupavgmask.(myroiname))
%             intersectionofthetwomask=sum(intersectedwithgroupavgmask.(region).(myroiname));
%             newmasksizes.(region).(myroiname)=newmasksize;
%             fprintf('in %s ROI:\n',myroiname)
%             fprintf('The size of the original mask (group level average):%i\n', [originalmasksize])
%             fprintf('The size of the new mask (indiviudal smoothed data):%i\n', [newmasksize])
%             fprintf('The size of the intersection of the two mask :%i\n', [intersectionofthetwomask])
%             fprintf('The size of the skipped vertices from the new mask (not intersected with the original) %i and percent of those in the new mask: %.4f%%\n', [newmasksize-intersectionofthetwomask, ...
%                 round((newmasksize-intersectionofthetwomask)/newmasksize*100,2)])
%             fprintf('The percentage of the original mask which is used: (100-The size of the skipped vertices from the original mask (not intersected with the new) %i and percent of those in the original mask: %.4f%%\n', [originalmasksize-intersectionofthetwomask, ...
%                 100-round((originalmasksize-intersectionofthetwomask)/originalmasksize*100,2)])
% %             fprintf('The percent difference between  of the intersection of the two mask :%i\n', [sum(intersectedwithgroupavgmask.(myroiname))])
%         end
        %we need a function to load the 
        %this would be the visualization of the 
%         if strcmp(visualizationtype,'no') || strcmp(visualizationtype,'ROI')
%             if strcmp(visualizationtype,'no')
%                 BK_displaytxt('No visualization of the ROI was selected but we DO the sampling.')
%             elseif strcmp(visualizationtype,'ROI')
                 
%             end
%             fprintf(sprintf('Starting layer sampling...\n'));
end
%             interimdata_columns_smootheddata ={intersectedwithgroupavgmask,newmasksizes};

% %             end
% %             


%             save(smoothedmaskfilenm,'interimdata_columns_smootheddata','-v7.3');
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


%% identify clusters in the surface - theoretically this function works as it is intended to work, I double checked
function clusters=clusteridentification(triangles,markedvertices)
%
% The function iterates through the vertices which are in the mask, and
% create clusters based on the edge information. It adds a vertex to a
% cluster and screen its neightbours, and add them if they are also in the
% mask. In the next step, it screen the newly added vertex (one by one) and
% add its neighbours to the cluster if they are in the mask.
% 
% t
% inputs: 
%   - triangles: defines three numbers, which form a triangle/face in
% the surface. Basically, the neighbours of each vertex is defined this
% way.
%   - markedvertices: defines the positin in the whole mesh(using the same
%   coding as the triangles) the different vertices in the mask
    if isempty(triangles) || size(triangles,2) ~= 3
        error('TRIP must be a T-by-3 matrix of vertex indices.');
    end
    marked = unique(markedvertices(:));
    if isempty(marked)
        clusters = {};
        return;
    end
    
    % Map vertex ID -> index within the marked list (0 if not marked)
    %define the total number of vertices (that is in the whole mesh, and
    %not only under the interest). Give an index within the mask (starting from 1 to end of the mask) 
    n = max(max(triangles(:)), max(marked));
    posInMarked = zeros(n,1);
    posInMarked(marked) = 1:numel(marked);
    %create a logical vector about the marked vertices, we use this later to
    %form the clusters and change the corresponding vertex(based on its ID)
    %to true,so it is in the cluster
    visitedMarked = false(numel(marked),1);
    clusters = {};
    %we start to iterate on the marked
    for idx = 1:numel(marked)
        %check if the vertex was already inlcuded in any cluster
        if visitedMarked(idx)
            continue;
        end
        v = marked(idx);
        comp = v;
        visitedMarked(idx) = true;
        %define the vertex for cluster initialization:
        queue = v;
        head = 1;
        while head <= numel(queue)
            u = queue(head);
            head = head + 1;
    
            % Find triangles containing u
            triRows = triangles(any(triangles == u, 2), :);
            %define the neighbours based on the triangle information, and
            %remove self connection
            neigh = unique(triRows(:));
            neigh(neigh == u) = [];
    
            % Restrict to marked and not yet visited
            %find out if the neighbouring vertices are also in the mask, if
            %so update the neighbour IDs only with those which are in the
            %mask.also keep the vertex ID(neigh) and its position in the
            %mask (neighPos)
            neighPos = posInMarked(neigh);
            mask = neighPos ~= 0;         % in marked
            neigh = neigh(mask);
            neighPos = neighPos(mask);
            %check if the neightoburs wich are in the mask are still
            %unvisited, if so visit the next condition. Important! the
            %new neighbours are checked if they have been already added to the cluster, to make sure that we do not add them two times 
            addMask = ~visitedMarked(neighPos);
            toAdd = neigh(addMask);
            toAddPos = neighPos(addMask);
            %if there are neighbours which are in the mask, visit the
            %condition and:
            if ~isempty(toAdd)
                %add morevertices to the the cluster
                queue = [queue; toAdd]; %#ok<AGROW>
                %update the information about vertex visits(it will skip than the vertex in
                %the initialization of a cluster with that specific verte)
                visitedMarked(toAddPos) = true;
                %update also the indices (in the whole mesh) of the cluster
                %elements
                comp = [comp; toAdd]; %#ok<AGROW>
            end
        end
        clusters{end+1} = comp; %#ok<AGROW>
    end
end

%% identify all neighbouring vertices, and use the outputvertex indices for clusterchecking.
function newMarkedVertices=neigbourexpansion(triangles,markedvertices)
%
% The function iterates through the vertices which are in the mask, and
% expand the set of marked vertices by including all their neighbors % (vertices that share any triangle with a marked vertex).. 
% 
% t
% inputs: 
%   - triangles: defines three numbers, which form a triangle/face in
% the surface. Basically, the neighbours of each vertex is defined this
% way.
%   - markedvertices: defines the positin in the whole mesh(using the same
%   coding as the triangles) the different vertices in the mask

    if isempty(triangles) || size(triangles,2) ~= 3
        error('triangles must be a T-by-3 matrix of vertex indices.');
    end
    
    marked = unique(markedvertices(:));
    if isempty(marked)
        newMarkedVertices = marked;
        return;
    end
    
    % Find all triangles that contain at least one marked vertex
    triHasMarked = any(ismember(triangles, marked), 2);
    touchedTris = triangles(triHasMarked, :);
    
    % All vertices in those triangles are the marked vertices plus their neighbors
    neighborsPlus = unique(touchedTris(:));
    
    % Expand the mask (include original marked vertices and their neighbors)
    % If you prefer to preserve the original order of 'marked' and only append new ones,
    % use the two lines below instead of the single unique() call.
    % newNeighbors = setdiff(neighborsPlus, marked);
    % newMarkedVertices = [marked; newNeighbors];
    
    newMarkedVertices = neighborsPlus;

end

%% load layers function (checked -keep it)
function [layer_boundaries,faces,vertexcoordinates] = BK_load_layer_boundaries_sidespecific(subid,fspath)
%this is a modification of VPF's VPF_load_layer_boundaries function. The
%differences are the follwoing: it outputs the left and right hemispheres
%vertex coordinates separately, and also the faces list are outputted
%only load white matter surface faces list as it is the same for pial
%surface list (this is how vertex pairs on each side define a "colulmn")
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
    [layer_boundaries,faces.lh_w] = read_surf(fn,asciiload);
    vertexcoordinates.left=layer_boundaries;
    fn = [fspath '\' subid '\surf\rh.white']; 
    [tmp,faces.rh_w] = read_surf(fn,asciiload);
    vertexcoordinates.right=tmp;
    layer_boundaries  = cat(1,layer_boundaries ,tmp);
    
    
    fn = [fspath '\' subid '\surf\lh.pial'];
    [tmp,~] = read_surf(fn,asciiload);

    fn = [fspath '\' subid '\surf\rh.pial'];
    [maki,~]=read_surf(fn,asciiload);
    tmp = cat(1,tmp,maki);
    layer_boundaries = cat(3,layer_boundaries,tmp);

end

%% convert the ROI to surface (wb_command) - main function (checked -keep it)
%I removed the wb_command, so we do not overwrite the original files.
%Theoretically it does not, but just to being on the safe side...
%the pipeline will crash if the old one does not there...
function ROIs = BK_convert_load_ROIs_originalmask(subid,roipath,fspath,N_vertex,region,localizer)
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
%             ROIs(:,ROI) = VPF_average_mask_and_thres(data(1:end-4));
            ROIs(:,ROI) = VPF_average_mask_and_thres(data(1:end-7));
        catch %apparently, not created yet so create them first
%             hemis = {'lh','rh'};
%             for hemi = 1:2
% %                 cmd = ['wsl wb_command -volume-to-surface-mapping ' '"' wls_data '"'...
% %                     ' "' wls_fspath '/' subid '/surf/' hemis{hemi} '.white.gii' '"'...
% %                     ' "' wls_data(1:end-4) '_' hemis{hemi} '.shape.gii' '"'...
% %                     ' -cubic'];
%                 cmd = ['wsl wb_command -volume-to-surface-mapping ' '"' wls_data '"'...
%                     ' "' wls_fspath '/' subid '/surf/' hemis{hemi} '.white.gii' '"'...
%                     ' "' wls_data(1:end-7) '_' hemis{hemi} '.shape.gii' '"'...
%                     ' -cubic'];
%                 system(cmd);
%             end
%             ROIs(:,ROI) = VPF_average_mask_and_thres(data(1:end-4));
            ROIs(:,ROI) = VPF_average_mask_and_thres(data(1:end-7));
        end
    end
end

%% convert the ROI to surface (wb_command) - main function (checked -keep it)
function ROIs = BK_convert_load_ROIs(subid,roipath,fspath,N_vertex,region,localizer)
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
    
%     ROIs = zeros([N_vertex,N_total_ROIs],'logical');
    ROIs = cell(N_total_ROIs,1);
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
%             ROIs(:,ROI) = VPF_average_mask_and_thres(data(1:end-4));
%             ROIs(:,ROI) = VPF_average_mask_and_thres(data(1:end-7));
            ROIs{ROI} = BK_average_mask_and_thres(data(1:end-7));
        catch %apparently, not created yet so create them first
            hemis = {'lh','rh'};
            for hemi = 1:2
%                 cmd = ['wsl wb_command -volume-to-surface-mapping ' '"' wls_data '"'...
%                     ' "' wls_fspath '/' subid '/surf/' hemis{hemi} '.white.gii' '"'...
%                     ' "' wls_data(1:end-4) '_' hemis{hemi} '.shape.gii' '"'...
%                     ' -cubic'];
%                 cmd = ['wsl wb_command -volume-to-surface-mapping ' '"' wls_data '"'...
%                     ' "' wls_fspath '/' subid '/surf/' hemis{hemi} '.white.gii' '"'...
%                     ' "' wls_data(1:end-7) '_' hemis{hemi} '.shape.gii' '"'...
%                     ' -cubic'];
%                 system(cmd);
            end
%             ROIs(:,ROI) = VPF_average_mask_and_thres(data(1:end-4));
%             ROIs(:,ROI) = VPF_average_mask_and_thres(data(1:end-7));
            ROIs{ROI} = BK_average_mask_and_thres(data(1:end-7));
        end
    end
end

function alldata = BK_average_mask_and_thres(in)
    %this is a modification of the similarly named function from VPF. Here
    %we provide left and right hemishpere data separately.
    % loads the ROI surfaces. As they are transformed from voxel to surface space, 
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
    
%     cdata = logical(cat(1,cdata_lh,cdata_rh));
    alldata.lh=cdata_lh;
    alldata.rh=cdata_rh;
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