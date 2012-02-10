%% THIS CODE IS PROTECTED BY COPYRIGHT LAW.
%% The owners of this code are Abu Sadat Md. Sayem, Nathanlial Bryans, Nike
%% S. Dattani, and Ronghai Tu. The copyright is currently owned by Nike S. Dattani.

%% NO PART OF THIS CODE IS TO BE MODIFIED OUTSIDE OF GIT. LEGAL
%% CONSEQUENCES WILL APPLY 
%%
clear all;
%% A file containing all Matlab's information on each accession number
tic;load('data.mat');elapsedTimeLoadStructureArray=toc;  % file must contain the SSIM distances between all organisms
disp(['Loading the gene structure array is done ! It took: ' num2str(elapsedTimeLoadStructureArray) ' seconds']);
%% Each sheet of the excel file generates a different figure, with different taxa included
for sheet=1:5
disp(['FIGURE  ' num2str(sheet) ' !!!'])   
    clearvars -except sheet allSpecies
    tic;[~,taxaToInclude,~]=xlsread('setsOfTaxaAndColors.xlsx',sheet);elapsedTimeLoadExcelSheet=toc;
disp(['Loading the excel sheet containing the specifications for figure ' num2str(sheet) ' is done ! It took: ' num2str(elapsedTimeLoadExcelSheet) ' seconds']);
    taxaToInclude=taxaToInclude';
    colors=taxaToInclude(3,:);
    %% If there's a loaded data set, here we filter out the ones included in taxaToInclude. If there's no loaded data set, here we get the appropriate gene structure arrays using GETGENBANK    
    filteredDataSet=zeros(length(allSpecies),length(taxaToInclude));shouldWeKeepThisOrganism=1;chosenOrganismIndex=1;
    warning off all
    for i=1:length(allSpecies) %if this loop is slow, consider preallocating the structure. according to this http://blogs.mathworks.com/loren/2008/02/01/structure-initialization/ perhaps the best way to do this is the do the loop backwards, so that the last entry of gene is defined first, that way the hwole structure is built in advance.
        %temp=getgenbank(allSpecies{i});
        temp=allSpecies(i); % if data set is preloaded
        temp2=temp.SourceOrganism';
        for j=1:length(taxaToInclude)
            if isequal(taxaToInclude{1},'all')==0
                shouldWeKeepThisOrganism=shouldWeKeepThisOrganism+mod(isempty(strfind(temp2(1:end),taxaToInclude{1,j}))+1,2); %if isempty is 1, then the organism is not in one of the desired taxa, so we want 'shouldWeKeepThisOrganism' to remain as 0, so we turn 1 into 0 by adding 1 (to make it 2) then modding it with 2 (to get 0), if isempty gives 0, then we turn it into 1 by adding 1 to it, then modding it with 2 which keeps it at one. We could also just subtract 1 in both cases, and then do shouldWeKeepThisOrganism=abs(shouldWeKeepThisOrganism+(isEmpty... -1)
            else
                shouldWeKeepThisOrganism=1;
            end
        end
        if shouldWeKeepThisOrganism==1;indices(chosenOrganismIndex)=i;geneData(chosenOrganismIndex)=temp;speciesDataSet{chosenOrganismIndex}=geneData(chosenOrganismIndex).Accession;chosenOrganismIndex=chosenOrganismIndex+1;end
        shouldWeKeepThisOrganism=0;
    end
disp(['Total number of organisms in unfiltered dataset:' num2str(i)]) 
    %% Make a gene structure array called geneData which is the relavent subset of the gene structure array allSpecies, and assign the appopriate colors to organisms according to their taxa
    for i=1:length(geneData) % can't include this in the loop above because some structures will already have the fileds 'taxonForOurPurposes' and 'color' but we'll be assigning geneData(i+1)=allSpecies(j) where allSpecies(j) does not have those fields, so we'll get the error "subscripted assignment between dissimilar structures", look into a way to fix this! Perhaps just add those fields as empty shit to the whole allOrganisms structure
        for k=1:length(taxaToInclude)
            temp=geneData(i).SourceOrganism';
            if mod(isempty(strfind(temp(1:end),taxaToInclude{1,k}))+1,2)
                geneData(i).taxonForOurPurposes=taxaToInclude{1,k};geneData(i).color=colors{k};geneData(i).taxonForLegend=taxaToInclude{2,k};
            end
        end
    end
disp(['Total number of organisms in filtered dataset:' num2str(i)])
    %% Make an ssim distance matrix, using the ssim distances stored in the structure array allSpecies
    ssimDistances=zeros(length(indices));
    for i=1:length(indices)
        ssimDistances(i,:)=allSpecies(indices(i)).ssimDistances(indices); % could probably also do with the structure array geneData right ?
    end
    % Or if if you have matrix for ALL organisms already:
    %ssimDistances(indices,indices)=ssim(indices,indices);
    %% Classical Multidimensional Scaling
    tic;Y = cmdscale(squareform(ssimDistances));elapsedTimeMDS=toc;
disp(['MDS calculation is done ! It took: ' num2str(elapsedTimeMDS) ' seconds']);
    %% Rotations and Scaling
    desiredSize=1;
    currentLargestSizeY=max(Y(:,2));currentLargestSizeX=max(Y(:,1));
    currentSmallestSizeY=min(Y(:,2));currentSmallestSizeX=min(Y(:,1));
    
    %Flip the grid in the y direction (around the x axis)
    %Y(:,2)=-Y(:,2);
    
    %Used to readjust the scale for the y axis
    % offsetY = desiredSize/currentLargestSizeY;
    % Y(:,2)=Y(:,2)*offsetY;
    Y(:,2)=2*(Y(:,2)-currentSmallestSizeY)/(currentLargestSizeY-currentSmallestSizeY)-1; % scale y-values to be between -1 and 1
    
    %Flip the grid in the x direction (around the y axis)
    %Y(:,1)=-Y(:,1);
    
    %Used to readjust the scale for the x axis
    % offsetX = desiredSize/currentLargestSizeX;
    % Y(:,1)=Y(:,1)*offsetX;
    Y(:,1)=2*(Y(:,1)-currentSmallestSizeX)/(currentLargestSizeX-currentSmallestSizeX)-1; % scale x-values to be between -1 and 1
    %% Make the figure
    tic;
    figure(sheet);
    hold('on')
    for i = 1:length(geneData)
        %for i = 1:15
        for j=1:length(taxaToInclude)
            %        for j=[1 4 5]
            if strcmp(geneData(i).taxonForLegend,taxaToInclude{2,j})
                plotHandle(j)=plot(Y(i,1),Y(i,2),'o','MarkerSize',6,'MarkerFaceColor',geneData(i).color,'MarkerEdgeColor','k'); % with black outline
                %plotHandle(j)=plot3(Y(i,1),Y(i,2),Y(i,3),strcat(geneData(i).color,'o'),'MarkerSize',6,'MarkerFaceColor',geneData(i).color);
            end
        end
    end
    axis([-desiredSize desiredSize -desiredSize desiredSize]);axis square
    hold('on');box('on');
    set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',3,'FontSize',16);
    legendHandle=legend(plotHandle,taxaToInclude{2,:});set(legendHandle,'Interpreter','latex','FontSize',16,'LineWidth',3,'Position',[0.747672758188061 0.753355153875044 0.134706814580032 0.139318885448916])
    elapsedTimeFigure=toc;
disp(['Plotting the figure ' num2str(sheet) ' is done ! It took: ' num2str(elapsedTimeFigure) ' seconds']);
end % loop through each sheet (each choice of taxa for a particular figure) of excel file
