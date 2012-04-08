%% THIS CODE IS PROTECTED BY COPYRIGHT LAW. 
%% The owners of this code are Abu Sadat Md. Sayem, Nathanlial Bryans, Nike
%% S. Dattani, and Ronghai Tu. The copyright is currently owned by Nike S. Dattani. 

%% NO PART OF THIS CODE IS TO BE MODIFIED OUTSIDE OF GIT.
%%
clear all;
%%
outputFileName='datadatadatarenamelaytah';
%%
tic;accessionNumbers=importdata('2012-03-17-accessionNumbers.txt');elapsedTimeLoadAccessionNumbers=toc;  %A file containing: Line 1: the http://www.ncbi.nlm.nih.gov/nuccore search query used to get this dataset, Lines 2 to end: the accession numbers
disp(['Loading the file containing the accession numbers is done ! It took: ' num2str(elapsedTimeLoadAccessionNumbers) ' seconds']);
%%
warning off all
tic;
for i=1:length(accessionNumbers)-1 %if this loop is slow, consider preallocating the structure. according to this http://blogs.mathworks.com/loren/2008/02/01/structure-initialization/ perhaps the best way to do this is the do the loop backwards, so that the last entry of gene is defined first, that way the hwole structure is built in advance.
    allSpecies(i)=getgenbank(accessionNumbers{i+1}); %the plus 1 is due to the fact that the first cell element of accessionNumbers isn't an accession number but just the search query that was used to get these accession numbers
    disp(['Data has now been loaded for accession number ' num2str(i) ' of ' num2str(length(accessionNumbers)-1)]); 
end
save(outputFileName,'allSpecies'); 
elapsedTimeGetGeneBankInformationForAllAccessionNumbers=toc;
disp(['Getting gene bank information for all accession numbers is done ! It took: ' num2str(elapsedTimeGetGeneBankInformationForAllAccessionNumbers) ' seconds']);
%% 12 minutes for 578 species: NEED THE PICTURES TO BE VISIBLE, THEREFORE THE WINDOWS HAVE TO BE OPEN AND THE MONITOR HAS TO BE ON, IN ORDER FOR THIS TO WORK !
allSpecies(1).CGR2dPixelMap=[];
%matlabpool('open','12') %lines below didn't work with parfor because the figure didn't pop up, havne't tried to fix it
for i=1:length(allSpecies)
    clf; % commenting this slows things down
    allSpecies(i).CGR2dPixelMap= cgr_plot(allSpecies(i).Sequence);
    %     set(gcf,'Visible', 'off'); % commenting this and the one below speeds things up
    %     figure('Visible', 'off');
    allSpecies(i).CGR2dPixelMap = rgb2gray(allSpecies(i).CGR2dPixelMap); 
    %if mod(i,500)==0;save('dataForParticularDataSet.mat','allSpecies');end %saving the entire array takes less than a minute, so it's not a major speed issue to be doing this in case matlab crashes midway
disp(['Getting the pixel map for the sequence of accession number ' num2str(i) ' is done !']);
end
%matlabpool('close')
save(outputFileName,'allSpecies'); 
%% Calculate the SSIM distance matrix
ssimDistances = zeros(length(allSpecies));
matlabpool('open','8')
tic;
%allSpeciesPersistant=WorkerObjWrapper(allSpecies);
parfor idx= 1:length(allSpecies)^2
     [i, j] = ind2sub([length(allSpecies), length(allSpecies)], idx);
        if j > i
              ssimDistances(idx)=1-ssim_index(allSpecies(i).CGR2dPixelMap,allSpecies(j).CGR2dPixelMap);
        end
        strcat('(i,j)=',num2str(i),',',num2str(j))
end
toc;
matlabpool('close')
ssimDistances=ssimDistances+ssimDistances.';
%%
for i =1:length(allSpecies)
    allSpecies(i).ssimDistances=ssimDistances(i,:);
end
%% If ssimDistances field of allSpecies is already filled in - currently quite memory intensive because each worker grabs the 2GB structure array, delete things we don't need, like allSpecies.sequence
ssimDistances=zeros(length(indices));
for i=1:length(indices)
    ssimDistances(i,:)=allSpecies(indices(i)).ssimDistances(indices);
end
% Or if if you have matrix already:
%ssimDistances(indices,indices)=ssim(indices,indices);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create the 2D plot from the distanceonly.txt file

%Y will contain the points calculated by CMDSCALE
tic;Y = cmdscale(squareform(ssimDistances));toc;
%%
desiredSize=1;
currentLargestSizeY=max(Y(:,2));currentLargestSizeX=max(Y(:,1));
currentSmallestSizeY=min(Y(:,2));currentSmallestSizeX=min(Y(:,1));

%Flip the grid in the y direction (around the x axis)
%Y(:,2)=-Y(:,2);

%Used to readjust the scale for the y axis
% offsetY = desiredSize/currentLargestSizeY;
% Y(:,2)=Y(:,2)*offsetY;
Y(:,2)=2*(Y(:,2)-currentSmallestSizeY)/(currentLargestSizeY-currentSmallestSizeY)-1;

%Flip the grid in the x direction (around the y axis)
%Y(:,1)=-Y(:,1);

%Used to readjust the scale for the x axis
% offsetX = desiredSize/currentLargestSizeX;
% Y(:,1)=Y(:,1)*offsetX;
Y(:,1)=2*(Y(:,1)-currentSmallestSizeX)/(currentLargestSizeX-currentSmallestSizeX)-1;
%%
tic;
hold('on')
clear('plotHandle');
    for i = 1:length(allSpecies)
%for i = 1:15
    for j=1:length(taxaToInclude)
%        for j=[1 4 5]
        if strcmp(allSpecies(i).taxonForLegend,taxaToInclude{2,j})
           plotHandle(j)=plot(Y(i,1),Y(i,2),'o','MarkerSize',6,'MarkerFaceColor',allSpecies(i).color,'MarkerEdgeColor','k'); % with black outline
            %plotHandle(j)=plot3(Y(i,1),Y(i,2),Y(i,3),strcat(allSpecies(i).color,'o'),'MarkerSize',6,'MarkerFaceColor',allSpecies(i).color);
            end
        end
    end
axis([-desiredSize desiredSize -desiredSize desiredSize]);axis square
hold('on');box('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',3,'FontSize',16);
legendHandle=legend(plotHandle,taxaToInclude{2,:});set(legendHandle,'Interpreter','latex','FontSize',16,'LineWidth',3,'Position',[0.747672758188061 0.753355153875044 0.134706814580032 0.139318885448916])
toc;

