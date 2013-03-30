%% THIS CODE IS PROTECTED BY COPYRIGHT LAW.
%% The owners of this code are Abu Sadat Md. Sayem, Nathanlial Bryans, Nike S. Dattani, and Ronghai Tu. The copyright is currently owned by Nike S. Dattani.

%% NO PART OF THIS CODE IS TO BE MODIFIED OUTSIDE OF GIT. LEGAL CONSEQUENCES WILL APPLY
function genomicMDS(dataFile,vectorOfSheetsToUse,writeToExcel)
close('all'); % it's important to close the figures because in case hole('off') wasn't called, this will ensure we don't have the new points plotted over the old points (which can look bizarre since matlab's MDS doesn't seem to be deterministic, so plotting the same dataset twice could give two different things on top of each other
% clearvars -except dataFile firstFigureNumber lastFigureNumber; % might want to do because datafile might be 2GB and might be in the workspace
% disp(['All variables except for allSpecies cleared !']);
%% A file containing all Matlab's information on each accession number
if isstruct(dataFile);
    allSpecies=dataFile;clear('dataFile'); % don't need to worry about memory duplication because it does 'lazy copying'
else
    tic;load(dataFile);elapsedTimeLoadStructureArray=toc;  % file must contain the SSIM distances between all organisms
    disp(['Loading the gene structure array is done ! It took: ' num2str(elapsedTimeLoadStructureArray) ' seconds']);
end
%% Each sheet of the excel file generates a different figure, with different taxa included
for sheetIndex=1:length(vectorOfSheetsToUse);sheet=vectorOfSheetsToUse(sheetIndex);
    disp(['FIGURE  ' num2str(sheet) ' !!!'])
    clearvars -except sheet allSpecies writeToExcel vectorOfSheetsToUse
    tic;[whetherOrNotToPlotNumbers,taxaToInclude,~]=xlsread('setsOfTaxaAndColors.xlsx',sheet);elapsedTimeLoadExcelSheet=toc;
    disp(['Loading the excel sheet containing the specifications for figure ' num2str(sheet) ' is done ! It took: ' num2str(elapsedTimeLoadExcelSheet) ' seconds']);
    taxaToInclude=taxaToInclude';
    %% If there's a loaded data set, here we filter out the ones included in taxaToInclude. If there's no loaded data set, here we get the appropriate gene structure arrays using GETGENBANK
    filteredDataSet=zeros(length(allSpecies),length(taxaToInclude));numberOfOrganismsInEachTaxon=zeros(size(taxaToInclude,2),1);numberOfOrganismsInEachTaxonGroupForLegend=zeros(length(unique(taxaToInclude(2,:))),1);shouldWeKeepThisOrganism=0;chosenOrganismIndex=1;
    warning off all
    for i=1:length(allSpecies) %if this loop is slow, consider preallocating the structure. according to this http://blogs.mathworks.com/loren/2008/02/01/structure-initialization/ perhaps the best way to do this is the do the loop backwards, so that the last entry of gene is defined first, that way the hwole structure is built in advance.
        %temp=getgenbank(allSpecies{i});
        temp=allSpecies(i); % if data set is preloaded
        temp2=temp.SourceOrganism';
        for j=1:size(taxaToInclude,2) % using length will run into problems when there' less than three taxa to include, because taxaToInclude will have more rows than columns.
            if isequal(taxaToInclude{1},'all')==0 % If the excel file's entry is 'all' then we want ALL organisms, so we just skip to the ELSE statement which automatically includes the organism
                shouldWeKeepThisOrganism=shouldWeKeepThisOrganism+mod(isempty(strfind(temp2(1:end),taxaToInclude{1,j}))+1,2); %if isempty is 1, then the organism is not in one of the desired taxa, so we want 'shouldWeKeepThisOrganism' to remain as 0, so we turn 1 into 0 by adding 1 (to make it 2) then modding it with 2 (to get 0), if isempty gives 0, then we turn it into 1 by adding 1 to it, then modding it with 2 which keeps it at one. We could also just subtract 1 in both cases, and then do shouldWeKeepThisOrganism=abs(shouldWeKeepThisOrganism+(isEmpty... -1)
            else % why do we need a special case for ALL though ? can't we just make geneData=allSpecies ?
                shouldWeKeepThisOrganism=1; % if we want ALL organisms then we obviously want shouldWeKeepThisOrganism = 1
            end
            if shouldWeKeepThisOrganism==1;
                indices(chosenOrganismIndex)=i;    geneData(chosenOrganismIndex)=temp;
                numberOfOrganismsInEachTaxon(j)=numberOfOrganismsInEachTaxon(j)+1;  chosenOrganismIndex=chosenOrganismIndex+1;   shouldWeKeepThisOrganism=0;
            end
        end
    end
    totalNumberOfOrganisms=i;
    
    taxaGroupsToIncludeInLegend=unique(taxaToInclude(2,:),'stable'); % without setOrder='stable' , they come out alphabetical, and you'll get amphibians before fish, in the legend, despite it being the other way around in the excel file
    for i=1:length(numberOfOrganismsInEachTaxon)
        for j=1:length(numberOfOrganismsInEachTaxonGroupForLegend)
            if strcmp(taxaToInclude(2,i),taxaGroupsToIncludeInLegend(j)); numberOfOrganismsInEachTaxonGroupForLegend(j)=numberOfOrganismsInEachTaxonGroupForLegend(j)+numberOfOrganismsInEachTaxon(i); end
        end
    end
    
    disp(['Total number of organisms in unfiltered dataset:' num2str(totalNumberOfOrganisms)])
    %% Make a gene structure array called geneData which is the relavent subset of the gene structure array allSpecies, and assign the appopriate colors to organisms according to their taxa
    for i=1:length(geneData) % can't include this in the loop above because some structures will already have the fileds 'taxonForOurPurposes' and 'color' but we'll be assigning geneData(i+1)=allSpecies(j) where allSpecies(j) does not have those fields, so we'll get the error "subscripted assignment between dissimilar structures", look into a way to fix this! Perhaps just add those fields as empty shit to the whole allOrganisms structure
        for k=1:size(taxaToInclude,2) % using length will run into problems when there' less than three taxa to include, because taxaToInclude will have more rows than columns.
            temp=geneData(i).SourceOrganism';
            if mod(isempty(strfind(temp(1:end),taxaToInclude{1,k}))+1,2)
                geneData(i).taxonForOurPurposes=taxaToInclude{1,k};geneData(i).taxonForLegend=taxaToInclude{2,k};geneData(i).color=taxaToInclude{3,k};geneData(i).MarkerEdgeColor=taxaToInclude{4,k};
                geneData(i).whetherOrNotToPlotNumber=whetherOrNotToPlotNumbers(k);
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
    disp(['Creating the ssim distance matrix for figure ' num2str(sheet) ' is done !']);
    %% Classical Multidimensional Scaling
    tic;Y = cmdscale(squareform(ssimDistances));elapsedTimeMDS=toc;
    disp(['MDS calculation for figure ' num2str(sheet) ' is done ! It took: ' num2str(elapsedTimeMDS) ' seconds']);
    %% Rotations and Scaling
    desiredSize=1;
    desiredSizeOfBorder=desiredSize*1.1;
    
    %Flip the grid in the y direction (around the x axis)
    %Y(:,2)=-Y(:,2);
    
    %Flip the grid in the x direction (around the y axis)
    %Y(:,1)=-Y(:,1);
        
    currentLargestSizeY=max(Y(:,2));currentLargestSizeX=max(Y(:,1));
    currentSmallestSizeY=min(Y(:,2));currentSmallestSizeX=min(Y(:,1));
    
    %Used to readjust the scale for the y axis
    % offsetY = desiredSize/currentLargestSizeY;
    % Y(:,2)=Y(:,2)*offsetY;
    Y(:,2)=2*(Y(:,2)-currentSmallestSizeY)/(currentLargestSizeY-currentSmallestSizeY)-1; % scale y-values to be between -1 and 1
    
    %Used to readjust the scale for the x axis
    % offsetX = desiredSize/currentLargestSizeX;
    % Y(:,1)=Y(:,1)*offsetX;
    Y(:,1)=2*(Y(:,1)-currentSmallestSizeX)/(currentLargestSizeX-currentSmallestSizeX)-1; % scale x-values to be between -1 and 1
    disp(['The MDS data for figure ' num2str(sheet) ' has been rescaled for optimal display!']);
    %% Make the figure
    tic;
    figure(sheet);
    %set(figure(sheet),'position',[200 200 500 500]);
    hold('on')
    
    %organismsToExcludeFromPlot=506:562; % combined with line 94 to remove organisms form plot
    for i = 1:length(geneData)
        %if ~ismember(indices(i),organismsToExcludeFromPlot) % combined with line 92 to remove organisms from plot 
        for j=1:length(taxaGroupsToIncludeInLegend)
            if strcmp(geneData(i).taxonForLegend,taxaGroupsToIncludeInLegend(j))
                plotHandle(j)=plot(Y(i,1),Y(i,2),'o','MarkerSize',12,'MarkerFaceColor',rgb(geneData(i).color),'MarkerEdgeColor',rgb(geneData(i).MarkerEdgeColor),'LineWidth',0.5); % with black outline
%                 for ii=1:length(unique(taxaToInclude{2,j})) % # of unique groups for legend entries
%                     if strcmp(geneData(i).taxonForLegend,taxaToInclude{2,j})
%                 set(plotHandle(j),'Parent',childrenOfAxesForGroupingLegendEntries())
%                     end
%                end
                if geneData(i).whetherOrNotToPlotNumber==1;text(Y(i,1)+0.005, Y(i,2), int2str(indices(i)), 'HorizontalAlignment', 'left', 'Color', rgb(geneData(i).color));end;
                %plotHandle(j)=plot3(Y(i,1),Y(i,2),Y(i,3),strcat(geneData(i).color,'o'),'MarkerSize',6,'MarkerFaceColor',geneData(i).color);
            end
        end
    end
axis([-desiredSizeOfBorder desiredSizeOfBorder -desiredSizeOfBorder desiredSizeOfBorder]);axis square
%axis([-0.65 -0.3 -0.55 -0.25]);axis square

%     fill([-1.13 -1.13 -0.995 -0.995],[-1.13 -0.995 -0.995 -1.13],'k') % bottom left corner
%     fill([-1.13 -1.13 -0.995 -0.995],[0.995  1.13   1.13  0.995],'k') % top left corner
%     fill([0.995 0.995  1.13    1.13],[-1.13 -0.995 -0.995 -1.13],'k') % bottom right corner
%     fill([0.995 0.995  1.13    1.13],[0.995  1.13   1.13  0.995],'k') % top right corner
    hold('on');box('on');
    set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',3,'FontSize',25,'TickDir','out','TickLength',[0.03 0.03],'XTick',-1:0.2:1,'YTick',-1:0.2:1,'Position',[0.129375 0.0717437722419927 0.775 0.815],'FontWeight','light');
    
    legendHandle=legend(plotHandle,strcat('\textbf{',taxaGroupsToIncludeInLegend.',{' ('}, cellstr(num2str(numberOfOrganismsInEachTaxonGroupForLegend)),{')'},'}'));
    set(legendHandle,'Interpreter','latex','FontSize',30,'LineWidth',3,'Position',[0.747672758188061 0.753355153875044 0.134706814580032 0.139318885448916]);
    axesHandle=copyobj(gca,gcf);set(axesHandle,'xaxislocation','top','yaxislocation','right');
    legendHandleChildren=get(legendHandle,'children');set(legendHandleChildren(1:3:end),'MarkerSize',20) 
    uistack(plotHandle,'bottom');uistack(legendHandle,'top'); 
    hold('off');elapsedTimeFigure=toc;
    disp(['Plotting the figure ' num2str(sheet) ' is done ! It took: ' num2str(elapsedTimeFigure) ' seconds']);
    maximumError=max(abs(squareform(ssimDistances) - pdist(Y(:,1:2))));averageError=sum(abs(squareform(ssimDistances) - pdist(Y(:,1:2))))/(length(squareform(ssimDistances)));averageErrorAsPercentageOfTheMaximumPossibleLineLengthInAsquareOfLength2=(averageError/(2*sqrt(2)))*100; %longest possible line in a square of length 2 is 2*sqrt(2)
    disp(['The Maximum Error for the plot is ' num2str(maximumError)]);disp(['The Average Error for the plot is ' num2str(averageError)]);disp(['The Percentage of Error for the plot is ' num2str(averageErrorAsPercentageOfTheMaximumPossibleLineLengthInAsquareOfLength2)]);
    %% Print the excel file with information about each species   
    if writeToExcel==1; % it will take ages if you're working on a remote machine
    tic;
    columnHeader = {'NUMBER','X-COORDINATES', 'Y-CORDINTATES','ACCESSION NUMBER','NAME','SEQUENCE LENGTH','COLOR','MARKER-EDGE-COLOR','TAXA'};
    
    Species_Information_For_Spreadsheet=cell(length(geneData),27); % 27 is kind of arbitrary. I'm just hoping to make sure I have enough columns for all the taxonomic categories, some organisms have more of these than others and I don't know what the largest number of them is
    for i=1:length(geneData)
        Species_Information_For_Spreadsheet(i,1)=cellstr(int2str(i));                               % Numerical identity
        Species_Information_For_Spreadsheet(i,4)=cellstr(geneData(i).Accession);                    % Accession Number
        Species_Information_For_Spreadsheet(i,5)=cellstr(geneData(i).SourceOrganism(1,1:end));      % Name of species
        Species_Information_For_Spreadsheet(i,6)=cellstr(geneData(i).LocusSequenceLength);          % Sequence length
        Species_Information_For_Spreadsheet(i,7)=cellstr(geneData(i).color);                        % Color
        Species_Information_For_Spreadsheet(i,7)=cellstr(geneData(i).MarkerEdgeColor);              % Marker edge color
        Species_Taxa=strtrim(regexp(sprintf((geneData(i).SourceOrganism(2:end,:))'),';','split'));  % turn SourceOrganism for species i into a 1xN cell array of taxon names
        for j=1:length(Species_Taxa) 
        Species_Information_For_Spreadsheet(i,7+j)=Species_Taxa(j);                                 
        end
    end
    
    Species_Information_For_Spreadsheet(:,2)=cellstr(num2str(Y(:,1)));                                                % X-coordinate
    Species_Information_For_Spreadsheet(:,3)=cellstr(num2str(Y(:,2)));                                                % Y-coordinate
       
    xlswrite('Dataset.xlsx', columnHeader ,sheet,'A1'); % if you get the error:  Error using dlmwrite (line 118) The input cell array cannot be converted to a matrix. It means you don't have windows, don't have excel, or have the excel starter version
    xlswrite('Dataset.xlsx', Species_Information_For_Spreadsheet,sheet,'A2');
    
    elapsedTimeWritingToExcelFile=toc;
    disp(['Writing to the excel file ' num2str(sheet) ' is done ! It took: ' num2str(elapsedTimeWritingToExcelFile) ' seconds']);
    end
end % loop through each sheet (each choice of taxa for a particular figure) of excel file
end % function (this end statement is not necessary, but once when I didn't have it there and I had an extra end statement in the file, I didn't get an error because it assigned that extra end statement to end the function)

%% COLOR DEFINITIONS   Taken from Kristjan Jonasson's FEX "RGB triple of color name, version 2": http://www.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name-version-2, with ReptileGreen, BirdOrange and InsectBrown added manually by Nike Dattani on 2013/03/27/6:40PM GMT
function rgb = rgb(s)
  persistent num name
  if isempty(num) % First time rgb is called
    [num,name] = getcolors();
    name = lower(name);
    num = reshape(hex2dec(num), [], 3);
    % Divide most numbers by 256 for "aesthetic" reasons (green=[0 0.5 0])
    I = num < 240;  % (interpolate F0--FF linearly from 240/256 to 1.0)
    num(I) = num(I)/256;
    num(~I) = ((num(~I) - 240)/15 + 15)/16; + 240;
  end
  if strcmpi(s,'chart')
    showcolors()
  else
    k = find(strcmpi(s, name));
    if isempty(k)
      error(['Unknown color: ' s]);
    else
      rgb = num(k(1), :);
    end
  end
end
function [hex,name] = getcolors()
  css = {
    %White colors
    'FF','FF','FF', 'White'
    'FF','FA','FA', 'Snow'
    'F0','FF','F0', 'Honeydew'
    'F5','FF','FA', 'MintCream'
    'F0','FF','FF', 'Azure'
    'F0','F8','FF', 'AliceBlue'
    'F8','F8','FF', 'GhostWhite'
    'F5','F5','F5', 'WhiteSmoke'
    'FF','F5','EE', 'Seashell'
    'F5','F5','DC', 'Beige'
    'FD','F5','E6', 'OldLace'
    'FF','FA','F0', 'FloralWhite'
    'FF','FF','F0', 'Ivory'
    'FA','EB','D7', 'AntiqueWhite'
    'FA','F0','E6', 'Linen'
    'FF','F0','F5', 'LavenderBlush'
    'FF','E4','E1', 'MistyRose'
    %Grey colors'
    '80','80','80', 'Gray'
    'DC','DC','DC', 'Gainsboro'
    'D3','D3','D3', 'LightGray'
    'C0','C0','C0', 'Silver'
    'A9','A9','A9', 'DarkGray'
    '69','69','69', 'DimGray'
    '77','88','99', 'LightSlateGray'
    '70','80','90', 'SlateGray'
    '2F','4F','4F', 'DarkSlateGray'
    '00','00','00', 'Black'
    %Red colors
    'FF','00','00', 'Red'
    'FF','A0','7A', 'LightSalmon'
    'FA','80','72', 'Salmon'
    'E9','96','7A', 'DarkSalmon'
    'F0','80','80', 'LightCoral'
    'CD','5C','5C', 'IndianRed'
    'DC','14','3C', 'Crimson'
    'B2','22','22', 'FireBrick'
    '8B','00','00', 'DarkRed'
    %Pink colors
    'FF','C0','CB', 'Pink'
    'FF','B6','C1', 'LightPink'
    'FF','69','B4', 'HotPink'
    'FF','14','93', 'DeepPink'
    'DB','70','93', 'PaleVioletRed'
    'C7','15','85', 'MediumVioletRed'
    %Orange colors
    'FF','A5','00', 'Orange'
    'FF','8C','00', 'DarkOrange'
    'FF','7F','50', 'Coral'
    'FF','63','47', 'Tomato'
    'FF','45','00', 'OrangeRed'
    'FF','88','00', 'BirdOrange'
    %Yellow colors
    'FF','FF','00', 'Yellow'
    'FF','FF','E0', 'LightYellow'
    'FF','FA','CD', 'LemonChiffon'
    'FA','FA','D2', 'LightGoldenrodYellow'
    'FF','EF','D5', 'PapayaWhip'
    'FF','E4','B5', 'Moccasin'
    'FF','DA','B9', 'PeachPuff'
    'EE','E8','AA', 'PaleGoldenrod'
    'F0','E6','8C', 'Khaki'
    'BD','B7','6B', 'DarkKhaki'
    'FF','D7','00', 'Gold'
    %Brown colors
    'A5','2A','2A', 'Brown'
    'FF','F8','DC', 'Cornsilk'
    'FF','EB','CD', 'BlanchedAlmond'
    'FF','E4','C4', 'Bisque'
    'FF','DE','AD', 'NavajoWhite'
    'F5','DE','B3', 'Wheat'
    'DE','B8','87', 'BurlyWood'
    'D2','B4','8C', 'Tan'
    'BC','8F','8F', 'RosyBrown'
    'F4','A4','60', 'SandyBrown'
    'DA','A5','20', 'Goldenrod'
    'B8','86','0B', 'DarkGoldenrod'
    'CD','85','3F', 'Peru'
    'D2','69','1E', 'Chocolate'
    '8B','45','13', 'SaddleBrown'
    'A0','52','2D', 'Sienna'
    '80','00','00', 'Maroon'
    '66','33','00', 'InsectBrown'
    %Green colors
    '00','80','00', 'Green'
    '98','FB','98', 'PaleGreen'
    '90','EE','90', 'LightGreen'
    '9A','CD','32', 'YellowGreen'
    'AD','FF','2F', 'GreenYellow'
    '7F','FF','00', 'Chartreuse'
    '7C','FC','00', 'LawnGreen'
    '00','FF','00', 'Lime'
    '32','CD','32', 'LimeGreen'
    '00','FA','9A', 'MediumSpringGreen'
    '00','FF','7F', 'SpringGreen'
    '66','CD','AA', 'MediumAquamarine'
    '7F','FF','D4', 'Aquamarine'
    '20','B2','AA', 'LightSeaGreen'
    '3C','B3','71', 'MediumSeaGreen'
    '2E','8B','57', 'SeaGreen'
    '8F','BC','8F', 'DarkSeaGreen'
    '22','8B','22', 'ForestGreen'
    '00','64','00', 'DarkGreen'
    '6B','8E','23', 'OliveDrab'
    '80','80','00', 'Olive'
    '55','6B','2F', 'DarkOliveGreen'
    '00','80','80', 'Teal'
    '00','BB','00', 'ReptileGreen'
    '00','FF','CC', 'AmphibianTurquoise'
    
    %Blue colors
    '00','00','FF', 'Blue'
    'AD','D8','E6', 'LightBlue'
    'B0','E0','E6', 'PowderBlue'
    'AF','EE','EE', 'PaleTurquoise'
    '40','E0','D0', 'Turquoise'
    '48','D1','CC', 'MediumTurquoise'
    '00','CE','D1', 'DarkTurquoise'
    'E0','FF','FF', 'LightCyan'
    '00','FF','FF', 'Cyan'
    '00','FF','FF', 'Aqua'
    '00','8B','8B', 'DarkCyan'
    '5F','9E','A0', 'CadetBlue'
    'B0','C4','DE', 'LightSteelBlue'
    '46','82','B4', 'SteelBlue'
    '87','CE','FA', 'LightSkyBlue'
    '87','CE','EB', 'SkyBlue'
    '00','BF','FF', 'DeepSkyBlue'
    '1E','90','FF', 'DodgerBlue'
    '64','95','ED', 'CornflowerBlue'
    '41','69','E1', 'RoyalBlue'
    '00','00','CD', 'MediumBlue'
    '00','00','8B', 'DarkBlue'
    '00','00','80', 'Navy'
    '19','19','70', 'MidnightBlue'
    %Purple colors
    '80','00','80', 'Purple'
    'E6','E6','FA', 'Lavender'
    'D8','BF','D8', 'Thistle'
    'DD','A0','DD', 'Plum'
    'EE','82','EE', 'Violet'
    'DA','70','D6', 'Orchid'
    'FF','00','FF', 'Fuchsia'
    'FF','00','FF', 'Magenta'
    'BA','55','D3', 'MediumOrchid'
    '93','70','DB', 'MediumPurple'
    '99','66','CC', 'Amethyst'
    '8A','2B','E2', 'BlueViolet'
    '94','00','D3', 'DarkViolet'
    '99','32','CC', 'DarkOrchid'
    '8B','00','8B', 'DarkMagenta'
    '6A','5A','CD', 'SlateBlue'
    '48','3D','8B', 'DarkSlateBlue'
    '7B','68','EE', 'MediumSlateBlue'
    '4B','00','82', 'Indigo'
    %Gray repeated with spelling grey
    '80','80','80', 'Grey'
    'D3','D3','D3', 'LightGrey'
    'A9','A9','A9', 'DarkGrey'
    '69','69','69', 'DimGrey'
    '77','88','99', 'LightSlateGrey'
    '70','80','90', 'SlateGrey'
    '2F','4F','4F', 'DarkSlateGrey'
    };
  hex = css(:,1:3);
  name = css(:,4);
end
%% ALL FEX FILES FO CONVERTING TICK LABELS TO LATEX ARE CURRENTLY NOT USER-FRIENDLY ENOUGH FOR THIS PROGRAM, SO WE AWAIT MATLAB BUILDING IN THIS FUNCTION