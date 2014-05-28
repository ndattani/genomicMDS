function plot_output = cgr_plot(st)
tic;

A=[0,0];                                  % the definition of point A
C=[0,10];                                 % the definition of point C
G=[10,10];                                % the definition of point G
T=[10,0];                                 % the definition of point T

startNuc=1;n=length(st);endNuc=n;

point=zeros(endNuc-startNuc+6,2);         % the result matrix
temp_center=(A+C+G+T)/4;                  % the initial center

point(1,:)=temp_center;
point(2,:)=A;
point(3,:)=C;
point(4,:)=G;
point(5,:)=T;

for i=startNuc:endNuc;
    switch st(i)                          % determine which letter the current base is
        case 'a'
            temp_p=A;                     % assign the corresponding vertex
        case 'g'
            temp_p=G;
        case 'c'
            temp_p=C;
        case 't'
            temp_p=T;
    end
    
    point(i+5,:)=(temp_center+temp_p)/2;  % find the result.       % IF YOU WANT TO STORE EVERYTHING 
    temp_center=point(i+5,:);             % define the new center. % IF YOU WANT TO STORE EVERYTHING 
    
    %pointTemp=(temp_center+temp_p)/2;    % find the result.       % if you DONT want to store
    %temp_center=pointTemp;               % define the new center. % if you DONT want to store
    %plot(temp_center(1),temp_center(2));                          % if you DONT want to store. You need to plot it right away, rather than at the end. Might speed up calculation and save memory for huge sequences
end

axis('off')
hold('on')

% text(-0.25,0,'A');
% text(-0.25,10,'C');
% text(10,10.05,'G');
% text(10.05,0,'T');

axis('square')
scatter(point(:,1),point(:,2),1,'black')
line(A,C,'Color','k');                     % This makes the box surrounding the plot
line(C,G,'Color','k');
line(G,T,'Color','k');
line(T,A,'Color','k');

%% Method 1
% pos = get(gcf,'Position');
% set(gcf, 'color', 'white');
% I= getframe(gcf,[0 0 pos(3) pos(4)]);

%% Method 2
% set(gca,'Position',[0 0 1 1]);set(gcf,'Position',[0 0 1 1])
% ti = get(gca,'TightInset');               % http://www.mathworks.co.uk/matlabcentral/newsreader/view_thread/314017
% set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

%% Method 3
set(gca, 'LooseInset', [0,0,0,0]);          % http://undocumentedmatlab.com/blog/axes-looseinset-property/
% set(gcf,'Visible', 'off')
set(gcf,'Units','centimeters','Position',[2 2 15 15])
I= getframe;
plot_output=I.cdata;

%%
hold('off');

% saveas(gcf,'(filename)','jpg');
% print('-djpeg', filename(1:(end)-6), '-r 150');

% print(image,'-dpng',strcat(filename(1:(end)-6),'(',num2str(startNuc),'-',num2str(endNuc),')','.png'));
% crop(strcat(filename(1:(end)-6),'(',num2str(startNuc),'-',num2str(endNuc),')','.png'),0)
clear('image')
toc;
