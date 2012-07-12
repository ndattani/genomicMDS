function plot_output = cgr_plot(st)
% Generates a 2D CGR image from downloaded fasta data
% version 1.3
% revised on 2009/02/18 by Ronghai
% revised on 2009/05/26 by Nike
% revised on 2011/06/09 by Nathaniel
% Main function  read a fasta file according to a filename  and
% generates a 2D CGR image of a fasta DNA file. 

A=[0,0];     %the definition of point A
C=[0,10];    %the definition of point C
G=[10,10];   %the definition of point G
T=[10,0];    %the definition of point T    
    

% location=fullfile(directory);
% fid=fopen(fullfile(location,filename));
% st=fread(fid,'uint8=>char');  % read the strings
% 
% %if st(end)= blank
%     st(end)=[]; %remove last spot which is blank if you joined the sequence by %join! , which is the case for everything in chromosomal folder
% 
% fclose(fid);                  %close the file number
% 
% [m,n]=size(st);  %read the size of string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%st = org.Sequence;
[m, n] = size(st);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

startNuc=1;
endNuc=n;


tic;

tic;


point=zeros(endNuc-startNuc+6,2);   %the result matrix
temp_center=(A+C+G+T)/4;   % the initial center.
point(1,:)=temp_center;
point(2,:)=A;
point(3,:)=C;
point(4,:)=G;
point(5,:)=T;
% for i=startNuc:endNuc;
%      switch st(i)    % determine which base A,T,G,C it is.
%         case 'A'
%             temp_p=A;  %assign the corresponding vetex
%         case 'G'
%             temp_p=G;
%         case 'C'
%             temp_p=C;
%         case 'T'
%             temp_p=T;
%     end
%     point(i+5,:)=(temp_center+temp_p)/2; %find the result. %IF YOU WANT TO
% %    STORE EVERYTHING ---ALSO UNCOMMENT 89 !!!!!!!!
%      pointTemp=(temp_center+temp_p)/2; % if you DONT want to store
% %    temp_center=point(i,:);  %define the new center. %IF YOU WANT TO STORE EVERYTHING
%     temp_center=pointTemp;  %define the new center.    % if you DONT want to store
%     %plot(temp_center(1),temp_center(2));
% end
for i=startNuc:endNuc;
     switch st(i)    % determine which base A,T,G,C it is.
        case 'a'
            temp_p=A;  %assign the corresponding vetex
        case 'g'
            temp_p=G;
        case 'c'
            temp_p=C;
        case 't'
            temp_p=T;
    end
%     point(i+5,:)=(temp_center+temp_p)/2; %find the result. %IF YOU WANT TO
% %    STORE EVERYTHING ---ALSO UNCOMMENT 89 !!!!!!!!
%      pointTemp=(temp_center+temp_p)/2; % if you DONT want to store
% %    temp_center=point(i,:);  %define the new center. %IF YOU WANT TO STORE EVERYTHING
%     temp_center=pointTemp;  %define the new center.    % if you DONT want to store
%     %plot(temp_center(1),temp_center(2));
        
    point(i+5,:)=(temp_center+temp_p)/2; %find the result. %IF YOU WANT TO
    %STORE EVERYTHING ---ALSO UNCOMMENT 89 !!!!!!!!
    %pointTemp=(temp_center+temp_p)/2; % if you DONT want to store
    temp_center=point(i+5,:);  %define the new center. %IF YOU WANT TO STORE EVERYTHING
    %temp_center=pointTemp;  %define the new center.    % if you DONT want to store
    %plot(temp_center(1),temp_center(2));
end


axis off

hold on



% text(-0.25,0,'A');
% text(-0.25,10,'C');
% text(10,10.05,'G');
% text(10.05,0,'T');

%%%
scatter(point(:,1),point(:,2),1,'black')
line(A,C,'Color','k'); % This makes the box surrounding the plot
line(C,G,'Color','k');
line(G,T,'Color','k');
line(T,A,'Color','k');
%figure('Visible','off');
 %hold on
 %for i=1:length(point)
  %   plot(point(i,1),point(i,2),'k-');
 %end
% hf=figure;
%saveas(gcf,'(filename)','jpg'); 
%print('-djpeg', filename(1:(end)-6), '-r 150');

% pos = get(gcf,'Position'); 
% set(gcf, 'color', 'white');
% I= getframe(gcf,[0 0 pos(3) pos(4)]);

I= getframe; % same as above three lines but without the massive border

plot_output = I.cdata;

set(gca,'Position',[0 0 1 1])
hold off;

 %print(image,'-dpng',strcat(filename(1:(end)-6),'(',num2str(startNuc),'-',num2str(endNuc),')','.png'));
% crop(strcat(filename(1:(end)-6),'(',num2str(startNuc),'-',num2str(endNuc),')','.png'),0)

%figure('Visible','on');
clear('image')


toc;
%hold off

% close all
% clear all
% clc;                

%end
%% This is working, how come above is not ??
% image=figure('Visible','off');
% hold on
% for i=1:5    
%     plot(i,i+1,'*');
% end
% print(image,'-dpng','delete.png');
