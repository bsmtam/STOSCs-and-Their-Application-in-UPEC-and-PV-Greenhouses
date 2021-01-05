clear all
close all
clc



[FileName,PathName,FilterIndex] = uigetfile('*.csv','select the file');
Folder=PathName;
files=dir(fullfile(Folder,'*.csv'));
currentend = regexpi(FileName, '800');   % insesitive of cap
hyphenindex = regexpi(FileName, '-');
hyi = find(hyphenindex < currentend(end)) ;
currenthead = hyphenindex (hyi(end));
mmm= [1:1:currenthead-4];
%mmm= [hyphenindex(2):1:currenthead-4];
name = FileName(mmm) ;
mkdir(Folder,'results');
erg={zeros(1),zeros(1),zeros(1),zeros(1)};
i=0;
%   EpeakFileHandle=fopen([Folder,'results\',name,'.txt'],'a');
%         fprintf(EpeakFileHandle,'\t%1$s\t%2$s\n','Ev', 'QE');
%         fclose(EpeakFileHandle)
for k = 1:length(files)
    
    data=csvread(files(k).name,9,0);
    if (regexpi(files(k).name,name))
        i=i+1;
        lambda=data(:,1);
        E=1240./lambda;
        
        erg1 = [1240./lambda, data(:,3)];
        erg{i}=erg1;
        
    end
    
    
end

B=vertcat(erg{1},erg{2},erg{3});
C=sortrows(B);
%C(:,1)=smooth(C(:,1),10);
%C(:,2)=smooth(C(:,2),10);
save([Folder,'results\',name,'3','.csv'] ,'C' , '-ASCII');
% axis label

        data=load([Folder,'results\',name,'3.csv']);
       B=[data(:,1), data(:,2)];
C=sortrows(B);
     plot(C(:,1),C(:,2));

        ylabel('EQE \phi') % label left y-axis

        xlabel('energy E [eV]');

        title(name);

        saveas(gcf,[Folder,'results\',name,'.jpg'],'jpg');
