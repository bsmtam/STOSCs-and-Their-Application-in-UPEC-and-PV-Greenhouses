clear all
close all
clc

[FileName,PathName,FilterIndex] = uigetfile('*.csv','select the file')
Folder=PathName;
files=dir(fullfile(Folder,'*.csv'));
currentend = regexpi(FileName, '800');   % insesitive of cap
hyphenindex = regexpi(FileName, '-');
hyi = find(hyphenindex < currentend(end)) ;
currenthead = hyphenindex (hyi(end));
mmm= [hyphenindex(2):1:currenthead-4];
name = FileName(mmm) ;  
mkdir(Folder,'results');
%   EpeakFileHandle=fopen([Folder,'results\',name,'.txt'],'a');
%         fprintf(EpeakFileHandle,'\t%1$s\t%2$s\n','Ev', 'QE');
%         fclose(EpeakFileHandle)
for k = 1:length(files)
     
    data=csvread(files(k).name,9,0);
     if (regexpi(files(k).name,name))
         
          lambda=data(:,1);
          E=1240./lambda;
          
        erg1 = [1240./lambda, data(:,3)];
        
        save([Folder,'results\',name,'.csv'] ,'erg1' , '-ASCII','-append');
     end
end 

      
        
        % axis label
        
        data=load([Folder,'results\',name,'.csv']);
       erg1=[data(:,1), data(:,2)];
        B=sort(erg1);
     semilogy(B(:,1),B(:,2));
        
        ylabel('EQE \phi') % label left y-axis
        
        xlabel('energy E [eV]');
        
        title(name);
        
        saveas(gcf,[Folder,'results\',name,'.jpg'],'jpg');
 