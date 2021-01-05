Dvnrmin=interp1(resmin(:,1),resmin(:,14),0.8:0.1:2);
Dvnrmax=interp1(resmax(:,1),resmax(:,14),0.8:0.1:2);

ii=0;
for EG=0.8:0.1:2
    ii=1+ii;
[data,Y]=ShockleyQueisser(EG+0.2);
Jsclimit(ii)=data.Jsc*0.8;
FFlimit(ii)=data.FF;
Efflimitmin(ii)=(data.Voc-Dvnrmin(ii))*data.Jsc*0.8*0.8;
Efflimitmax(ii)=(data.Voc)*data.Jsc*data.FF;

end 
figure(1)
hold on
plot(0.8:0.1:2,Efflimitmin,'b-*',0.8:0.1:2,Efflimitmax,'r-*')
figure(2)
plot(0.8:0.1:2,FFlimit,'-*')
figure(3)
plot(0.8:0.1:2,Jsclimit,'-*')