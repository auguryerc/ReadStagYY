clear all
close all

%istep  time  F_top  F_bot   Tmin  Tmean Tmax  Vmin  Vrms Vmax  eta_min  eta_mean  eta_max   ra_eff Nu_top  Nu_bot  
% C_min  C_mean  C_max F_mean  F_max  erupt_rate  erupta erupt_heatflux entrainment Cmass_error H_int 
%---------
blue1 =[0.1 0.4 1];
blue2 =[0.1 0.7 1];
blue3 =[0.1 0.9 1];
cgreen = [0.1 0.7 0.1];
green = [0 0.9 0.2];
yellowy = [1 1 0];
orange_light = [1 0.7 0];
cred = [1 0 0];
cblue = [0 0 1];
orange = [1 0.5 0];
purple = [0.7 0 1];
purple2 = [1 0.2 1];
blacky = [0 0 0];
col=1;
c(col,:)= cblue; col=col+1;
c(col,:)= blue2; col=col+1;
c(col,:)= blue3; col=col+1;
c(col,:)= cred; col=col+1;
c(col,:)= cgreen; col=col+1;
%c(col,:)= yellowy; col=col+1;
c(col,:)= orange_light; col=col+1;
%c(col,:)= orange; col=col+1;
c(col,:)= cred; col=col+1;
c(col,:)= purple2; col=col+1;
c(col,:)= purple; col=col+1;
c(col,:)= blacky;col=col+1;
clear col

%----------------------------------
netalaw   = 4;         %number of viscosity laws used
directory='/Users/coltice/Desktop/Work/COURS/newM2/2015/export_calcul/';%'/Users/coltice/Desktop/Work/Projets/Lea/';%'/Users/Nicolas/Work/Research/PROJETS/Plates/PJB2/Temp/';%
cd(directory)
%sed2='grep -l "istep" *_time.dat | xargs -n 1 sed -i.2 "1d" ;';
%remove='rm -f *.dat.2';
%system(sed2);
%system(remove);

vsmean=1100;
ttrans=2900000/39560.;
ti=0;
        
figure(1)

for etalaw=4:netalaw
      
  if (etalaw==1) 
        titi='YS= 0.5e4', 
        refname ='PJB1_YS05_Rh32_Ts06';
    end
    
   if (etalaw==2)
        titi='YS= 1.e4', 
        refname ='PJB1_YS1_Rh32_Ts06';
   end
  
   if (etalaw==3)
        titi='YS= 1.5e4', 
        refname ='PJB1_YS15_Rh32_Ts06';
   end

      if (etalaw==4)
        titi='PJB1_E0', 
        refname ='2D_jump_martina';
      end
      
   if (etalaw==5)
        titi='PJB1_E9', 
        refname ='Ra6-YS2';
   end
   
     if (etalaw==6) 
        titi='PJB1_E13', 
        refname ='E8_YS1.5';
    end
    
   if (etalaw==7)
        titi='PJB1_E18', 
        refname ='PJB1_E18';
   end
  
   if (etalaw==8)
        titi='PJB1_E15', 
        refname ='PJB1_E15';
   end

      if (etalaw==9)
        titi='YS= 4.5e4', 
        refname ='PJB1_YS45_Rh32_Ts06';
      end
      
   if (etalaw==10)
        titi='YS= 5.e4', 
        refname ='PJB1_YS5_Rh32_Ts06';
   end
   
titlename(etalaw)=cellstr(titi);
   
fname=eval([sprintf('refname')]); 
dataref=load(sprintf('%s%s_time.dat',directory,fname));

    for nst=1:size(dataref,1)
        nstep(nst)=dataref(nst,1)-dataref(1,1);
        top_hf(nst,etalaw)=dataref(nst,3);
        bot_hf(nst,etalaw)=dataref(nst,4);
        THF(nst)=dataref(nst,3);
        BHF(nst)=dataref(nst,4);
        Tmean(nst) = dataref(nst,6);
        vrms(nst)=dataref(nst,9);
        total_time_ref(nst)=-(ti-dataref(nst,2).*ttrans*vsmean);
        time(nst,etalaw)=total_time_ref(nst);  
    end
   
subplot(2,2,1)
plot(total_time_ref,vrms,'color',c(etalaw,:),'LineWidth',1.4)
set(gca,'FontSize',16,'fontweight','bold')
xlabel('Time (Ma)');
ylabel('Vmean (non dim)')
title('velocity')
hold on

subplot(2,2,2)
plot(total_time_ref,THF,'color',c(etalaw,:),'LineWidth',1.4)
set(gca,'FontSize',16,'fontweight','bold')
xlabel('Time (Ma)')
ylabel('mean top HF (non dim)')
title('Surface heat Flux')
hold on

subplot(2,2,3)
plot(total_time_ref,Tmean,'color',c(etalaw,:),'LineWidth',1.4)
set(gca,'FontSize',16,'fontweight','bold')
ylim([0 1])
xlabel('Time (Ma)')
ylabel('Tmean (non dim)')
title('Temperature')
hold on

subplot(2,2,4)
plot(total_time_ref,BHF,'color',c(etalaw,:),'LineWidth',1.4) 
set(gca,'FontSize',16,'fontweight','bold')
xlabel('Time (Ma)')
ylabel( 'mean bot HF (non dim)')
title('Bottom heat flux')
hold on

clear nstep THF BHF Tmean vrms total_time_ref titi

end
legend(titlename)

figure(2)
ratio=bot_hf./top_hf.*100;
for etalaw=1:netalaw
plot(time(:,etalaw),ratio(:,etalaw),'color',c(etalaw,:),'LineWidth',1.4)
hold on
end 
set(gca,'FontSize',16,'fontweight','bold')
ylim([0 100])
xlabel('Time (Ma)')
title( 'Ratio bottom/top heat flux (%)')
legend(titlename)
