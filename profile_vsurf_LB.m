clear all
close all

% Program to plot depth-profiles read in from the files rprof.dat.
% For example horizantally average temperature vs. depth
% There is a loop over the timesteps which creates the time evolution.

% ------Parameters------------------------------------------------------
slowdown=0.005;     % time delay between plotting of frames for easier viewing

netalaw=5;
directory='/Users/coltice/Desktop/Work/Projets/Predict/Today/';%'/Users/Nicolas/Work/Research/PROJETS/Plates/PJB2/';%'/Users/Nicolas/Work/Research/PROJETS/Plates/PJB2';%'/Users/Nicolas/Work/Research/PROJETS/MARIE/Longrun/';%'/Users/Nicolas/Work/Research/PROJETS/Plates/PJB2/';%'/Users/Nicolas/Work/Research/PROJETS/FASTOPT/StagAge';%'/Users/coltice/Desktop/Work/Projets/Lea';
number_start=1;
number_end=1;
imax=96 ; 

option_video=0;

%--------Edit files _rprof.dat before first -------------------------------------------
cd(directory);
sed1='sed -i".2" -e "s/0-/0E-/g" -e "s/0+/0E+/g" -e "s/1-/1E-/g" -e "s/1+/1E+/g" -e "s/2-/2E-/g" -e "s/2+/2E+/g" -e "s/3-/3E-/g" -e "s/3+/3E+/g" -e "s/4-/4E-/g" -e "s/4+/4E+/g" -e "s/5-/5E-/g" -e "s/5+/5E+/g" -e "s/6-/6E-/g" -e "s/6+/6E+/g" -e "s/7-/7E-/g" -e "s/7+/7E+/g" -e "s/8-/8E-/g" -e "s/8+/8E+/g" -e "s/9-/9E-/g" -e "s/9+/9E+/g" *_rprof.dat ;';
sed2='grep -l "istep" *_time.dat | xargs -n 1 sed -i.2 "1d" ;';
remove='rm -f *.dat.2';
system(sed1);
system(sed2);
system(remove);

% --------calcul solidus-------------------------------------------------
deltaT_dimensional  = 1600.0;
g_dimensional       = 9.81;
D_dimensional       = 2890.e3;
D_dim_top           = 0.0; %default D_dim_top= 0.0
    %if (ddim>660)
    %    simple_solidus = (2760. + 0.45*ddim + 1700*(erfnr(ddim/1000.0)-1)) /...
    %    deltaT_dimensional;
    %else
   % simple_solidus = (2050. + 0.62*ddim +  660*(erfnr(ddim/220.0)-1)) /...
    %         deltaT_dimensional;
    %end

% ------------------------------------------------------------------------
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

% ---------------------------------------------------------------------- 
it=1;
for etalaw=5:netalaw
    
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
        refname ='PJB1_E0';
      end
      
   if (etalaw==5)
        titi='Inst-today', 
        refname ='Inst-today'
   end
   
     if (etalaw==6) 
        titi='PJB1_E13', 
        refname ='PJB_wcont_YS2';
    end
    
   if (etalaw==7)
        titi='PJB1_E18', 
        refname ='PJB_wcont_YS2-E3';
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
        
if (option_video==1)
    writerObj = VideoWriter(sprintf('%sProfiles%i.avi',directory,etalaw));
    writerObj.FrameRate = 4;
    open(writerObj);
end

meltT(1)=0.65;   % Melting line for viscosity reduction
meltT(2)=1.5;   % This one- T=0.6 + 2*(1-z)
meltz(1)=1.0;
meltz(2)=0.75;

cd(directory);
fid=fopen(sprintf('%s_rprof.dat',refname));


for j=number_start:number_end  % timestep-loop

fscanf(fid,'%s',1);         % ******step:

timestep(j)=fscanf(fid,'%d',1); % 0
fscanf(fid,'%s',1);         % ;
fscanf(fid,'%s',1);         % time
fscanf(fid,'%s',1);         % =
simtime=fscanf(fid,'%e',1);  % 0   

for i=1:imax
    
z(i)=fscanf(fid,'%e',1);  % z-coordinate
tave(j,i)=fscanf(fid,'%e',1);% average Temperature     
tmin(j,i)=fscanf(fid,'%e',1);% minimum temp 
tmax(j,i)=fscanf(fid,'%e',1);% maximum temp 
vave(j,i)=fscanf(fid,'%e',1);% v 
vmin(j,i)=fscanf(fid,'%e',1);
vmax(j,i)=fscanf(fid,'%e',1);
fscanf(fid,'%e',1);        % vz 
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);        % vhor 
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
etaave(j,i)=fscanf(fid,'%e',1);   % Viscosity 
etamin(j,i)=fscanf(fid,'%e',1);
etamax(j,i)=fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);
fscanf(fid,'%e',1);

end  % ----------------------------------- end profile-loop
vsurf(j)=vave(j,imax);

vave(j,:);


if (option_video==1) 
    figure(1)
    subplot(1,2,1)
    plot(vave(j,:),z,'r','LineWidth',1.2)%,meltT,meltz,'g')  % temperature
%    plot(tmin(j,:),z,'b',tave(j,:),z,'k',tmax(j,:),z,'r','LineWidth',1.2)%,meltT,meltz,'g')  % temperature
    legend('min','avg','max','Location','East')
    ti=sprintf('%s : Temperature profile',titlename(etalaw,:));
    title(ti,'FontSize',14,'fontweight','bold')

    subplot(1,2,2)
    semilogx(etaave(j,:),z,'k',etamin(j,:),z,'b--',etamax(j,:),z,'r--')
    xlim([0.01 10000])
    legend('average','min','max','Location','east')
    ti2=sprintf('and Viscosity profile at %i Ma',j);
    title(ti2,'FontSize',14,'fontweight','bold')

    picname=sprintf('Profiles%iMa',j);
    print('-dpng',picname)
    Im = imread(sprintf('Profiles%iMa.png',j),'BackgroundColor','none');
    writeVideo(writerObj,Im);
end 

end  %-------------------------------------- end time-loop

if (option_video==1)
    close(writerObj);
    remove=strcat({'rm '},'Profiles**Ma.png');
    system(remove{:});
end

figure(1)

subplot(1,2,1)

%if (etalaw==5)
    d=transpose(1-z);
    ddim =(d*D_dimensional+D_dim_top) * (g_dimensional/9.81) /1000. ; %in km 
    linear_solidus2=0.7+7*d;
    linear_solidus1=0.6+2*d;
    linear_solidus=0.6+7.5*d;
%end    

etadim=etaave.*3.1e22;
tdim=tave.*1400.0 ;
d
tave(j,:)'
etaave(j,:)'
z(:);
etadim=etaave.*3.1e22;
tdim=tave.*1400.0 ;

plot(tave(j,:),-ddim,'color',c(etalaw,:),'LineWidth',1.2); hold on
titlename(it)=cellstr(titi); it=it+1;
plot(tmax(j,:),-ddim,'--','color',c(etalaw,:),'LineWidth',1.2); hold on
titi2=sprintf('%s extremum',titi);
titlename(it)=cellstr(titi2);it=it+1;

subplot(1,2,2)
semilogx(etaave(j,:),-ddim,etamin(j,:),-ddim,'--','color',c(etalaw,:),'LineWidth',1.2); hold on

vsurfmoy(etalaw)=mean(vsurf);
clear vsurf titi titi2 vmin vave vmax tave tmin tmax  etaave etamin etamax


end  % ---------------------------------------- end etalaw loop

subplot(1,2,1)
plot(linear_solidus,-ddim,'r--','LineWidth',1.2); hold on 
titlename(it)=cellstr('Linear solidus T=0.6+7.5*d ');
titlename(it+1)=cellstr('Linear solidus Tsol2');
titlename(it+2)=cellstr('Linear solidus Paul');

ylim([-3000 0])
xlim([0 1.3])
plot([1 1],[-3000 0],'k-')
set(gca,'xtick',0:0.1:1.2,'FontSize',14)
set(gca,'ytick',-3000:200:0,'FontSize',14)
ylabel('Depth (km)','FontSize',16)
title('Temperature profile','FontSize',16,'fontweight','bold')
legend(titlename,'location','best')

subplot(1,2,2)
semilogx([1 1],[-3000 0],'k-')
ylim([-3000 0])
xlim([0.001 1000])
set(gca,'ytick',-3000:200:0,'FontSize',14)
title('Viscosity profile','FontSize',16,'fontweight','bold')


vsurfmoy*number_end/(number_end-number_start);

%what is in rprof.dat:

%write(1,'(19(1pe12.4),$)') zg(1,iz),(tprof(j,iz),j=1,3),(vprof(j,iz),j=1,3),(vzprof(j,iz),j=1,3),&
%                      (vhprof(j,iz),j=1,3),(etaprof(j,iz),j=1,3),(eprof(j,iz),j=1,3)
%write(1,'(29(1pe12.4))') (sprof(j,iz),j=1,3),((wprof(j,iz,k),j=1,3),k=1,2), &
%                      (dprof(j,iz),j=1,3),(enprof(j,iz),j=1,5),(cprof(j,iz),j=1,3),(denprof(j,iz),j=1,3), &
%                      (airprof(j,iz),j=1,3),(primprof(j,iz),j=1,3)

