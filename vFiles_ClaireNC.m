clear all
close all
             
directory='/Users/nicolas/Work/Research/PROJETS/Plates/PJB2/';%'/Users/coltice/Desktop/Work/Projets/Lea/';%'/Users/Nicolas/Work/Research/PROJETS/Plates/PJB2/Temp/Claire-YS1.5/';%
refname='Rectibet';%'try10-nowc';%'try40-YS2';%'Inst00-256';%
vsmean=1100;
number_start=8;
number_end=8;
increment=1;
nr_save         = [2 2 1];
nrstart        = 1;
size_reduction = true;
nz=96
nzd=70
vscale=3.956e-2/vsmean ; 

for step_number=number_start:increment:number_end
        
  step_number
%[~,~,~,~,~,~,~, rcmb] = ReadStag3D_Lea(directory,refname,step_number,'viscosity');
[~, ~, ~, ~,~, ~, ~, rcmb] = ReadStag3D2_reduce_size_LB(directory, refname, step_number,  'viscosity',size_reduction,nr_save,nrstart)

%[~, X_3D,Y_3D,Z_3D, VX_3D_1,VY_3D_1,VZ_3D_1,~, VX_3D_2,VY_3D_2,VZ_3D_2, ~,~] = ReadStag3D_Lea(directory,refname,step_number,'velocity');  
%ReadStag3D_Lea(directory,refname,step_number,'velocity')  
[nnb, X_3D, Y_3D, Z_3D,VX_3D_1, VY_3D_1, VZ_3D_1, P_3D_1, VX_3D_2, VY_3D_2, VZ_3D_2, P_3D_2, time] = ReadStag3D2_reduce_size_LB(directory, refname, step_number,   'velocity',size_reduction,nr_save,nrstart);
%[~, X_3D, Y_3D, Z_3D,VX_3D_1, VY_3D_1, VZ_3D_1, ~, VX_3D_2, VY_3D_2, VZ_3D_2, ~, ~] = ReadStag3D2_reduce_size_LB(directory, refname, step_number,   'velocity',size_reduction,nr_save,nrstart);
         
% Transform coordinates for Yin & Yang grids
    R  = Z_3D+rcmb;
    lat = pi/4-X_3D;%+pi/2;
    lon = Y_3D-3*pi/4;           
% Yin grid
    X_3D_1 = R.*cos(lat).*cos(lon);   
    Y_3D_1 = R.*cos(lat).*sin(lon);
    Z_3D_1 = R.*sin(lat);             
%  Yang grid            
    X_3D_2 = -X_3D_1;       
    Y_3D_2 =  Z_3D_1;
    Z_3D_2 =  Y_3D_1;
    
    size(X_3D_2,1)
    size(X_3D_2,2)
    size(X_3D_2,3)


                
% Transform velocities
% on Yin grid
    Vtheta = VX_3D_1; 
    Vphi = VY_3D_1; 
    Vr = VZ_3D_1;                                 
    VX_3D_1 =  Vtheta.*sin(lat).*cos(lon) - Vphi.*sin(lon) + Vr.*cos(lat).*cos(lon)*vscale;
    VY_3D_1 =  Vtheta.*sin(lat).*sin(lon) + Vphi.*cos(lon) + Vr.*cos(lat).*sin(lon)*vscale;
    VZ_3D_1 = -Vtheta.*cos(lat) + Vr.*sin(lat)*vscale;
    Vr1 = VZ_3D_1; 
% on Yang grid
    Vtheta = VX_3D_2; 
    Vphi = VY_3D_2; 
    Vr = VZ_3D_2;                                 
    VX_3D_2 = -(Vtheta.*sin(lat).*cos(lon) - Vphi.*sin(lon) + Vr.*cos(lat).*cos(lon));
    VZ_3D_2 =  Vtheta.*sin(lat).*sin(lon) + Vphi.*cos(lon) + Vr.*cos(lat).*sin(lon);    
    VY_3D_2 = -Vtheta.*cos(lat) + Vr.*sin(lat);
    Vr2 = VZ_3D_2;    

%conversion to spherical coordinates (from Fabio's script)
lat1    = atan2(sqrt(X_3D_1.^2 + Y_3D_1.^2),Z_3D_1);  % cos-1 (z/r)
long1   = atan2(Y_3D_1,X_3D_1);  % tan-1 (y/x)
lat2    = atan2(sqrt(X_3D_2.^2 + Y_3D_2.^2),Z_3D_2);  % cos-1 (z/r)
long2   = atan2(Y_3D_2,X_3D_2);  % tan-1 (y/x)


%[Az_1,Elev_1,vecnorm_1] = cart2sph(VX_3D_1,VY_3D_1,VZ_3D_1);
%[Az_2,Elev_2,vecnorm_2] = cart2sph(VX_3D_2,VY_3D_2,VZ_3D_2);

 Vlat1   = VX_3D_1.*(cos(long1).*cos(lat1)) + VY_3D_1.*(sin(long1).*cos(lat1)) - VZ_3D_1.*(sin(lat1));
 Vlong1	= -VX_3D_1.*(sin(long1)) + VY_3D_1.*(cos(long1));
 Vlat2   = VX_3D_2.*(cos(long2).*cos(lat2)) + VY_3D_2.*(sin(long2).*cos(lat2)) - VZ_3D_2.*(sin(lat2));
 Vlong2  = -VX_3D_2.*(sin(long2)) + VY_3D_2.*(cos(long2));

%VXYY(:,:,1) = Az_1(:,:,nz)/pi*180.;
%VXYY(:,:,2) = Az_2(:,:,nz)/pi*180.;

 VXYY(:,:,1) = -Vlat1(:,:,nzd); % Yin  grid data
 VXYY(:,:,2) = -Vlat2(:,:,nzd); % Yang grid data        
 [VX2D,~,~]  = YYtoMap2(VXYY);

 Vrr(:,:,1) = Vr1(:,:,nzd)
 Vrr(:,:,2) = Vr2(:,:,nzd)
 [Vr2D,~,~]  = YYtoMap2(Vrr);

%VYYY(:,:,1) = vecnorm_1(:,:,nz);
%VYYY(:,:,2) = vecnorm_2(:,:,nz);
 VYYY(:,:,1) = Vlong1(:,:,nzd); % Yin  grid data
 VYYY(:,:,2) = Vlong2(:,:,nzd); % Yang grid data        
[VY2D,theta,phi]  = YYtoMap2(VYYY); 
    
% Ecriture fichier Claire avec pas tous les points   
  
  X=[X_3D_1 X_3D_2];
  Y=[Y_3D_1 Y_3D_2];
  Z=[Z_3D_1 Z_3D_2];
  
  nxmax=size(X,1)
  nymax=size(X,2)
  nzmax=size(X,3)
 
  tot=(size(X,1)/4*size(X,2)/4);

  %Vitesse ? transformer en m/a : vtoday_Earth(m/a) / vsmean_model(nondim) 

  thetad             = 90.-theta*180./3.14159;
  phid               = phi*180./3.14159;
  [DAT]    = savedatavec(VX2D,VY2D,thetad,phid);
  itt=num2str(step_number);
  save_str = strcat(directory,refname,'_t',itt,'_projectV.dat');
  save(save_str,'DAT','-ascii'); 
  
  [DATx]    = savedata(VX2D,thetad,phid);
  itt=num2str(step_number);
  save_str = strcat(directory,refname,'_t',itt,'_projectTv.dat');
  save(save_str,'DATx','-ascii'); 

  [DATy]    = savedata(VY2D,thetad,phid);
  itt=num2str(step_number);
  save_str = strcat(directory,refname,'_t',itt,'_projectRv.dat');
  save(save_str,'DATy','-ascii'); 
  
  [DATz]    = savedata(Vr2D,thetad,phid);
  itt=num2str(step_number);
  save_str = strcat(directory,refname,'_t',itt,'_projectZv.dat');
  save(save_str,'DATz','-ascii'); 
  
  
%     k=1; 
%   while k<=tot 
%     for i=1:4:nxmax
%         for j=1:4:nymax
%             Mforsave(k,1)=  X(i,j,nzmax) ;
%             Mforsave(k,2)=  Y(i,j,nzmax) ;
%             Mforsave(k,3)=  Z(i,j,nzmax) ;
%             Mforsave(k,4)= VX2D(i,j) * vscale;
%             Mforsave(k,5)= VY2D(i,j) * vscale; 
%             k=k+1;
%         end
%     end
%  end
% k
%                 
%     currentFileR = sprintf('%s%s%iXYZ.xy',directory,refname,step_number);
%     save(currentFileR,'Mforsave','-ascii');
%     clear Mforsave 
%     clear k
end

