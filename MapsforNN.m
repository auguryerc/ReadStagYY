% Routine for calculation of surface heatflow, seafloor age and ocean
% bathymetry. This routine calculates the surface heat flow via Fourier's
% law from the binary temperature output of the mantle convection code
% STAG-YY (Tackley,PEPI,2008) in YinYang geometry. It can do the
% partitioning between oceanic and continental regions as well as the
% differentiation between different continents.The heat flow distribution 
% is converted into adistribution of seafloor ages by using the half-space 
% cooling approximation. The age distribution can then be used to estimate 
% the bathymetry of any point in the ocean and with that the total volume
% of ocean basins can be track in time. Additional features are the
% production rate of new seafloor, the mean age of the lithosphere and
% mollweid projections of surface temperature, composition, heatflow, age,
% bathymetry, viscosity and velocity.

% Major rewriting of old surface_age_bathymetry script, aiming for easier
% handling and possible distribution to others

% written by Tobias Rolf, ETH Zurich, October 2012
% ...tested in YinYang geometry, 2D spherical annullus is implemented and
% ...heat flow calculation seems to work fine. However the age distribution
% ...can look quite strange. It can be kind of flat of triangular, but it
% ...also often has some contribution at larger age. This might be because
% ...the number of age bins is very small in 2D (512 usually). So, the very
% ...old lithosphere, which usually sticks to the continent-ocean boundary
% ...has a much larger influence on the age distribution than in 3D.
% ...Bathymetry looks even worse, but is anyway based on the age distribution.

% ...24/10/2012 - before the volumetrically-averaged rms-velocity was used
% ...to estimate the model transit time. Now changed to using the
% ...surface-averaged rms-velocity! Should give larger age, smaller C0!

% ...22/10/2013 - added divergence and vorticity projection switches

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; tinit = cputime; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DECIDE WHICH DATA YOU WANT TO PROCESS
dir          = '/Volumes/HD1/AI-Plates/';%'/Users/Nicolas/Work/Research/PROJETS/Plates/PJB2/Temp/Claire-YS1.5/';%'/Users/Nicolas/Work/Research/PROJETS/Lea/Plates/PJB2/StegmanEurasia/'; 
fname        = 'today-2';%'try10-nowc';%'try40-YS2';%'Inst10-YS1';%'Inst00';%                 % File name stem 
file_stem    = strcat(dir,fname);
GridType     = 'YinYang';                   % Grid geometry, either YinYang or Annullus
ncont        = 0;                           % Number of continents
t0           = 1.14976E-01;                   % Point of time origin, transition from fixed to mobile continentd
number_trans = 1;                         % Frame number that belongs to time t0
transit_M    = 1/1100.;                    % Model Transit time - UNTIL RECENTLY THIS CONSIDERED the V-averaged rms_velocity, but it should consider surface_rms
Tscale       = 1.000;                       % Non-dim scaling factor for temperature dimensionalization (mean subocean-T of ROSA case: 6x5_ys1e4_int20)
Eta_Ref      = 1.0e+22;                     % Reference viscosity in Pa s
nstart       = 645;                         % first output frame to process
nend         = 1810;                         % last output frame to process
increment    = 1;                           % frame increment, set to 1 if you want to process all frames
n            = (nend-nstart)/increment + 1; % total number of frames to process -- MUST BE AN INTEGER!
Rsurf        = 2.2;   % Non-dim surface radius
ndepth       = 64     % Number of the layer to map

% LOAD SUBOCEANIC TEMPERATURE OF THIS CASE 
%data_str   = strcat(dir,fname,'_subtemp.dat');
%data       = load(data_str);
%subTo      = data(:,4);
subTo       = ones(nend,1)*.86;     %%% FOR TESTING ONLY %%%

% LOAD SURFACE VELOCITIES OF THIS CASE (IN NON-DIM UNITS) TO ESTIMATE THE
% TRANSIT TIME OF THIS CASE
data        = ones(nend,2);     %%% FOR TESTING ONLY %%%
%data_str    = strcat(dir,fname,'_nondim_velocities.dat');
%data        = load(data_str);
% using global average (as before)
vrms_g     = data(:,2);
vmean      = mean(vrms_g); 
%transit_M  = 1./vmean;
% using oceanic average only
%vrms_o     = data(:,3);
%vmean      = mean(vrms_o); 
%transit_M = 1./vmean;

% SPECIFY PHYSICAL PROPERTIES
k0         = 3.15;                      % thermal conductivity in W/m/K
deltaT     = 1300;                      % temperature drop across lithosphere in K
Tsurf      = 273;                       % surface temperature in K
kappa0     = 1.0e-06;                   % thermal diffusivity in m/s^2
alpha0     = 3.0e-05;                   % thermal expansivity in 1/K
rho_m      = 3300;                      % Density of mantle rock in kg/m^3
rho_w      = 1000;                      % Density of ocean water in kg/m^3
drho       = rho_m-rho_w;               % Density contrast
ocean_area = 3.57e14;                   % Present-day oceanic area in m^2
Myr        = 1./(365*24*60*60*1.0e6);   % conversion from s to Ma 
transit_E  =  3.6e6*3.141e7; %2.68056e+15;               % Earth's transit time in seconds
Radius     = 6.371e6;                   % Earth's surface radius in m
Dscale     = 2.890e6;                   % Earth's mantle thickness in m

% OTHER PARAMETERS SPECIFIC FOR CALCULATION

Cthresh             = 0.50;    % >continent, <ocean, default=0.5
WriteComposition    = false;    % continental roots or not, default=true
WriteBelts          = false;   % mobile belts or sutures, default=false, NOT YET TESTED 
WriteContinent      = false;    % multiple continents (nrc field?), default=true
WriteAge            = false;
WriteEdot           = false;
Age_thresh          = 0.05;    % non-dim ages larger than that are neglected,default=0.05
incr                = 0.00025; % non-dim size of age bins,default=0.00025
incr_dim            = 2.5e+14/21; % dim size of age bins in seconds,default=2.5e14

% PARAMETER SPECIFIC FOR AGE/BATHYMETRY CALCULATION

min_age             = 0.0;    % non-dim minimum age
max_age             = 0.1;    % non_dim maximum age
min_age_dim         = 0.0;    % dimensional minimum age in seconds
max_age_dim         = 1.0e17; % dimensional maximum age in seconds
b_incr              = 100;    % non-dim size of bathymetry bin,default=100
bath_min            = 0.0;    % minimum bathymetry,default=0.0
bath_max            = 80000.0; % maximum bathymetry,default=8000

% DECIDE WHAT YOU WANT TO PLOT/SAVE
dimensional_units   = false;
plotany             = false;
projectT            = false;     save_projectT_data    = false;     % project temperature field and save it ?
projectC            = false;      save_projectC_data    = false;     % project composition field and save it ?
projectQ            = false;     save_projectQ_data    = false;     % project heatflow field and save it ?
projectAGE          = false;      save_projectAGE_data  = false;     % project seafloor age field and save it ?
projectBATH         = false;     save_projectBATH_data = false;     % project bathymetry field and save it ?
projectETA          = false;     save_projectETA_data  = false;     % project viscosity field and save it ?
projectDIV          = true;      save_projectDIV_data  = true;     % project divergence field and save it ?
projectVOR          = false;     save_projectVOR_data  = false;     % project vorticity field and save it ?
projectVEL          = false;     % add surface velocity arrows to the switched-on field projections 
projectEDOT         = false;     save_projectEdot_data  = false;

save_age_hist_file  = false;      save_age_hist_plot  = false;      % save age histogram as figure and/or data files
save_bath_hist_file = false;     save_bath_hist_plot  = false;      % save bathymetry histogram and/or data files
save_C0_CX          = false;
save_hf_series      = false;

% SPECIFY PARAMETERS NEEDED FOR PLOTTING

modulo              = 20;             % keep every 'modulo' arrow in quiver plots, neglect the others
arrowscale          = 4.0;            % Scaling factor for the arrows
ncontour            = 20;             % number of contour lines
FontSize            = 18;             % FontSize in all plots
FigSize             = [1 200 1280 640]; % Dimension of rhe projection figures

% INITIALIZE ARRAYS/FIELDS

% ... time series
t               = zeros(n,1);          % Time 
qgt             = zeros(n,1);          % Global surface heat flow
qot             = zeros(n,1);          % Oceanic surface heat flow
qct             = zeros(ncont+1,n);    % Continental surface heat flow
C0              = zeros(n,1);          % Seafloor production rate from linear fit
CX              = zeros(n,1);          % Actual seafloor production rate (no fit)
mean_age        = zeros(n,1);          % Mean age of the seafloor
real_area       = zeros(n,1);          % Fraction of used seafloor to get mean age
ocean_basin     = zeros(n,1);          % Total volume of ocean basins
Tref            = zeros(n,1);          % Reference (suboceanic) temperature
F               = zeros(n,1);          % Time series of figure frame (needed for movies)

% ... non-dim age histograms

xhist           = min_age:incr:max_age+incr;     % non-dim age bin array
stack_nhist     = zeros(length(xhist),2);        % stacked age bin array
stack_harea     = zeros(length(xhist),2);        % stacked area per unit age array
MH              = zeros(length(xhist),n);        % Something ???

% ... dimensional age histograms

xhist_dim       = min_age_dim:incr_dim:max_age_dim+incr_dim;  
stack_nhist_dim = zeros(length(xhist_dim),2);               
stack_harea_dim = zeros(length(xhist_dim),2);       
MH_dim          = zeros(length(xhist_dim),n);

% ... bathymetry histograms

bhist           = bath_min:b_incr:bath_max+b_incr;  % dimensional bathymetry bin array
na              = zeros(length(xhist_dim),1);       % like stack_nhist_dim for bathymetry
nb              = zeros(length(bhist),1);           % Something ???

% ... some counters
num             = 1;
nstack1         = 0;
nstack2         = 0;

% START PROCESSING - MAIN TIME STEPPING LOOP

for itf=nstart:increment:nend 
    
%     if (itf == 452)
%         projectAGE = true;     save_projectAGE_data = true;
%     else
%         projectAGE = false;    save_projectAGE_data = false;
%     end
    
    itf
    nfn = num2str(itf);
    if (itf<10)                    
        nfn = strcat('0000',nfn);
    elseif (itf >= 10 && itf < 100)
        nfn = strcat('000',nfn);
    elseif (itf >= 100 && itf < 1000)
        nfn = strcat('00',nfn);
    elseif (itf >= 1000 && itf < 10000)
        nfn = strcat('0',nfn);
    elseif (itf >= 99999)
        error('Too many input files for unique saving!')
    end
    
    % SWITCH TO APPROPRIATE GEOMETRY
    
    switch GridType
        case 'YinYang'
        
        % FIRST TIME STEP REQUIRES SOME ADDITIONAL INITIALIZATION - SKIP
        % THAT IN FOLLOWING STEPS FOR SPEED
            
        if (itf==nstart)
            [nnb,~,~,Z_3D,T_3D(:,:,:,1),T_3D(:,:,:,2),time,c] = ReadStag3D2_251012(file_stem,itf,'temperature',false ,[1 1 1]);
            
            nx  = size(T_3D,1)
            ny  = size(T_3D,2)
            nz  = size(T_3D,3)     
           
            if     (nx==64  && ny==192);  Grid = load('Cell_area_YY_64x192x32_T');
            elseif (nx==96  && ny==288);  Grid = load('Cell_area_YY_96x288x48_T'); 
            elseif (nx==128 && ny==384);  Grid = load('Cell_area_YY_128_384');
            elseif (nx==192 && ny==576);  Grid = load('Cell_area_YY_192_576');    
           %Nicolas else
           %     error('Cannot get grid data - missing file')
            end
            
                               
            % INITIALIZE COMPOSITIONAL FIELDS
            
            C_3D   = zeros(nx,ny,nz,nnb); % cratonic root 
            C1_3D  = zeros(nx,ny,nz,nnb); % belt or suture                           
            NRC_3D = zeros(nx,ny,nz,nnb); % continent number
            CS_3D  = zeros(nx,ny,nz,nnb); % Strain rate, remplacer nz par 2 pour topo
        else
            [~,~,~,~,T_3D(:,:,:,1),T_3D(:,:,:,2),time,~] = ReadStag3D2_251012(file_stem,itf,'temperature',false ,[1 1 1]);
        end    
        
 %       if WriteComposition;  [~,~,~,~,C_3D(:,:,:,1),C_3D(:,:,:,2),~,~]     = ReadStag3D2_251012(file_stem,itf,'cont root',false ,[1 1 1]);  end
 %       if WriteBelts;        [~,~,~,~,C1_3D(:,:,:,1),C1_3D(:,:,:,2),~,~]   = ReadStag3D2_251012(file_stem,itf,'cont crust',false,[1 1 1]);  end    
        if WriteContinent;    [~,~,~,~,NRC_3D(:,:,:,1),NRC_3D(:,:,:,2),~,~] = ReadStag3D2_251012(file_stem,itf,'continent',false ,[1 1 1]);  end
        if WriteEdot; [~,~,~,~,CS_3D(:,:,:,1),CS_3D(:,:,:,2),~,~] = ReadStag3D2_251012(file_stem,itf,'strain rate',false ,[1 1 1]); end
        if WriteAge;    [~,~,~,~,AGE_3D(:,:,:,1),AGE_3D(:,:,:,2),~,~] = ReadStag3D2_251012(file_stem,itf,'age',false ,[1 1 1]);  end
        
        
        % COMBINE COMPOSITINOAL FIELDS TO ONE FIELD
        sumC(:,:,1,1) = C_3D(:,:,nz,1) + C1_3D(:,:,nz,1);
        sumC(:,:,1,2) = C_3D(:,:,nz,2) + C1_3D(:,:,nz,2);
        
        if (projectETA)
            [~,~,~,~,ETA_3D(:,:,:,1),ETA_3D(:,:,:,2),~,~] = ReadStag3D2_251012(file_stem,itf,'viscosity',false ,[1 1 1]); 
        end 
        
        if (projectVEL==false)
           VX2D = 0.0; % dummy argument during call of function_project3Dfield
           VY2D = 0.0; % dummy argument during call of function_project3Dfield
        end
        
%Nico Topo - remplac par Edot
if (projectEDOT)
            DATAYY(:,:,:)      = CS_3D(:,:,ndepth,:); 
            [DATA2D,theta,phi] = YYtoMap2(DATAYY);
            thetad             = 90.-theta*180./3.14159;
            phid               = phi*180./3.14159;
            magic              = 1;
            cminmax            = [0 1];
            time1              = time;
            log                = false;
          
            if (projectVEL) % add quiver?
                [VX2D,VY2D] = transform_VEL(dir,fname,itf);  
            end 
            
            if (dimensional_units)
                % dimensionallzation using transit time approach
                % v [m/s]  = v' * transit_M*Dscale/transit_E
                % v [cm/a] = v [m/s] * 100 cm/m * 3.1536e7 s/a 
                if (projectVEL)
                    VX2D    = Dscale*transit_M/transit_E*3.1536e9*VX2D;
                    VY2D    = Dscale*transit_M/transit_E*3.1536e9*VY2D;
                end
                time1   = (time-t0)*transit_E/transit_M*3.170979e-17; % in Gyr
                % Composition is always non-dimensional
            end
            
            [h,time_str] = function_project3Dfield(DATA2D,thetad,phid,FontSize,ncontour,time1,itf,magic,FigSize,VX2D,VY2D,arrowscale,modulo,projectVEL,cminmax,dimensional_units,log,1);
            fig_num  = 1000000*magic+itf;
            save_pdf = true; 
            save_eps = false;
            save_jpg = false;
            save_png = false;
            save_str = strcat(dir,fname,'_n',nfn,'_project_Edot');  
% NICO            [ierr]   = savecoolplots(save_str,fig_num,save_pdf,save_eps,save_jpg,save_png);
            close(h)
                     
            if (save_projectEdot_data)
                [DAT] = savedata(DATA2D,thetad,phid);
                save_string = strcat(dir,fname,'_n',nfn,'_project_Edot.dat');
                save(save_string,'DAT','-ascii');        
            end
end       
        
%Nico Topo


        if (projectC)
            DATAYY(:,:,:)      = NRC_3D(:,:,nz,:); 
            [DATA2D,theta,phi] = YYtoMap2(DATAYY);
            thetad             = 90.-theta*180./3.14159;
            phid               = phi*180./3.14159;
            magic              = 1;
            cminmax            = [0 1];
            time1              = time;
            log                = false;
          
            if (projectVEL) % add quiver?
                [VX2D,VY2D] = transform_VEL(dir,fname,itf);  
            end 
            
            if (dimensional_units)
                % dimensionallzation using transit time approach
                % v [m/s]  = v' * transit_M*Dscale/transit_E
                % v [cm/a] = v [m/s] * 100 cm/m * 3.1536e7 s/a 
                if (projectVEL)
                    VX2D    = Dscale*transit_M/transit_E*3.1536e9*VX2D;
                    VY2D    = Dscale*transit_M/transit_E*3.1536e9*VY2D;
                end
                time1   = (time-t0)*transit_E/transit_M*3.170979e-17; % in Gyr
                % Composition is always non-dimensional
            end
            
            [h,time_str] = function_project3Dfield(DATA2D,thetad,phid,FontSize,ncontour,time1,itf,magic,FigSize,VX2D,VY2D,arrowscale,modulo,projectVEL,cminmax,dimensional_units,log,1);
            fig_num  = 1000000*magic+itf;
            save_pdf = true; 
            save_eps = false;
            save_jpg = false;
            save_png = false;
            save_str = strcat(dir,fname,'_n',nfn,'_projectC');  
% NICO            [ierr]   = savecoolplots(save_str,fig_num,save_pdf,save_eps,save_jpg,save_png);
            close(h)
                     
            if (save_projectC_data)
                [DAT] = savedata(DATA2D,thetad,phid);
                save_string = strcat(dir,fname,'_n',nfn,'_t',time_str,'_projectC.dat');
                save(save_string,'DAT','-ascii');        
            end
        end       
        
        if (projectT)
            
            resT_3D = zeros(nx,ny,nz,2);
            
            
            hormeanT = 0.0;
            horvol   = 0.0;
            %Load volumes of celles
            for ib = 1:2
                for iz = 1:nz
                    for iy = 1:ny
                        for ix = 1:nx
                            l = (ib-1)*nx*ny*nz + (iz-1)*nx*ny + (iy-1)*nx + ix;
                            dvol(ix,iy,iz)    = 1.;%Grid(l,1);
                        end
                    end
                end
            end
                      
            iz = ndepth;
            
            for iy = 1:ny
                for ix = 1:nx
                    hormeanT = hormeanT + (T_3D(ix,iy,iz,1) + T_3D(ix,iy,iz,2))*dvol(ix,iy,iz);
                    horvol   = horvol + 2*dvol(ix,iy,iz);
                end
            end
            
            hormeanT = hormeanT / horvol
            resT_3D(:,:,iz,1) = (T_3D(:,:,iz,1) - hormeanT)/hormeanT;
            resT_3D(:,:,iz,2) = (T_3D(:,:,iz,2) - hormeanT)/hormeanT;
                

            
            %DATAYY(:,:,:)      = T_3D(:,:,ndepth,:); 
            DATAYY(:,:,:)      = resT_3D(:,:,iz,:);
            [DATA2D,theta,phi] = YYtoMap2(DATAYY);
            thetad             = 90.-theta*180./3.14159;
            phid               = phi*180./3.14159;
            magic              = 2;
            cminmax            = [0 max(max(DATA2D))]; 
            time1              = time;
            log                = false;
            
            if (projectVEL) % add quiver?
                [VX2D,VY2D] = transform_VEL(dir,fname,itf);                
            end            
            
            if (dimensional_units)
                DATA2D  = DATA2D/Tscale*deltaT+Tsurf;
                if (projectVEL)
                    VX2D    = Dscale*transit_M/transit_E*3.1536e9*VX2D;
                    VY2D    = Dscale*transit_M/transit_E*3.1536e9*VY2D;
                end
                time1   = (time-t0)*transit_E/transit_M*3.170979e-17;
                cminmax = [Tsurf max(max(DATA2D))];
            end
            
            [h,time_str] = function_project3Dfield(DATA2D,thetad,phid,FontSize,ncontour,time1,itf,magic,FigSize,VX2D,VY2D,arrowscale,modulo,projectVEL,cminmax,dimensional_units,log,1);
            fig_num  = 1000000*magic+itf;
            save_pdf = true; 
            save_eps = false;
            save_jpg = false;
            save_png = false;
            save_str = strcat(dir,fname,'_n',nfn,'_projectResT');  
% NICO            [ierr]   = savecoolplots(save_str,fig_num,save_pdf,save_eps,save_jpg,save_png);
            close(h)
         
            if (save_projectT_data)
                ddepth= num2str(ndepth);
                [DAT] = savedata(DATA2D,thetad,phid);
                save_string = strcat(dir,fname,'_d',ddepth,'_project_ResT.dat');
                if (dimensional_units); save_string = strcat(dir,fname,'_d',ddepth,'_project_dim_ResT.dat'); end
                save(save_string,'DAT','-ascii');        
            end
        end
        
        if (projectQ)
            DATAYY(:,:,:)      = fg2D(:,:,:); 
            [DATA2D,theta,phi] = YYtoMap2(DATAYY);
            thetad             = 90.-theta*180./3.14159;
            phid               = phi*180./3.14159;
            magic              = 3;
            cminmax            = [0 max(max(DATA2D))];
            time1              = time;
            log                = false;
            
            if (projectVEL) % add quiver?
                [VX2D,VY2D] = transform_VEL(dir,fname,itf);                      
            end
            
            if (dimensional_units)
                if (projectVEL)
                    VX2D    = Dscale*transit_M/transit_E*3.1536e9*VX2D;
                    VY2D    = Dscale*transit_M/transit_E*3.1536e9*VY2D;
                end
                DATA2D  = DATA2D*k0*deltaT/mean(subTo(nstart:nend))/sqrt(kappa0*transit_E/transit_M)*1.e3;
                time1   = (time-t0)*transit_E/transit_M*3.170979e-17;
                cminmax = [0 200];
            end
            
            [h,time_str] = function_project3Dfield(DATA2D,thetad,phid,FontSize,ncontour,time1,itf,magic,FigSize,VX2D,VY2D,arrowscale,modulo,projectVEL,cminmax,dimensional_units,log,1);            
            fig_num  = 1000000*magic+itf;
            save_pdf = true; 
            save_eps = false;
            save_jpg = false;
            save_png = false;
            save_str = strcat(dir,fname,'_n',nfn,'_projectQ'); 
% NICO            [ierr]   = savecoolplots(save_str,fig_num,save_pdf,save_eps,save_jpg,save_png);
            close(h)
         
            if (save_projectQ_data)
                [DAT] = savedata(DATA2D,thetad,phid);
                save_string = strcat(dir,fname,'_t',time_str,'_project_Q.dat');
                if (dimensional_units); save_string = strcat(dir,fname,'_t',time_str,'_project_dim_Q.dat');  end
                save(save_string,'DAT','-ascii');        
            end
        end
             
        if (projectAGE)
            DATAYY(:,:,:)      = AGE_3D(:,:,nz,:);%AGEYY(:,:,:); 
            [DATA2D,theta,phi] = YYtoMap2(DATAYY);
            thetad             = 90.-theta*180./3.14159; %theta*180./3.14159;
            phid               = phi*180./3.14159; %phi*180./3.14159;
            magic              = 4;
            cminmax            = [0 300];
            time1              = time;
            log                = false;
            
            if (projectVEL) % add quiver?
                [VX2D,VY2D] = transform_VEL(dir,fname,itf);               
            end
            
            if (dimensional_units)
                % Ages are dimensional already
                if (projectVEL)
                    VX2D    = Dscale*transit_M/transit_E*3.1536e9*VX2D;
                    VY2D    = Dscale*transit_M/transit_E*3.1536e9*VY2D;
                end
                time1   = (time-t0)*transit_E/transit_M*3.170979e-17;
            end
            
            [h,time_str] = function_project3Dfield(DATA2D,thetad,phid,FontSize,ncontour,time1,itf,magic,FigSize,VX2D,VY2D,arrowscale,modulo,projectVEL,cminmax,dimensional_units,log,1);                           
            fig_num  = 1000000*magic+itf;
            save_pdf = true; 
            save_eps = false;
            save_jpg = false;
            save_png = false;
            save_str = strcat(dir,fname,'_n',nfn,'_projectAGE');  
 % NICO           [ierr]   = savecoolplots(save_str,fig_num,save_pdf,save_eps,save_jpg,save_png);
            close(h)
         
            if (save_projectAGE_data)
                [DAT]    = savedata(DATA2D,thetad,phid);
                time1s   = num2str(time1);
                save_str = strcat(dir,fname,'-age',nfn,'.dat');
                save(save_str,'DAT','-ascii');        
            end           
        end
        
        if (projectBATH)
            DATAYY(:,:,:)      = bathYY(:,:,:); 
            [DATA2D,theta,phi] = YYtoMap2(DATAYY);
            thetad             = 90.-theta*180./3.14159;
            phid               = phi*180./3.14159;
            magic              = 5;
            cminmax            = [bath_min bath_max];
            time1              = time;
            log                = false;
            
            if (projectVEL) % add quiver?
                [VX2D,VY2D] = transform_VEL(dir,fname,itf);
            end
            
            if (dimensional_units)
                % Bathymetry is dimensional already
                if (projectVEL)
                    VX2D    = Dscale*transit_M/transit_E*3.1536e9*VX2D;
                    VY2D    = Dscale*transit_M/transit_E*3.1536e9*VY2D;
                end
                time1   = (time-t0)*transit_E/transit_M*3.170979e-17;
            end
                        
            [h,time_str] = function_project3Dfield(DATA2D,thetad,phid,FontSize,ncontour,time1,itf,magic,FigSize,VX2D,VY2D,arrowscale,modulo,projectVEL,cminmax,dimensional_units,log,1);
            fig_num  = 1000000*magic+itf;
            save_pdf = true; 
            save_eps = false;
            save_jpg = false;
            save_png = false;
            save_str = strcat(dir,fname,'_n',nfn,'_projectBATH');  
 % NICO           [ierr]   = savecoolplots(save_str,fig_num,save_pdf,save_eps,save_jpg,save_png);
            close(h)
         
            if (save_projectBATH_data)
                [DAT]    = savedata(DATA2D,thetad,phid);
                save_str = strcat(dir,fname,'_t',time_str,'_project_dim_Bath.dat');
                save(save_str,'DAT','-ascii');       
            end
        end
        
        if (projectETA)
            DATAYY(:,:,:)      = ETA_3D(:,:,nz,:); 
            [DATA2D,theta,phi] = YYtoMap2(DATAYY);
            thetad             = 90.-theta*180./3.14159;
            phid               = phi*180./3.14159;
            magic              = 6;
            cminmax            = [0.1 1e+6];
            time1              = time;
            log                = true;
            
            if (projectVEL) % add quiver?
                [VX2D,VY2D] = transform_VEL(dir,fname,itf);
            end
            
            if (dimensional_units)
                DATA2D  = DATA2D*Eta_Ref;
          %Nicolas      if (projectVEL)                  
          %          VX2D    = Dscale*transit_M/transit_E*3.1536e9*VX2D;
          %          VY2D    = Dscale*transit_M/transit_E*3.1536e9*VY2D;
          %     end
                cminmax = [1.e+21 1.e+28];
                time1   = (time-t0)*transit_E/transit_M*3.170979e-17;
            end
                     
         %Nicolas   [h,time_str] = function_project3Dfield(DATA2D,thetad,phid,FontSize,ncontour,time1,itf,magic,FigSize,VX2D,VY2D,arrowscale,modulo,projectVEL,cminmax,dimensional_units,log,1);
          %  fig_num  = 1000000*magic+itf;
          %  save_pdf = true; 
          %  save_eps = false;
          %  save_jpg = false;
          %  save_png = false;
           % save_str = strcat(dir,fname,'_n',nfn,'_projectETA');  
 % NICO           [ierr]   = savecoolplots(save_str,fig_num,save_pdf,save_eps,save_jpg,save_png);
  %          close(h)
         
            if (save_projectETA_data)
                [DAT] = savedata(DATA2D,thetad,phid);
                save_string = strcat(dir,fname,nfn,'_project_Eta.dat');
                if (dimensional_units); save_string = strcat(dir,fname,'_project_dim_Eta.dat'); end
                save(save_string,'DAT','-ascii');
            end
        end
        
        if (projectDIV)
            [~,~,~,~,DIV_3D(:,:,:,1),DIV_3D(:,:,:,2),~,~] = ReadStag3D2_251012(file_stem,itf,'divergence',false ,[1 1 1]);
            DATAYY(:,:,:)      = DIV_3D(:,:,nz,:); % surface layer only 
            [DATA2D,theta,phi] = YYtoMap2(DATAYY);
            thetad             = 90.-theta/3.14159*180. %theta*180./3.14159;
            phid               = phi/3.14159*180.;%phi*180./3.14159;
            magic              = 13;
            cminmax            = [-0.1/kappa0 0.1/kappa0]; % might need adjustment
            time1              = time;
            log                = false;
            
            if (projectVEL) % add quiver?
                [VX2D,VY2D] = transform_VEL(dir,fname,itf);
                %Nicolas save data
                [THD PHD] = meshgrid(thetad,phid);
                ntheta    = length(thetad);
                nphi      = length(phid);

                DAT  = zeros(ntheta*nphi,4);
                for ith = 1:ntheta
                    for iph = 1:nphi
                        ilg         = (ith-1)*nphi + iph;
                        DAT(ilg,1)  = PHD(iph,ith);
                        DAT(ilg,2)  = THD(iph,ith);
                        DAT(ilg,3)  = VX2D(ith,iph);
                        DAT(ilg,4)  = VY2D(ith,iph);
                    end
                end
 
                save_string = strcat(dir,fname,nfn,'_project_VEL.dat');
                if (dimensional_units); save_string = strcat(dir,fname,'-vel',nfn,'.dat'); end
                save(save_string,'DAT','-ascii');               
    
            end
            
            if (dimensional_units)
                DATA2D  = DATA2D*kappa0;
               if (projectVEL)                  
                    VX2D    = Dscale*transit_M/transit_E*3.1536e9*VX2D;
                    VY2D    = Dscale*transit_M/transit_E*3.1536e9*VY2D;
                end
                cminmax = [-0.1 0.1]; % might need adjustment
                time1   = (time-t0)*transit_E/transit_M*3.170979e-17;
            end
                     
        %Nicolas    [h,time_str] = function_project3Dfield(DATA2D,thetad,phid,FontSize,ncontour,time1,itf,magic,FigSize,VX2D,VY2D,arrowscale,modulo,projectVEL,cminmax,dimensional_units,log,1);
        %    fig_num  = 1000000*magic+itf;
        %    save_pdf = true; 
        %    save_eps = false;
        %    save_jpg = false;
        %    save_png = false;
        %    save_str = strcat(dir,fname,'_n',nfn,'_projectDIV');  
   % NICO             [ierr]   = savecoolplots(save_str,fig_num,save_pdf,save_eps,save_jpg,save_png);
    %    close(h)
         
            if (save_projectDIV_data)
                [DAT] = savedata(DATA2D,thetad,phid);
                save_string = strcat(dir,'serie1',nfn,'_project_DIV.dat');
                if (dimensional_units); save_string = strcat(dir,fname,'-div',nfn,'.dat'); end
                DATNN=DAT(:,3);
                save(save_string,'DATNN','-ascii');
            end
        end
        
        if (projectVOR)
            [~,~,~,~,VOR_3D(:,:,:,1),VOR_3D(:,:,:,2),~,~] = ReadStag3D2_251012(file_stem,itf,'vorticity',false ,[1 1 1]);
            DATAYY(:,:,:)      = VOR_3D(:,:,nz,:); % surface layer only 
            [DATA2D,theta,phi] = YYtoMap2(DATAYY);
            thetad             = 90.-theta*180./3.14159;%-theta*180./3.14159;
            phid               = phi*180./3.14159;%phi*180./3.14159;
            magic              = 14;
            cminmax            = [-0.1/kappa0 0.1/kappa0]; % might need adjustment
            time1              = time;
            log                = false;
            
            if (projectVEL) % add quiver?
                [VX2D,VY2D] = transform_VEL(dir,fname,itf);
            end
            
            if (dimensional_units)
                DATA2D  = DATA2D*kappa0;
                if (projectVEL)                  
                    VX2D    = Dscale*transit_M/transit_E*3.1536e9*VX2D;
                    VY2D    = Dscale*transit_M/transit_E*3.1536e9*VY2D;
                end
                cminmax = [-0.1 0.1]; % might need adjustment
                time1   = (time-t0)*transit_E/transit_M*3.170979e-17; % in Gyr
            end
                     
            [h,time_str] = function_project3Dfield(DATA2D,thetad,phid,FontSize,ncontour,time1,itf,magic,FigSize,VX2D,VY2D,arrowscale,modulo,projectVEL,cminmax,dimensional_units,log,1);
            fig_num  = 1000000*magic+itf;
            save_pdf = true; 
            save_eps = false;
            save_jpg = false;
            save_png = false;
            save_str = strcat(dir,fname,'_n',nfn,'_projectVOR');  
 % NICO           [ierr]   = savecoolplots(save_str,fig_num,save_pdf,save_eps,save_jpg,save_png);
            close(h)
         
            if (save_projectVOR_data)
                [DAT] = savedata(DATA2D,thetad,phid);
                save_string = strcat(dir,fname,'_t',time_str,'_project_VOR.dat');
                if (dimensional_units); save_string = strcat(dir,fname,'-vor',nfn,'.dat'); end

                save(save_string,'DAT','-ascii');
            end
        end
                           
        case 'Annullus'
            
            if (itf == nstart)
                warning('Age Distribution is probably biased in 2D!')
                
                % 2D benchmark of heat flow
                % load_str = strcat(dir,fname,'_time.dat');
                % HFDAT    = load(load_str);
                % hftime   = HFDAT(:,2);
                % hfdata   = HFDAT(:,3);
                               
                [nnb,~,~,Z_3D,T_3D,time,rcmb] = ReadStag3D2_251012(file_stem,itf,'temperature',false,[1 1 1]);
                 
                nx = 1;
                ny = size(T_3D,2);
                nz = size(T_3D,3);         
            
                % CALCULATE SURFCE CELL AREAS - EASY, ALL CELLS HAVE THE SAME SIZE
                
                darea      = zeros(nx,ny);           
                darea(:,:) = 2*pi*Radius/ny;
                dimarea    = Radius;
                
                % INITIALIZE COMPOSITIONAL FIELDS 
                
                C_3D   = zeros(nx,ny,nz,nnb);
                C1_3D  = zeros(nx,ny,nz,nnb);             
                NRC_3D = zeros(nx,ny,nz,nnb);
                sumC   = zeros(nx,ny,1,nnb);
                
            else
                tic
               [~,~,~,~,T_3D,time,~] = ReadStag3D2_290611(dir,fname,itf,'temperature');
            end
                 
            if WriteComposition;  [~,~,~,~,C_3D,~,~]   = ReadStag3D2_251012(file_stem,itf,'cont root' ,false ,[1 1 1]);  end             
            if WriteBelts;        [~,~,~,~,C1_3D,~,~]  = ReadStag3D2_251012(file_stem,itf,'cont crust',false ,[1 1 1]);  end             
            if WriteContinent;    [~,~,~,~,NRC_3D,~,~] = ReadStag3D2_251012(file_stem,itf,'continent' ,false ,[1 1 1]);  end 
            
            
            sumC(:,:,1,:) = C_3D(:,:,nz) + C1_3D(:,:,nz,:);
     
            % CALCULATION OF SURFACE HEAT FLUX VIA FOURIERS LAW

            qg           = 0.0;
            qo           = 0.0;
            qc           = zeros(ncont+1,1);

            fo           = 0.0;
            fg           = 0.0;
            fc           = zeros(ncont+1,1);
             
            qg_area      = 0.0;  
            qo_area      = 0.0;  
            qc_area      = zeros(ncont+1,1);
        
            Ts           = zeros(nx,ny,2);    
        
            barea        = zeros(length(bhist),1);
            harea        = zeros(length(xhist),1);
            harea_dim    = zeros(length(xhist_dim),1);
   
            domain_depth = Rsurf-rcmb;
            dzf          = 2*(domain_depth - Z_3D(1,1,nz));
     
            for iy = 1:ny
                for ix = 1:nx
                    % surface temperature (do it like Paul in STAG's subroutine topbotflux
                    Ts(ix,iy,1) = - T_3D(ix,iy,nz);
                    Ts(ix,iy,2) = + T_3D(ix,iy,nz);                         

                    % global heatflux
                    fg            = (Ts(ix,iy,2) - Ts(ix,iy,1))/dzf;
                    qg            = qg           + fg*darea(ix,iy);
                    qg_area       = qg_area      + darea(ix,iy);

                    % oceanic heatflux
                    if (sumC(ix,iy,1,1) <= Cthresh)
                        fo          = (Ts(ix,iy,2)-Ts(ix,iy,1))/dzf;

                        if (fo>0. && darea(ix,iy)>0.);                 
                            qo         = qo + fo*darea(ix,iy);
                            qo_area    = qo_area + darea(ix,iy);

                            if (itf > 0)
                                age(1) = subTo(itf)^2/(pi*fo*fo);
                                Toc    = subTo(itf)*deltaT/Tscale + Tsurf;
                            else
                                age(1) = subTo(itf+1)^2/(pi*fo*fo);
                                Toc    = subTo(itf+1)*deltaT/Tscale + Tsurf; 
                            end

                            age(2) = darea(ix,iy);
                            age(3) = age(1)*transit_E/transit_M;

                            % Calculate the bathymetry in that cell from the age 
                            AGEYY(ix,iy)  = age(3)*Myr;
                            bathYY(ix,iy) = 2*alpha0*rho_m/drho*(Toc-Tsurf)*sqrt(kappa0/Myr/pi)*sqrt(AGEYY(ix,iy));

                            for ixh = 1:length(xhist)-1 % non-dimensional
                                if (age(1) > xhist(ixh) && age(1) <= xhist(ixh+1))
                                    harea(ixh)   = harea(ixh) + age(2);
                                end
                            end

                            for ixh = 1:length(xhist_dim)-1 % dimensional                        
                                if (age(3) > xhist_dim(ixh) && age(3) <= xhist_dim(ixh+1))
                                    harea_dim(ixh) = harea_dim(ixh) + age(2)*dimarea;
                                    na(ixh)        = na(ixh) + 1;
                                end
                            end                          

                            for ibh = 1:length(bhist)-1  % same for bathymetry
                                if (bathYY(ix,iy) > bhist(ibh) && bathYY(ix,iy) <= bhist(ibh+1))
                                    barea(ibh) = barea(ibh) + darea(ix,iy)*dimarea;
                                    nb(ibh)    = nb(ibh) + 1;
                                end
                            end
                        end
                    end           
                    % continental heatflux
                    if (sumC(ix,iy,1,1) >= Cthresh)
                       fc(ncont+1)              = (Ts(ix,iy,2)     - Ts(ix,iy,1))/dzf;
                       qc(ncont+1)              = qc(ncont+1)      + fc(ncont+1)*darea(ix,iy);
                       qc_area(ncont+1)         = qc_area(ncont+1) + darea(ix,iy);

                        for icont = 1:ncont
                            if (NRC_3D(ix,iy,nz) == icont) 
                                fc(icont)       = (Ts(ix,iy,2)   - Ts(ix,iy,1))/dzf;
                                qc(icont)       = qc(icont)      + fc(icont)*darea(ix,iy);
                                qc_area(icont)  = qc_area(icont) + darea(ix,iy);
                            end
                        end
                    end           
                end
            end
     
            qg     = qg / qg_area;
            qo     = qo / qo_area;
        
            for icont = 1:ncont+1
                qc(icont) = qc(icont) / qc_area(icont);
            end
        
       case 'Cartesian'
            error('Cartesian geometry not yet implemented!')
   end
   
    % UPDATE COUNTER
    num  = num+1;
    tic
end      % END-OF_TIMESTEP LOOP

