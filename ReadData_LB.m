%ReadData
%
% Read binary output files of STAG3D and write them into Paraview files
%
% Boris Kaus, 26.2.2008

% Modified by Lea to include total_time and dimensional time in the pvd
% files.
%


% PJT edits August 2010:
% - spherical: move coordinate transform just before writing, otherwise reset
% - spherical: centre domain around equator, rotate velocity components
% - YinYang: transform velocity components
% - new Annulus mode (replaced previous 'Cylindrical' mode)

% TR edits June 2011:
% - Continents

clear

directory='';%
fname ='';%

%GridType       =   'Cartesian';
%GridType       =   'Annulus';         
GridType       =   'YinYang';

number_start   = 100;
number_end     = 100;                                                                          
incr           = 1;
size_reduction = true; 
nr_save         = [4 4 1];
nrstart        = 63 %starting at which depth?

if (size_reduction)
    display('Size reduction of *.vtu files is switched on, not all data points will be stored')
end

% Specify which fields we want to read/write in addition to T:
WriteResidualT    =   false;  % Residual temp, i.e. after subtracting horiz. mean T
WriteViscosity    =   true;  % Viscosity                         
WriteVelocity     =   false;  % Velocity and Pressure 
WriteStress       =   false;  % Convective Stresses
WriteEdot         =   false;  % Strain rate
WriteHeatflux     =   false;
WriteContRoot     =   false;  % continental root     ...always true
WriteComposition  =   false;  % some compositional field
WriteTopography   =   false;  % Togography 
WriteGeoid        =   false;  % Geoid from sph. harmonics expansion
WriteContCrust    =   false;  % continental crust    ...actually it can be belt/suture as well, it is just the 2nd C-Field in STAG
WriteContBelt     =   false;  % continental margin   ...actually it can be be belt or suture
WriteContSuture   =   false;  % continental suture   ...always true
WriteContPlot     =   false;  % allows for visualizing all compositions at the same time
WriteContinent    =   false;  % continental number
WriteTopoSG       =   false;  % Self-gravitatad Topography
WriteMeltFrac     =   false;  % Melt fraction
WriteResidue      =   false;  % Residue
WriteDamage       =   false;  % Damage Evolution
WriteAge          =   false;
WriteDiv          =   false;
WriteVor          =   false;

num = 1; % number of frame for name array

for fname_number=number_start:incr:number_end
  fname_number
    switch GridType
%%        
        case 'YinYang'

            % YING YANG GRID
            disp(['Creating a YIN-YANG grid ...'])

             if (size_reduction)
                [nnb, X_3D, Y_3D, Z_3D,T_3D_1, T_3D_2, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,  'temperature',size_reduction,nr_save,nrstart);
             else
                [nnb, X_3D, Y_3D, Z_3D,T_3D_1, T_3D_2, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'temperature');
             end
             
             nnb
             rcmb
             nx = size(X_3D,1);
             ny = size(Y_3D,2);
             nz = size(Z_3D,3);      

             if nnb==-999
                % The file does not exist, and we have finished processing all
                % data

                % Create a PVD file, which contains all data including
                % time-information
                Create_PVD_file(FileNames,fname,directory);

                error(['finished processing all files in directory'])

             else
                 if WriteResidualT
                %    if (size_reduction)
                %        error('Grid size reduction not implemented for residualT')
                %    end
                    
%                     if (nx==64 && ny==192 && nz==32)
%                         load 'Cell_vol_YY_64_192_32_T'
%                         Grid = Cell_vol_YY_64_192_32_T;
%                     elseif (nx==64 && ny==192 && nz==64)
%                         load 'Cell_vol_YY_64_192_64_T'
%                         Grid = Cell_vol_YY_64_192_64_T;
%                     elseif (nx==96 && ny==288 && nz==48)
%                         load 'Cell_vol_YY_96_288_48_T'
%                         Grid = Cell_vol_YY_96_288_48_T;
%                     elseif (nx==128 && ny==384 && nz==64)
%                         load 'Cell_vol_128_384_64_T'
%                         Grid = Cell_vol_YY_128_384_64_T;
%                     else
%                         error('Cannot load grid cell data - misssing file')
%                     end
            
                    % reloading the cell volumes for the YinYang case
                    % to be improved - better is to calculate those from  
                    % X,Y,Z coordinates via volume integration
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
                    
                    % calculate horizontal mean and residual temperature
                    
                    resT_3D_1 = zeros(nx,ny,nz);
                    resT_3D_2 = zeros(nx,ny,nz);
                    
                    for iz = 1:nz                       
                        hormeanT = 0.0;
                        horvol   = 0.0;
                        
                        for iy = 1:ny                            
                            for ix = 1:nx                                
                                hormeanT = hormeanT + (T_3D_1(ix,iy,iz) + T_3D_2(ix,iy,iz))*dvol(ix,iy,iz);
                                horvol   = horvol + 2*dvol(ix,iy,iz);
                            end                          
                        end
                        
                        hormeanT = hormeanT / horvol;
                        resT_3D_1(:,:,iz) = T_3D_1(:,:,iz) - hormeanT;
                        resT_3D_2(:,:,iz) = T_3D_2(:,:,iz) - hormeanT;
                    end

                end
                
                if WriteVelocity
                    % Read pressure & velocity information
                    if (size_reduction) 
                       [nnb, X_3D, Y_3D, Z_3D,VX_3D_1, VY_3D_1, VZ_3D_1, P_3D_1, VX_3D_2, VY_3D_2, VZ_3D_2, P_3D_2, time] = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'velocity',size_reduction,nr_save,nrstart); 
                    else
                       [nnb, X_3D, Y_3D, Z_3D,VX_3D_1, VY_3D_1, VZ_3D_1, P_3D_1, VX_3D_2, VY_3D_2, VZ_3D_2, P_3D_2, time] = ReadStag3D_LeaV(directory, fname, fname_number,  'velocity');                           
                    end

                end
                
                if WriteResidue
                    % Read Residue
                    if (size_reduction)
                       [nnb, X_3D, Y_3D, Z_3D,RX_3D_1, RY_3D_1, RZ_3D_1, RP_3D_1, RX_3D_2, RY_3D_2, RZ_3D_2, RP_3D_2, time] = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'residue',size_reduction,nr_save,nrstart);
                    else                     
                       [nnb, X_3D, Y_3D, Z_3D,RX_3D_1, RY_3D_1, RZ_3D_1, RP_3D_1, RX_3D_2, RY_3D_2, RZ_3D_2, RP_3D_2, time] = ReadStag3D_LeaV(directory, fname, fname_number,  'residue');                           
                
                    end
                end
                
                if WriteViscosity
                    if (size_reduction)
                        [nnb, X_3D, Y_3D, Z_3D,ETA_3D_1, ETA_3D_2, time, rcmb] = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'viscosity',size_reduction,nr_save,nrstart);
                    else
                        [nnb, X_3D, Y_3D, Z_3D,ETA_3D_1, ETA_3D_2, time, rcmb] = ReadStag3D_LeaV(directory, fname, fname_number,  'viscosity');
                    end
                end
                
                if WriteComposition
                    if (size_reduction)
                        [nnb, X_3D, Y_3D, Z_3D, C_3D_1, C_3D_2, time, rcmb]   = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'composition',size_reduction,nr_save,nrstart);  
                    else
                        [nnb, X_3D, Y_3D, Z_3D, C_3D_1, C_3D_2, time, rcmb]   = ReadStag3D_LeaV(directory, fname, fname_number,  'composition');
                    end
                end
                
                if WriteContRoot
                    if (size_reduction)
                        [nnb, X_3D, Y_3D, Z_3D, C1_3D_1, C1_3D_2, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'cont root',size_reduction,nr_save,nrstart);
                    else
                        [nnb, X_3D, Y_3D, Z_3D, C1_3D_1, C1_3D_2, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'cont root');
                    end
                end
                
                if WriteContCrust
                    if (size_reduction)
                        [nnb, X_3D, Y_3D, Z_3D, C2_3D_1, C2_3D_2, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'cont crust',size_reduction,nr_save,nrstart);
                    else
                        [nnb, X_3D, Y_3D, Z_3D, C2_3D_1, C2_3D_2, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'cont crust');
                    end
                end
                
                if WriteContBelt
                    if (size_reduction)
                        [nnb, X_3D, Y_3D, Z_3D, C3_3D_1, C3_3D_2, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'cont belt',size_reduction,nr_save,nrstart);
                    else
                        [nnb, X_3D, Y_3D, Z_3D, C3_3D_1, C3_3D_2, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'cont belt');
                    end
                end
                
                if WriteContSuture
                    if (size_reduction)
                        [nnb, X_3D, Y_3D, Z_3D, C4_3D_1, C4_3D_2, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'cont suture',size_reduction,nr_save,nrstart);
                    else
                        [nnb, X_3D, Y_3D, Z_3D, C4_3D_1, C4_3D_2, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'cont suture');
                    end
                end
                
                if WriteContPlot
                    if (size_reduction)
                        [nnb, X_3D, Y_3D, Z_3D, CPL_3D_1, CPL_3D_2, time, rcmb]   = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'cont plot',size_reduction,nr_save);
                    else
                        [nnb, X_3D, Y_3D, Z_3D, CPL_3D_1, CPL_3D_2, time, rcmb]   = ReadStag3D_LeaV(directory, fname, fname_number,  'cont plot');
                    end
                end
                
                if WriteContinent
                    if (size_reduction)                  
                        [nnb, X_3D, Y_3D, Z_3D, NRC_3D_1, NRC_3D_2, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'continent',size_reduction,nr_save,nrstart);
                    else
                        [nnb, X_3D, Y_3D, Z_3D, NRC_3D_1, NRC_3D_2, time, rcmb]   = ReadStag3D_LeaV(directory, fname, fname_number,  'continent');
                    end
                end
                
                if WriteTopography
                    % Topography is a 2D Field
                    if (size_reduction)
                        % update dimensions
                        nx = size(X_3D,1);
                        ny = size(Y_3D,2);
                        nz = size(Z_3D,3);
                    end
                    TP_3D_1 = zeros(nx,ny,nz);
                    TP_3D_2 = zeros(nx,ny,nz);
                    
                    if (size_reduction)
                        [ ~, ~, ~, ~,TP1_3D_1, TP1_3D_2, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'topography',size_reduction,nr_save,nrstart);
                    else
                        [ ~, ~, ~, ~,TP1_3D_1, TP1_3D_2, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'topography'); 
                    end
                    % CAUTION: check if surface and cmb topography are in
                    % the right order ...same for 2D cases

                    TP_3D_1(:,:,1)  = TP1_3D_1(:,:,1);
                    TP_3D_1(:,:,nz) = TP1_3D_1(:,:,2);
                    TP_3D_2(:,:,1)  = TP1_3D_2(:,:,1);
                    TP_3D_2(:,:,nz) = TP1_3D_2(:,:,2);
                end                
                
                if WriteTopoSG
                    % Topography Self-grav is a 2D Field
                    TPSG_3D_1 = zeros(nx,ny,nz);
                    TPSG_3D_2 = zeros(nx,ny,nz);
                    if (size_reduction)
                        [ ~, ~, ~, ~,TPSG1_3D_1, TPSG1_3D_2, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'topography self-grav',size_reduction,nr_save,nrstart); 
                    else
                        [ ~, ~, ~, ~,TPSG1_3D_1, TPSG1_3D_2, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'topography self-grav'); 
                    end
                    % CAUTION: check if surface and cmb topography are in
                    % the right order ...same for 2D cases
                    TPSG_3D_1(:,:,1)  = TPSG1_3D_1(:,:,1);
                    TPSG_3D_1(:,:,nz) = TPSG1_3D_1(:,:,2);
                    TPSG_3D_2(:,:,1)  = TPSG1_3D_2(:,:,1);
                    TPSG_3D_2(:,:,nz) = TPSG1_3D_2(:,:,2);
                end
                
                
                if WriteGeoid
                    % Topography is a 2D Field
                    G_3D_1 = zeros(nx,ny,nz);
                    G_3D_2 = zeros(nx,ny,nz);
                    if (size_reduction)
                        [ ~, ~, ~, ~,G1_3D_1, G1_3D_2, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'geoid',size_reduction,nr_save,nrstart);
                    else
                        [ ~, ~, ~, ~,G1_3D_1, G1_3D_2, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'geoid');
                    end
                    % CAUTION: check if surface and cmb topography are in
                    % the right order ...same for 2D cases
                    G_3D_1(:,:,1)  = G1_3D_1(:,:,1);
                    G_3D_1(:,:,nz) = G1_3D_1(:,:,2);
                    G_3D_2(:,:,1)  = G1_3D_2(:,:,1);
                    G_3D_2(:,:,nz) = G1_3D_2(:,:,2);
                end
                
                if WriteHeatflux
                    % Topography is a 2D Field
                    HF_3D_1 = zeros(nx,ny,nz);
                    HF_3D_2 = zeros(nx,ny,nz);
                    if (size_reduction)
                        [ ~, ~, ~, ~,HF1_3D_1, HF1_3D_2, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'heat flux',size_reduction,nr_save,nrstart);
                    else
                        [ ~, ~, ~, ~,HF1_3D_1, HF1_3D_2, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'heat flux');
                    end
                    % CAUTION: check if surface and cmb topography are in
                    % the right order ...same for 2D cases
                    HF_3D_1(:,:,1)  = HF1_3D_1(:,:,1);
                    HF_3D_1(:,:,nz) = HF1_3D_1(:,:,2);
                    HF_3D_2(:,:,1)  = HF1_3D_2(:,:,1);
                    HF_3D_2(:,:,nz) = HF1_3D_2(:,:,2);
                end
                
                        
                if WriteStress
                    if (size_reduction)
                        [nnb, X_3D, Y_3D, Z_3D, S_3D_1, S_3D_2, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'stress',size_reduction,nr_save,nrstart);
                    else
                        [nnb, X_3D, Y_3D, Z_3D,S_3D_1, S_3D_2, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'stress');
                    end
                end
                
                if WriteEdot
                    if (size_reduction)
                        [nnb, X_3D, Y_3D, Z_3D, E_3D_1, E_3D_2, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'strain rate',size_reduction,nr_save,nrstart);
                    else
                        [nnb, X_3D, Y_3D, Z_3D, E_3D_1, E_3D_2, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'strain rate');
                    end
                end
                
                if WriteDamage
                    if (size_reduction)
                        [nnb, X_3D, Y_3D, Z_3D, D_3D_1, D_3D_2, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'damage',size_reduction,nr_save,nrstart);
                    else
                        [nnb, X_3D, Y_3D, Z_3D, D_3D_1, D_3D_2, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'damage');
                    end
                end

                %xs = X_3D(:,:,end)+pi/4;
                %ys = Y_3D(:,:,end)-pi/4;
                %zs = Z_3D(:,:,end);

                % Transform coordinates for Yin & Yang grids
                R  = Z_3D+rcmb;     lat = pi/4-X_3D;  lon = Y_3D-3*pi/4;
                
                X_3D_1 = R.*cos(lat).*cos(lon);   % Yin grid
                Y_3D_1 = R.*cos(lat).*sin(lon);
                Z_3D_1 = R.*sin(lat);             
                
                X_3D_2 = -X_3D_1;       %  Yang grid
                Y_3D_2 =  Z_3D_1;
                Z_3D_2 =  Y_3D_1;
                
                Xs1=double(X_3D_1);
                Ys1=double(Y_3D_1);
                Zs1=double(Z_3D_1);
                Ts1=double(T_3D_1);

                
                % Transform velocities, if needed
                if WriteVelocity
                    Vtheta = VX_3D_1; Vphi = VY_3D_1; Vr = VZ_3D_1;                                 % on Yin grid
                    VX_3D_1 =  Vtheta.*sin(lat).*cos(lon) - Vphi.*sin(lon) + Vr.*cos(lat).*cos(lon);
                    VY_3D_1 =  Vtheta.*sin(lat).*sin(lon) + Vphi.*cos(lon) + Vr.*cos(lat).*sin(lon);
                    VZ_3D_1 =  -Vtheta.*cos(lat)                            + Vr.*sin(lat);
                    Vtheta = VX_3D_2; Vphi = VY_3D_2; Vr = VZ_3D_2;                                 % on Yang grid
                    VX_3D_2 = -( Vtheta.*sin(lat).*cos(lon) - Vphi.*sin(lon) + Vr.*cos(lat).*cos(lon) );
                    VZ_3D_2 =  Vtheta.*sin(lat).*sin(lon) + Vphi.*cos(lon) + Vr.*cos(lat).*sin(lon);
                    VY_3D_2 =  -Vtheta.*cos(lat)                            + Vr.*sin(lat);
                    
                   
                end
                
                if WriteAge
                    % Read Residue
                    if (size_reduction)
                       [nnb, X_3D, Y_3D, Z_3D, AGE_3D_1, AGE_3D_2, time, rcmb] = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'age',size_reduction,nr_save,nrstart);
                    else                     
                       [nnb, X_3D, Y_3D, Z_3D, AGE_3D_1, AGE_3D_2, time, rcmb] = ReadStag3D_LeaV(directory, fname, fname_number,   'age');                           
                
                    end
                end
                
                
                if WriteDiv
                    % Read Residue
                    if (size_reduction)
                       [nnb, X_3D, Y_3D, Z_3D, DIV_3D_1, DIV_3D_2, time, rcmb] = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'divergence',size_reduction,nr_save,nrstart);
                    else                     
                       [nnb, X_3D, Y_3D, Z_3D, DIV_3D_1, DIV_3D_2, time, rcmb] = ReadStag3D_LeaV(directory, fname, fname_number,   'divergence');                           
                
                    end
                end
                
                
                if WriteVor
                    % Read Residue
                    if (size_reduction)
                       [nnb, X_3D, Y_3D, Z_3D, VOR_3D_1, VOR_3D_2, time, rcmb] = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,   'vorticity',size_reduction,nr_save,nrstart);
                    else                     
                       [nnb, X_3D, Y_3D, Z_3D, VOR_3D_1, VOR_3D_2, time, rcmb] = ReadStag3D_LeaV(directory, fname, fname_number,   'vorticity');                           
                
                    end
                end
                  
                % Transform to dimensional time
                %time_dim(num)=200-time*vsmean*2900e3/(vrecmean*1e4);
                time_dim(num)=time;
                
                % Write a ying-yang grid
                WriteStag3D_VTK_YinYang_LB;    % Write VTK XML file
                fname_vtk       = [fname,'_',num2str(1000000+fname_number),'.vtu'];

                % Store the filename and time in a structure, which is later
                % used to create a pvd file, that contains all timesteps and
                % times of the files
                FileNames{num} = {time_dim(num),1,fname_vtk};
                %FileNames{num} = {time,fname_number,fname_vtk}; % ????

            end
%%
         case 'Cartesian'
            % CARTESIAN GRID
            
            % Read temperature information
            [nnb, X_3D, Y_3D, Z_3D,T_3D, time]     = ReadStag3D_LeaV(directory, fname, fname_number,  'temperature');

             if nnb==-999
                % The file does not exist, and we have finished processing all
                % data

                % Create a PVD file, which contains all data including
                % time-information
                Create_PVD_file(FileNames,fname,directory);

                error(['finished processing all files in directory'])

            else

                if WriteVelocity
                    % Read pressure & velocity information
                    [nnb,X_3D, Y_3D, Z_3D, VX_3D, VY_3D, VZ_3D, P_3D, time] = ReadStag3D_LeaV(directory, fname, fname_number,  'velocity');
                end

                if WriteViscosity
                    % Read viscosity information
                    [nnb,X_3D, Y_3D, Z_3D,ETA_3D, time] = ReadStag3D_LeaV(directory, fname, fname_number,  'viscosity');
                end
                
                if WriteComposition
                    [nnb,X_3D, Y_3D, Z_3D,C_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'composition');
                end
                
                if WriteContRoot
                    [nnb,X_3D, Y_3D, Z_3D,C1_3D, time]  = ReadStag3D_LeaV(directory, fname, fname_number,  'cont root');
                end
               
                if WriteContCrust
                    [nnb,X_3D, Y_3D, Z_3D,C2_3D, time]  = ReadStag3D_LeaV(directory, fname, fname_number,  'cont crust');
                end
                
                if WriteContBelt
                    [nnb,X_3D, Y_3D, Z_3D,C3_3D, time]  = ReadStag3D_LeaV(directory, fname, fname_number,  'cont belt');
                end
                
                if WriteContSuture
                    [nnb,X_3D, Y_3D, Z_3D,C4_3D, time]  = ReadStag3D_LeaV(directory, fname, fname_number,  'cont suture');
                end
                
                if WriteContPlot
                    [nnb,X_3D, Y_3D, Z_3D, CPL_3D, time]  = ReadStag3D_LeaV(directory, fname, fname_number,  'cont plot');
                end
                
                if WriteContinent
                    [nnb,X_3D, Y_3D, Z_3D, NRC_3D, time]  = ReadStag3D_LeaV(directory, fname, fname_number,  'continent');
                end
                
                if WriteTopography 
                    % Read topography information
                    % CAUTION: 2D field - dimensions doesn't match with the
                    % other fields --> add zeros
                    nx = size(X_3D,1);
                    ny = size(Y_3D,2);
                    nz = size(Z_3D,3);
                    TP_3D = zeros(nx,ny,nz);
                    [nnb,X_3D, Y_3D, Z_3D,TP1_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'topography');
                    for j = 1:nz
                        if (j <= 0.5*nz) 
                            TP_3D(1,:,j)     =  TP1_3D(1,:,1); % surface topography
                        else           
                            TP_3D(1,:,j)     =  TP1_3D(1,:,2); % cmb topography
                        end
                    end             
                end
               
                if WriteStress
                    [nnb,X_3D, Y_3D, Z_3D,S_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'stress');
                end
                
                if WriteEdot
                    [nnb,X_3D, Y_3D, Z_3D,E_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'strain rate');
                end
                
                if WriteDamage
                    [nnb,X_3D, Y_3D, Z_3D,D_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'damage');
                end
                
                %WriteStag3D_LegacyVTK_Cartesian;   %-> Writes VTK legacy files
                WriteStag3D_VTK_Cartesian_290611;          %-> Writes VTK-XML files

                % Store the filename and time in a structure, which is later
                % used to create a pvd file, that contains all timesteps and
                % times of the files
               
                FileNames{num} = {time,fname_number,fname_vtk};

          end
%%           
          case 'Spherical'
            % SPHERICAL GRID
            
            % Read temperature information
            [nnb, X_3D, Y_3D, Z_3D,T_3D, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'temperature');

             if nnb==-999
                % The file does not exist, and we have finished processing all
                % data

                % Create a PVD file, which contains all data including
                % time-information
                Create_PVD_file(FileNames,fname,directory);

                error(['finished processing all files in directory'])

             else

                if WriteVelocity
                    % Read pressure & velocity information
                    [nnb,X_3D, Y_3D, Z_3D, VX_3D, VY_3D, VZ_3D, P_3D, time] = ReadStag3D_LeaV(directory, fname, fname_number,  'velocity');
                end

                if WriteViscosity
                    % Read viscosity information
                    [nnb,X_3D, Y_3D, Z_3D,ETA_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'viscosity');
                end

                if WriteComposition
                    % Read viscosity information
                    [nnb,X_3D, Y_3D, Z_3D,C_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'composition');
                end

                if WriteStress
                    % Read viscosity information
                    [nnb,X_3D, Y_3D, Z_3D,S_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'stress');
                end

                % Transform coordinates for Yin grid
                latmean = (max(max(max(X_3D))) + min(min(min(X_3D)))) / 2;       % spherical patch is centred at equator
                lonmean = (max(max(max(Y_3D))) + min(min(min(Y_3D)))) / 2;       
                R    = Z_3D + rcmb;  lat = latmean - X_3D;  lon = Y_3D - lonmean;
                X_3D = R.*cos(lat).*cos(lon);
                Y_3D = R.*cos(lat).*sin(lon);
                Z_3D = R.*sin(lat);
                
                % Transform velocities if they are needed
                if WriteVelocity
                    Vtheta = VX_3D; Vphi = VY_3D; Vr = VZ_3D;
                    VX_3D =  Vtheta.*sin(lat).*cos(lon) - Vphi.*sin(lon) + Vr.*cos(lat).*cos(lon);
                    VY_3D =  Vtheta.*sin(lat).*sin(lon) + Vphi.*cos(lon) + Vr.*cos(lat).*sin(lon);
                    VZ_3D = -Vtheta.*cos(lat)                            + Vr.*sin(lat);
                end

                %WriteStag3D_LegacyVTK_Cartesian;   %-> Writes VTK legacy files
                WriteStag3D_VTK_Cartesian_290611;          %-> Writes VTK-XML files

                % Store the filename and time in a structure, which is later
                % used to create a pvd file, that contains all timesteps and
                % times of the files
                FileNames{num} = {time,fname_number,fname_vtk};

             end
%%         
         case 'Annulus'
            % SPHERICAL ANNULUS GRID
            % same as spherical but add extra points to avoid a gap in the
            % visualisation

            % Read temperature information
            if (size_reduction)
                [nnb, X_3D, Y_3D, Z_3D,T_3D, time, rcmb]     = ReadStag3D2_reduce_size_LB(directory, fname, fname_number,  'temperature',size_reduction,nr_save,nrstart);
             else
                [nnb, X_3D, Y_3D, Z_3D,T_3D, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'temperature');
            end
             
            %[nnb, X_3D, Y_3D, Z_3D,T_3D, time, rcmb]     = ReadStag3D_LeaV(directory, fname, fname_number,  'temperature');
            ny = size(Y_3D,2);
            nz = size(Z_3D,3);
           
            
            % add one row for nicer visualization
            T_3D(1,end+1,:) = T_3D(1,1,:);

             if fname_number>number_end
                % The file does not exist, and we have finished processing all
                % data

                % Create a PVD file, which contains all data including
                % time-information
                Create_PVD_file(FileNames,fname,directory);

                error(['finished processing all files in directory'])

             else
                

               if WriteVelocity
                    % Read pressure & velocity information
                    [nnb,X_3D, Y_3D, Z_3D, VX_3D, VY_3D, VZ_3D, P_3D, time] = ReadStag3D_LeaV(directory, fname, fname_number,  'velocity');
                    VX_3D(1,end+1,:) = VX_3D(1,1,:);
                    VY_3D(1,end+1,:) = VY_3D(1,1,:);
                    VZ_3D(1,end+1,:) = VZ_3D(1,1,:);
                     P_3D(1,end+1,:) =  P_3D(1,1,:);
               end

                if WriteViscosity
                    % Read viscosity information
                    [nnb,X_3D, Y_3D, Z_3D,ETA_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'viscosity');
                    ETA_3D(1,end+1,:) = ETA_3D(1,1,:);
                end

                if WriteComposition
                    % Read viscosity information
                    [nnb,X_3D, Y_3D, Z_3D,C_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'composition');
                    C_3D(1,end+1,:) = C_3D(1,1,:);
                end
                
                if WriteContRoot
                    [nnb,X_3D, Y_3D, Z_3D,C1_3D, time]  = ReadStag3D_LeaV(directory, fname, fname_number,  'cont root');
                    C1_3D(1,end+1,:) = C1_3D(1,1,:);
                end
               
                if WriteContCrust
                    [nnb,X_3D, Y_3D, Z_3D,C2_3D, time]  = ReadStag3D_LeaV(directory, fname, fname_number,  'cont crust');
                    C2_3D(1,end+1,:) = C2_3D(1,1,:);
                end
                
                if WriteContBelt
                    [nnb,X_3D, Y_3D, Z_3D,C3_3D, time]  = ReadStag3D_LeaV(directory, fname, fname_number,  'cont belt');
                    C3_3D(1,end+1,:) = C3_3D(1,1,:);
                end
                
                if WriteContSuture
                    [nnb,X_3D, Y_3D, Z_3D,C4_3D, time]  = ReadStag3D_LeaV(directory, fname, fname_number,  'cont suture');
                    C4_3D(1,end+1,:) = C4_3D(1,1,:);
                end
                
                if WriteContPlot
                    [nnb,X_3D, Y_3D, Z_3D, CPL_3D, time]  = ReadStag3D_LeaV(directory, fname, fname_number,  'cont plot');
                    CPL_3D(1,end+1,:) = CPL_3D(1,1,:);
                end
                
                if WriteContinent
                    [nnb,X_3D, Y_3D, Z_3D,NRC_3D, time]  = ReadStag3D_LeaV(directory, fname, fname_number,  'continent');
                    NRC_3D(1,end+1,:) = NRC_3D(1,1,:);
                end
                
                if WriteAge
                    [nnb,X_3D, Y_3D, Z_3D,AGE_3D, time]  = ReadStag3D_LeaV(directory, fname, fname_number,  'age');
                    AGE_3D(1,end+1,:) = AGE_3D(1,1,:);
                end
                
                if WriteTopography 
                    % Read topography information
                    % CAUTION: 2D field - dimensions doesn't match with the
                    % other fields --> add zeros
                    TP_3D = zeros(1,ny,nz);
                    [nnb,X_3D, Y_3D, Z_3D,TP1_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'topography');
                    for j = 1:nz
                        if (j <= 0.5*nz) 
                            TP_3D(1,:,j)     =  TP1_3D(1,:,1); % surface topography
                        else           
                            TP_3D(1,:,j)     =  TP1_3D(1,:,2); % cmb topography
                        end
                    end
                    TP_3D(1,end+1,:) = TP_3D(1,1,:);        
                end
                
                if WriteTopoSG
                    % Read topography(selfgrav) information
                    % CAUTION: (N-1)DIM field - dimensions doesn't match with the
                    % other fields --> add zeros
                    TPSG_3D = zeros(1,ny,nz);
                    [nnb,X_3D, Y_3D, Z_3D,TPSG1_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'topography self-grav');
                    for j = 1:nz
                        if (j <= 0.5*nz) 
                            TPSG_3D(1,:,j)     =  TPSG1_3D(1,:,1); % surface topography
                        else           
                            TPSG_3D(1,:,j)     =  TPSG1_3D(1,:,2); % cmb topography
                        end
                    end
                    TPSG_3D(1,end+1,:) = TPSG_3D(1,1,:);        
                end
                
                if WriteGeoid
                    % Read geoid information
                    % CAUTION: DIM(N-1) field - dimensions doesn't match with the
                    % other fields --> add zeros
                    G_3D = zeros(1,ny,nz);
                    [nnb,X_3D, Y_3D, Z_3D,G1_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'geoid');
                    for j = 1:nz
                        if (j <= 0.5*nz) 
                            G_3D(1,:,j)     =  G1_3D(1,:,1); % surface topography
                        else           
                            G_3D(1,:,j)     =  G1_3D(1,:,2); % cmb topography
                        end
                    end
                    G_3D(1,end+1,:) = G_3D(1,1,:);
             
                end

                if WriteStress
                    % Read stress information
                    [nnb,X_3D, Y_3D, Z_3D,S_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'stress');
                    S_3D(1,end+1,:) = S_3D(1,1,:);
                end
                
                if WriteEdot
                    % Read strainrate information
                    [nnb,X_3D, Y_3D, Z_3D,E_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'strain rate');
                    E_3D(1,end+1,:) = E_3D(1,1,:);
                end
                
                if WriteDamage
                    % Read damage information
                    [nnb,X_3D, Y_3D, Z_3D,D_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'damage');
                    D_3D(1,end+1,:) = D_3D(1,1,:);
                end
                
                if WriteAge
                    % Read viscosity information
                    [nnb,X_3D, Y_3D, Z_3D,AGE_3D, time]   = ReadStag3D_LeaV(directory, fname, fname_number,  'age');
                    AGE_3D(1,end+1,:) = AGE_3D(1,1,:);
                end
                
                Y_3D(1,end+1,:) = Y_3D(1,1,:) + 2*pi;
                X_3D(1,end+1,:) = X_3D(1,1,:);
                Z_3D(1,end+1,:) = Z_3D(1,1,:);
                
               
                
                % Transform coordinates for Yin grid
                R    = Z_3D + rcmb;  lat = zeros(size(X_3D));  lon = Y_3D - pi;
                X_3D = R.*cos(lat).*cos(lon);
                Y_3D = R.*cos(lat).*sin(lon);
                Z_3D = R.*sin(lat);
                
                % Transform velocities if they are needed
                if WriteVelocity
                    Vtheta = VX_3D; Vphi = VY_3D; Vr = VZ_3D;
                    VX_3D =  Vtheta.*sin(lat).*cos(lon) - Vphi.*sin(lon) + Vr.*cos(lat).*cos(lon);
                    VY_3D =  Vtheta.*sin(lat).*sin(lon) + Vphi.*cos(lon) + Vr.*cos(lat).*sin(lon);
                    VZ_3D = -Vtheta.*cos(lat)                            + Vr.*sin(lat);
                end

                %WriteStag3D_LegacyVTK_Cartesian;   %-> Writes VTK legacy files
                WriteStag3D_VTK_Cartesian_290611;          %-> Writes VTK-XML files

                % Store the filename and time in a structure, which is later
                % used to create a pvd file, that contains all timesteps and
                % times of the files
                FileNames{num} = {time,fname_number,fname_vtk};
             end
%%
         otherwise

             error('Unknown gridtype')
     end
 
     num = num+1;

end

Create_PVD_file(FileNames,fname,directory);
disp('Finished')











