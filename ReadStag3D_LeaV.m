function [varargout] = ReadStag3D_Lea(directory, fname_input, fname_number,Type)% Reads a STAG3D output file and transforms it into a MATLAB file
%
% Tested for 3D cartesian
%
%
% Syntax:
%     For scalar fields:
%
%       [X_3D, Y_3D, Z_3D, DATA_3D] = ReadStag3D(fname_input, fname_number, {'viscosity','temperature'});
%
%   For vector fields:
%
%       [X_3D, Y_3D, Z_3D, VX_3D, VY_3D, VZ_3D, P_3D] = ReadStag3D(fname_input, fname_number, {'viscosity','temperature'});
%
%
% Boris Kaus
%
% Modifications by Paul Tackley August 2010:
%  - extra fields: topography, topography self-grav, crustal
%  thickness, age, strain rate, geoid
%  - unimplemented 2Dfield stuff deleted: 3D routine works for 2D fields
%  - velocities & p now read into 4D arrays: yy split is done on output
%  - correction to the storage order of v* and p
%  - checks magic whether #components matches nval
%  - compatibility with magic=9 storing extra x & y points for vectors
%  - added missing scaling factor read for vectors 
%
%
% Adapted by Lea Bello and Marie Bocher to be used in Error_Temp (for forecasting time calculation)
%

start_dir = pwd;
cd(directory);

FileFormat = 'n';            % native - default
%FileFormat = 'l';           % Little Endian
%FileFormat = 'b';           % Big    Endian

if fname_number<10000
    number_string = num2str(10000+fname_number);
    number_string(1)='0';
else
    number_string = num2str(10000+fname_number);
end

switch Type
    case 'velocity'
        fname       = [fname_input,'_vp',number_string];
        scalardata  = false;
    case 'residue'  
        fname       = [fname_input,'_rs',number_string];
        scalardata  = false;
    case 'temperature'
        fname       = [fname_input,'_t',number_string];
        scalardata  = true;
    case 'viscosity'
        fname       = [fname_input,'_eta',number_string];
        scalardata  = true;
    case 'composition'
        fname       = [fname_input,'_comp',number_string];
        scalardata  = true;
    case 'cont root'
        fname       = [fname_input,'_c',number_string];
        scalardata  = true;
    case 'cont crust'
        fname       = [fname_input,'_cc',number_string];
        scalardata  = true;
    case 'cont belt'
        fname       = [fname_input,'_ccc',number_string];
        scalardata  = true;
    case 'cont suture'
        fname       = [fname_input,'_cccc',number_string];
        scalardata  = true;
    case 'cont plot'
        fname       = [fname_input,'_cpl',number_string];
        scalardata  = true;
    case 'continent'
        fname       = [fname_input,'_nrc',number_string];
        scalardata  = true;
    case 'melt fraction'
        fname       = [fname_input,'_f',number_string];
        scalardata  = true;
    case 'topography'
        fname       = [fname_input,'_cs',number_string];
        scalardata  = true;
    case 'topography self-grav'
        fname       = [fname_input,'_csg',number_string];
        scalardata  = true;
    case 'crustal thickness'
        fname       = [fname_input,'_cr',number_string];
        scalardata  = true;
    case 'age'
        fname       = [fname_input,'_age',number_string];
        scalardata  = true;
    case 'strain rate'
        fname       = [fname_input,'_ed',number_string];
        scalardata  = true;
    case 'geoid'
        fname       = [fname_input,'_g',number_string];
        scalardata  = true;
    case 'stress'
        fname       = [fname_input,'_str',number_string];
        scalardata  = true; 
    case 'damage'
        fname       = [fname_input,'_d',number_string];
        scalardata  = true; 
    case 'heat flux'
        fname       = [fname_input,'_hf',number_string];
        scalardata  = true;
    otherwise
        error('Unknown property')
end
if ~exist(fname)
    % The file does not exist and we should stop processing data

    cd(start_dir);
    for i=1:10
        varargout{i}    = -999;
    end

    return
end

if scalardata
    nval    =   1;      % temperature has only one value
else
    nval    =   4;      % assumed that we have a velocity-pressure file
end

fid         =   fopen(fname,'r',FileFormat);               % Open File

%==========================================================================
% READ HEADER
%==========================================================================
magic       = fread(fid,1,'int32');         % Version
if (magic<100 && nval>1) || (magic>300 && nval==1) % check #components
    error('wrong number of components in field')
end
magic       = mod(magic,100);
if magic>=9 && nval==4
    xyp = 1;     % extra ghost point in x & y direction
else
    xyp = 0;
end

nxtot       = fread(fid,1,'int32');         % nx total
nytot       = fread(fid,1,'int32');         % ny total
nztot       = fread(fid,1,'int32');         % nz total
nblocks     = fread(fid,1,'int32');         % # of blocks, 2 for yinyang
Aspect      = fread(fid,2,'single');        % Aspect ratio
nnx         = fread(fid,1,'int32');         % Number of parallel subdomains
nny         = fread(fid,1,'int32');         %  in the x,y,z and b directions
nnz         = fread(fid,1,'int32');         %
nnb         = fread(fid,1,'int32');         %

nz2         = nztot*2 + 1;
zg          = fread(fid,nz2,'single');      % z-coordinates

% compute nx, ny, nz and nb PER CPU
nx          =   nxtot/nnx;
ny          =   nytot/nny;
nz          =   nztot/nnz;
nb          =   nblocks/nnb;
npi         =   (nx+xyp)*(ny+xyp)*nz*nb*nval;      % the number of values per 'read' block


rcmb        = fread(fid,1,'single');
istep       = fread(fid,1,'int32');
time        = fread(fid,1,'single');
erupta_total= fread(fid,1,'single');
botT_val    = fread(fid,1,'single');


x           = fread(fid,nxtot,'single');      % x-coordinates
y           = fread(fid,nytot,'single');      % y-coordinates
z           = fread(fid,nztot,'single');      % z-coordinates

% read the parallel blocks
if scalardata
    DATA_3D = zeros(nxtot,nytot,nztot);
else
    scalefac= fread(fid,1,'single');            % scale factor
    VX_3D   = zeros(nxtot,nytot,nztot);         %   Vx
    VY_3D   = zeros(nxtot,nytot,nztot);         %   Vy
    VZ_3D   = zeros(nxtot,nytot,nztot);         %   Vz
    P_3D    = zeros(nxtot,nytot,nztot);         %   Pressure
end

for ibc=1:nnb       % loop over parallel subdomains
    for izc =1:nnz
        for iyc =1:nny
            for ixc =1:nnx
                data_CPU    =   fread(fid,npi,'single');      % read the data for this CPU

                % Create a 3D matrix from these data
                if scalardata
                    data_CPU_3D = reshape(data_CPU, [nx ny nz nb]) ;                  
                else
                    data_CPU_3D = reshape(data_CPU*scalefac, [nval nx+xyp ny+xyp nz nb]);
                end

                % Add local 3D matrix to global matrix
                if scalardata
                    % Scalar data
                    DATA_3D((ixc-1)*nx + (1:nx), (iyc-1)*ny + (1:ny), (izc-1)*nz + (1:nz),(ibc-1)*nb + (1:nb)) = data_CPU_3D;
                    
                else
                    % velocity-pressure data
                    VX_3D((ixc-1)*nx + (1:nx), (iyc-1)*ny +  (1:ny), (izc-1)*nz + (1:nz), (ibc-1)*nb + (1:nb)) = squeeze(data_CPU_3D(1,1:nx,1:ny,:,:));
                    VY_3D((ixc-1)*nx + (1:nx), (iyc-1)*ny +  (1:ny), (izc-1)*nz + (1:nz), (ibc-1)*nb + (1:nb)) = squeeze(data_CPU_3D(2,1:nx,1:ny,:,:));
                    VZ_3D((ixc-1)*nx + (1:nx), (iyc-1)*ny +  (1:ny), (izc-1)*nz + (1:nz), (ibc-1)*nb + (1:nb)) = squeeze(data_CPU_3D(3,1:nx,1:ny,:,:));
                    P_3D( (ixc-1)*nx + (1:nx), (iyc-1)*ny +  (1:ny), (izc-1)*nz + (1:nz), (ibc-1)*nb + (1:nb)) = squeeze(data_CPU_3D(4,1:nx,1:ny,:,:));
                
                    RX_3D((ixc-1)*nx + (1:nx), (iyc-1)*ny +  (1:ny), (izc-1)*nz + (1:nz), (ibc-1)*nb + (1:nb)) = squeeze(data_CPU_3D(1,1:nx,1:ny,:,:));
                    RY_3D((ixc-1)*nx + (1:nx), (iyc-1)*ny +  (1:ny), (izc-1)*nz + (1:nz), (ibc-1)*nb + (1:nb)) = squeeze(data_CPU_3D(2,1:nx,1:ny,:,:));
                    RZ_3D((ixc-1)*nx + (1:nx), (iyc-1)*ny +  (1:ny), (izc-1)*nz + (1:nz), (ibc-1)*nb + (1:nb)) = squeeze(data_CPU_3D(3,1:nx,1:ny,:,:));
                    RP_3D( (ixc-1)*nx + (1:nx), (iyc-1)*ny +  (1:ny), (izc-1)*nz + (1:nz), (ibc-1)*nb + (1:nb)) = squeeze(data_CPU_3D(4,1:nx,1:ny,:,:));
                end


            end
        end
    end
end

fclose(fid);                                % close file


[Y_3D, X_3D, Z_3D]  = meshgrid(y,x,z);


% Transform coordinates for Yin & Yang grids

RYin  = Z_3D+rcmb;
ThYin = pi/4-X_3D;
PhiYin = Y_3D-3*pi/4;

XYin=RYin.*cos(ThYin).*cos(PhiYin);
YYin=RYin.*cos(ThYin).*sin(PhiYin);
ZYin = RYin.*sin(ThYin);

XYang=-XYin;
YYang=ZYin;
ZYang=YYin;

RedFlagYin=zeros(size(XYin));
RedFlagYang=zeros(size(XYang));

for kray=1:nztot
    minYYang=min(min(YYang(:,:,kray)));
    maxYYang=max(max(YYang(:,:,kray)));
    YinDel=XYin(:,:,kray)<=0&YYin(:,:,kray)<=maxYYang&YYin(:,:,kray)>=minYYang;
    RedFlagYin(:,:,kray)=YinDel;
    YangDel=XYang(:,:,kray)>=0&ZYang(:,:,kray)<=max(max(ZYin(:,:,kray)))&ZYang(:,:,kray)>=min(min(ZYin(:,:,kray)));
    RedFlagYang(:,:,kray)=YangDel;
end


% erase overlapping points

RedFlag=[RedFlagYin RedFlagYang];
RYinYang=[RYin RYin];
CoordX=[XYin XYang];
CoordY=[YYin YYang];
CoordZ=[ZYin ZYang];

%% volume
VolYin=zeros(nxtot,nytot,nztot);
Vol=zeros(nxtot,2*nytot,nztot);
for i=1:nxtot
for k=1:nztot
    dz(k)= zg(2*k+1)-zg(2*k-1);
    VolYin(i,:,k)= dz(k) .* sin(pi/4+X_3D(i,:,k)).*(pi/(2*nxtot))^2 .* (rcmb+zg(2*k)).^2;
end
end
Vol=[VolYin VolYin]; %same evolution of the volume on yin and yang, decrasing with colatitude=pi/4+X_3D(i,:,k)
%Vol=Vol(~logical(RedFlag));       
        
        
%% prepare output data

varargout{1} = Vol; % old version: varargout{1} = nblocks;
varargout{2} =CoordX; % ThYin;
varargout{3} =CoordY; %PhiYin;
varargout{4} =CoordZ;  %CoordZ;

if nblocks==1
    % no ying-yang
    switch Type
        case 'velocity'
            varargout{5} = VX_3D;
            varargout{6} = VY_3D;
            varargout{7} = VZ_3D;
            varargout{8} = P_3D;
            varargout{9} = time;
            varargout{10}= rcmb;
        otherwise
            varargout{5} = DATA_3D;
            varargout{6} = time;
            varargout{7} = rcmb;
    end
else
    % ying-yang grid
    switch Type
        case 'velocity'
            varargout{5} = VX_3D(:,:,:,1);
            varargout{6} = VY_3D(:,:,:,1);
            varargout{7} = VZ_3D(:,:,:,1);
            varargout{8} = RedFlag; %previous version: P_3D (:,:,:,1);

            varargout{9}  = VX_3D(:,:,:,2);
            varargout{10} = VY_3D(:,:,:,2);
            varargout{11} = VZ_3D(:,:,:,2);
            varargout{12} = CoordZ(:,:,:); %previous version: P_3D (:,:,:,2);

            varargout{13} = time;
            varargout{14} = rcmb;
            varargout{15} = z;            
            
%% calcul of vsurf

    RedFlagsurf=[RedFlagYin(:,:,size(RedFlagYin,3)) RedFlagYang(:,:,size(RedFlagYang,3))];
    VXsurf = [VX_3D(:,:,size(VX_3D,3),1) VX_3D(:,:,size(VX_3D,3),2)];
    VYsurf = [VY_3D(:,:,size(VY_3D,3),1) VY_3D(:,:,size(VY_3D,3),2)];
    VZsurf = [VZ_3D(:,:,size(VZ_3D,3),1) VZ_3D(:,:,size(VZ_3D,3),2)];
    vsurf_tot=(VXsurf.^2 + VYsurf.^2);
    vsurf_tot_real=vsurf_tot(~logical(RedFlagsurf));

    %bottom velocity
    %RedFlagbot=[RedFlagYin(:,:,1) RedFlagYang(:,:,1)];
    %VXsurf = [VX_3D(:,:,1,1) VX_3D(:,:,1,2)];
    %VYsurf = [VY_3D(:,:,1,1) VY_3D(:,:,1,2)];
    %vsurf_tot=(VXsurf.^2 + VYsurf.^2);
    %vsurf_tot_real=vsurf_tot(~logical(RedFlagbot));
          
            varargout{16} = mean(mean(vsurf_tot_real)).^(1/2);  
                                    
% calcul of v

    VX = [VX_3D(:,:,:,1) VX_3D(:,:,:,2)];
    VY = [VY_3D(:,:,:,1) VY_3D(:,:,:,2)];
    VZ = [VZ_3D(:,:,:,1) VZ_3D(:,:,:,2)];
    v_tot=(VX.^2 + VY.^2 + VZ.^2);
    
             varargout{17} = v_tot;  %CAUTION: v_tot=velocity*Volume of cell
             varargout{18} = VXsurf;
             varargout{19} = VYsurf;
             varargout{20} = RedFlagsurf;
             
             
        case 'residue'
            varargout{5} = RX_3D(:,:,:,1);
            varargout{6} = RY_3D(:,:,:,1);
            varargout{7} = RZ_3D(:,:,:,1);
            varargout{8} = RP_3D (:,:,:,1);

            varargout{9}  = RX_3D(:,:,:,2);
            varargout{10} = RY_3D(:,:,:,2);
            varargout{11} = RZ_3D(:,:,:,2);
            varargout{12} = RP_3D (:,:,:,2);

            varargout{13} = time;
            varargout{14} = zg;
            
        otherwise  
            varargout{5} = DATA_3D(:,:,:,1);
            varargout{6} = DATA_3D(:,:,:,2);
            
            DATA_tot = [DATA_3D(:,:,:,1) DATA_3D(:,:,:,2)];
            varargout{7} = DATA_tot(~logical(RedFlag));
            varargout{8} = RedFlag;
            
            varargout{9}  = XYin;
            varargout{10} = YYin;
            varargout{11} = ZYin;
            varargout{12} = CoordZ(:,:,:); %previous version: P_3D (:,:,:,2);
            
            varargout{13} = time;
            varargout{14} = rcmb;
            varargout{15} = z;  
    end

end

cd(start_dir);



