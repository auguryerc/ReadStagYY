% Write STAG_3D format from matlab in XML-VTK format, using ASCII or BINARY file
% format, in case we are dealing with a YING-YANG GRID.
%
%
%
% Boris Kaus, 26 Feb. 2008
% 8.13.2008  Bugfix that caused compositional fields etc. not to be written
%            correctly to file

% Paul Tackley 3.08.2010: corrections to pressure and velocity writing

% The ying-yang grid consists of two blocks, with overlapping corners.
% The corners are removed, and the grid is triangulated. A single unstructured
% mesh is generated from this, that has VT_WEDGE elements. 

ASCII           = false;                       % ASCII or BINARY?

% Change to correct directory
start_dir       = pwd;
cd(directory);

%==========================================================================
% 1) Take the surface of the 2 grids, patch together and triangulate
%==========================================================================
xs1                 =   X_3D_1(:,:,end);
ys1                 =   Y_3D_1(:,:,end);
zs1                 =   Z_3D_1(:,:,end);

xs2                 =   X_3D_2(:,:,end);
ys2                 =   Y_3D_2(:,:,end);
zs2                 =   Z_3D_2(:,:,end);

r1                  =   sqrt(xs1.^2 + ys1.^2 + zs1.^2);
theta1              =   atan2(sqrt(xs1.^2 + ys1.^2),zs1);
phi1                =   atan2(ys1,xs1);


r2                  =   sqrt(xs2.^2 + ys2.^2 + zs2.^2);
theta2              =   atan2(sqrt(xs2.^2 + ys2.^2),zs2);
phi2                =   atan2(ys2,xs2);


% Cut off the corners from grid #1:
%radius_critical     =   ((phi2(1,1)+pi/2)^2+(theta2(1,1)-pi/2)^2);
%radius_grid1        =   ((phi1+pi/2).^2+(theta1-pi/2).^2);
%ind_corner1      	=   find(radius_grid1>=radius_critical & (phi1+pi/2)>0 & (phi1+pi/2)<pi/2);
%
% Cut off the corners from grid #1:
%radius_critical     =   ((phi2(end,end)-pi/2)^2+(theta2(end,end)-pi/2)^2);
%radius_grid2        =   ((phi1-pi/2).^2+(theta1-pi/2).^2);
%ind_corner2         =   find(radius_grid2>=radius_critical & (phi1-pi/2)<0 & (phi1-pi/2)>-pi/2 );

% Cut off the corners from grid #1, which seems to do #2 as well (PJT):
theta12             = acos(sin(theta1).*sin(phi1));  % theta coords of grid 1 in coord system of grid 2
ind_corner          = find( (theta12>pi/4 & phi1>pi/2) | (theta12<3*pi/4 & phi1<-pi/2 ) );


% Form indices of remaining 
ind                           = 1:prod(size(phi1));
ind([ind_corner])=[];
%ind([ind_corner1 ind_corner2])=[];

% Create an array with unique r, theta & phi values
numberYin = ones(size(r1));
numberYin(find(r1==r1))           =   find(r1==r1);
numberYang = ones(size(r2));
numberYang(find(r2==r2))          =   find(r2==r2)+ max(numberYin(:));

%Closed Yin grid, w/o corners and completely closed
r                   =   [r1(ind),      ]';
theta               =   [theta1(ind)   ]';
phi                 =   [phi1(ind),    ]';
NumYin              =   [numberYin(ind)]';


% x,y,z coordinates of complete grid:
x                   =   [xs1(ind) -xs1(ind)];
y                   =   [ys1(ind)  zs1(ind)];
z                   =   [zs1(ind)  ys1(ind)];
NumYang             =   [NumYin + max(numberYin(:)) ];
Number              =   [NumYin(:);    NumYang(:)  ]';

%   remove redundant coordinates
tri                 =   convhulln([x(:), y(:),z(:)]);       % simple way to grid it
triYingYang         =   Number(tri);

x                   =   [xs1(:); xs2(:)];
y                   =   [ys1(:); ys2(:)];
z                   =   [zs1(:); zs2(:)];

%==========================================================================
% triYingYang now contains the numbers of all triangles 
%==========================================================================

%==========================================================================
% 2) Create a 3D grid with tetrahedron elements
%==========================================================================

% Number all gridpoints we have
NUMBER_1                 = ones(size(X_3D_1));
NUMBER_2                 = ones(size(X_3D_2));
NUMBER_1(find(NUMBER_1)) = find(NUMBER_1);
NUMBER_2(find(NUMBER_2)) = find(NUMBER_2) + max(NUMBER_1(:));


% Make a loop over all levels
ElementNumbers      = [];
for iz=1:size(X_3D_2,3)-1

    num_upper1   	= NUMBER_1(:,:,iz+1);
    num_upper2      = NUMBER_2(:,:,iz+1);
    num_upper       = [num_upper1(:); num_upper2(:)];

    num_lower1      = NUMBER_1(:,:,iz);
    num_lower2      = NUMBER_2(:,:,iz);
    num_lower       = [num_lower1(:); num_lower2(:)];
    ElementNumbers  = [ElementNumbers; num_upper(triYingYang), num_lower(triYingYang)];
end


%--------------------------------------------------------------------------
% Convert data into correct vector format
%--------------------------------------------------------------------------
Points                  = zeros(max(NUMBER_2(:)),3);
Points(NUMBER_1(:),1)       = X_3D_1(:);
Points(NUMBER_2(:),1)       = X_3D_2(:);
Points(NUMBER_1(:),2)       = Y_3D_1(:);
Points(NUMBER_2(:),2)       = Y_3D_2(:);
Points(NUMBER_1(:),3)       = Z_3D_1(:);
Points(NUMBER_2(:),3)       = Z_3D_2(:);

Radius                      = sqrt(sum(Points.^2,2));


Temperature              = zeros(max(NUMBER_2(:)),1);
Temperature(NUMBER_1(:)) = T_3D_1(:);
Temperature(NUMBER_2(:)) = T_3D_2(:);

if WriteVelocity
    Velocity_vec                = zeros(max(NUMBER_2(:)),3);
    Velocity_vec(NUMBER_1(:),1) = VX_3D_1(:);
    Velocity_vec(NUMBER_2(:),1) = VX_3D_2(:);
    Velocity_vec(NUMBER_1(:),2) = VY_3D_1(:);
    Velocity_vec(NUMBER_2(:),2) = VY_3D_2(:);
    Velocity_vec(NUMBER_1(:),3) = VZ_3D_1(:);
    Velocity_vec(NUMBER_2(:),3) = VZ_3D_2(:);
    
    Pressure_vec    = zeros(max(NUMBER_2(:)),1);
    Pressure_vec(NUMBER_1(:),1) = P_3D_1(:);
    Pressure_vec(NUMBER_2(:),1) = P_3D_2(:);
end

if WriteResidue
    Residue_vec                = zeros(max(NUMBER_2(:)),3);
    Residue_vec(NUMBER_1(:),1) = RX_3D_1(:);
    Residue_vec(NUMBER_2(:),1) = RX_3D_2(:);
    Residue_vec(NUMBER_1(:),2) = RY_3D_1(:);
    Residue_vec(NUMBER_2(:),2) = RY_3D_2(:);
    Residue_vec(NUMBER_1(:),3) = RZ_3D_1(:);
    Residue_vec(NUMBER_2(:),3) = RZ_3D_2(:);
   
    ResidueP_vec    = zeros(max(NUMBER_2(:)),1);
    ResidueP_vec(NUMBER_1(:),1) = RP_3D_1(:);
    ResidueP_vec(NUMBER_2(:),1) = RP_3D_2(:);
end

if WriteResidualT
    resT_vec                = zeros(max(NUMBER_2(:)),1);
    resT_vec(NUMBER_1(:),1) = resT_3D_1(:);
    resT_vec(NUMBER_2(:),1) = resT_3D_2(:);
end

if WriteViscosity
    Viscosity_vec                = zeros(max(NUMBER_2(:)),1);
    Viscosity_vec(NUMBER_1(:),1) = ETA_3D_1(:);
    Viscosity_vec(NUMBER_2(:),1) = ETA_3D_2(:);
end

if WriteComposition
    Composition_vec                = zeros(max(NUMBER_2(:)),1);
    Composition_vec(NUMBER_1(:),1) = C_3D_1(:);
    Composition_vec(NUMBER_2(:),1) = C_3D_2(:);
end

if WriteContRoot
    Composition1_vec                = zeros(max(NUMBER_2(:)),1);
    Composition1_vec(NUMBER_1(:),1) = C1_3D_1(:);
    Composition1_vec(NUMBER_2(:),1) = C1_3D_2(:);
end

if WriteContCrust
    Composition2_vec                = zeros(max(NUMBER_2(:)),1);
    Composition2_vec(NUMBER_1(:),1) = C2_3D_1(:);
    Composition2_vec(NUMBER_2(:),1) = C2_3D_2(:);
end

if WriteContBelt
    Composition3_vec                = zeros(max(NUMBER_2(:)),1);
    Composition3_vec(NUMBER_1(:),1) = C3_3D_1(:);
    Composition3_vec(NUMBER_2(:),1) = C3_3D_2(:);
end

if WriteContSuture
    Composition4_vec                = zeros(max(NUMBER_2(:)),1);
    Composition4_vec(NUMBER_1(:),1) = C4_3D_1(:);
    Composition4_vec(NUMBER_2(:),1) = C4_3D_2(:);
end

if WriteTopography
    Topography_vec                = zeros(max(NUMBER_2(:)),1);
    Topography_vec(NUMBER_1(:),1) = TP_3D_1(:);
    Topography_vec(NUMBER_2(:),1) = TP_3D_2(:);
end

if WriteTopoSG
    TopoSG_vec                = zeros(max(NUMBER_2(:)),1);
    TopoSG_vec(NUMBER_1(:),1) = TPSG_3D_1(:);
    TopoSG_vec(NUMBER_2(:),1) = TPSG_3D_2(:);
end

if WriteGeoid
    Geoid_vec                = zeros(max(NUMBER_2(:)),1);
    Geoid_vec(NUMBER_1(:),1) = G_3D_1(:);
    Geoid_vec(NUMBER_2(:),1) = G_3D_2(:);
end

if WriteHeatflux
    HF_vec                = zeros(max(NUMBER_2(:)),1);
    HF_vec(NUMBER_1(:),1) = HF_3D_1(:);
    HF_vec(NUMBER_2(:),1) = HF_3D_2(:);
end

if WriteStress
    Stress_vec                = zeros(max(NUMBER_2(:)),1);
    Stress_vec(NUMBER_1(:),1) = S_3D_1(:);
    Stress_vec(NUMBER_2(:),1) = S_3D_2(:);
end

if WriteEdot
    Strainrate_vec                = zeros(max(NUMBER_2(:)),1);
    Strainrate_vec(NUMBER_1(:),1) = E_3D_1(:);
    Strainrate_vec(NUMBER_2(:),1) = E_3D_2(:);
end

if WriteDamage
    Damage_vec                = zeros(max(NUMBER_2(:)),1);
    Damage_vec(NUMBER_1(:),1) = D_3D_1(:);
    Damage_vec(NUMBER_2(:),1) = D_3D_2(:);
end

if WriteContinent
    Continent_vec                = zeros(max(NUMBER_2(:)),1);
    Continent_vec(NUMBER_1(:),1) = NRC_3D_1(:);
    Continent_vec(NUMBER_2(:),1) = NRC_3D_2(:);
end

if WriteContPlot
    ContPlot_vec                = zeros(max(NUMBER_2(:)),1);
    ContPlot_vec(NUMBER_1(:),1) = CPL_3D_1(:);
    ContPlot_vec(NUMBER_2(:),1) = CPL_3D_2(:);
end

if WriteAge
    Age_vec                = zeros(max(NUMBER_2(:)),1);
    Age_vec(NUMBER_1(:),1) = AGE_3D_1(:)*75.*1100.;% multipli? ici par le transit time ; A CHANGER !!
    Age_vec(NUMBER_2(:),1) = AGE_3D_2(:)*75*1100.;%
end

if WriteDiv
    Div_vec                = zeros(max(NUMBER_2(:)),1);
    Div_vec(NUMBER_1(:),1) = DIV_3D_1(:);% multipli? ici par le transit time ; A CHANGER !!
    Div_vec(NUMBER_2(:),1) = DIV_3D_2(:);%
end

if WriteVor
    Vor_vec                = zeros(max(NUMBER_2(:)),1);
    Vor_vec(NUMBER_1(:),1) = VOR_3D_1(:);% multipli? ici par le transit time ; A CHANGER !!
    Vor_vec(NUMBER_2(:),1) = VOR_3D_2(:);%
end
%==========================================================================
% 3) Write VTK file (unstructured mesh)
%==========================================================================
ElementNumbers  = ElementNumbers-1;      % VTK is zero-based 


%==========================================================================
% Definitions and initialization
sizeof_Float32  =   4;      
sizeof_Float64  =   4;     
sizeof_UInt32   =   4; 
sizeof_UInt8    =   1; 

Offset          =   0;      % Initial offset

%--------------------------------------------------------------------------
% Write the header for a structured grid:
%--------------------------------------------------------------------------
fname_vtk       = [fname,'_',num2str(1000000+fname_number),'.vtu'];
fid             = fopen(fname_vtk,'w','b');           % note the 'b': not doing BigEndian does not work with MATLAB!
fprintf(fid,'<?xml version="1.0"?> \n');
fprintf(fid,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian" >\n');
fprintf(fid,'  <UnstructuredGrid>\n');
fprintf(fid,'    <Piece NumberOfPoints="%i"  NumberOfCells="%i">\n', int32(size(Points,1)), int32(size(ElementNumbers,1)));



%--------------------------------------------------------------------------
% Add point-wise data
%--------------------------------------------------------------------------
fprintf(fid,'    <PointData Scalars="T" Vectors="Velocity"  >\n');

% TEMPERATURE -----------
if ASCII
    % ASCII:
    fprintf(fid,'      <DataArray type="Float32" Name="T" format="ascii">\n');
    for i=1:length(Temperature)
        fprintf(fid,'        %g \n',single(Temperature(i)));
    end
else
    % BINARY:
    fprintf(fid,'      <DataArray type="Float32" Name="T" format="appended" offset="%i">\n', int32(Offset));
    Offset = Offset + length(Temperature(:))*sizeof_Float32 + 1*sizeof_UInt32;
end
fprintf(fid,'      </DataArray>\n');
% -----------------------

if WriteResidualT
    % residual Temperature---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="ResidualT" format="ascii">\n');
        for i=1:length(resT_vec)
            fprintf(fid,'        %g \n',single(resT_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="ResidualT" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(resT_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteViscosity
    % VISCOSITY---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Viscosity" format="ascii">\n');
        for i=1:length(Viscosity_vec)
            fprintf(fid,'        %g \n',single(Viscosity_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Viscosity" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Viscosity_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteComposition
    % COMPOSITION---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Composition" format="ascii">\n');
        for i=1:length(Composition_vec)
            fprintf(fid,'        %g \n',single(Composition_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Composition" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Composition_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end


if WriteContRoot
    % COMPOSITION1---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Root" format="ascii">\n');
        for i=1:length(Composition1_vec)
            fprintf(fid,'        %g \n',single(Composition1_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Root" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Composition1_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteContCrust
    % COMPOSITION2---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Crust" format="ascii">\n');
        for i=1:length(Composition2_vec)
            fprintf(fid,'        %g \n',single(Composition2_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Crust" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Composition2_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteContBelt
    % COMPOSITION3---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Belt" format="ascii">\n');
        for i=1:length(Composition3_vec)
            fprintf(fid,'        %g \n',single(Composition3_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Belt" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Composition3_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteContSuture
    % COMPOSITION3---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Suture" format="ascii">\n');
        for i=1:length(Composition4_vec)
            fprintf(fid,'        %g \n',single(Composition4_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Suture" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Composition4_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteContPlot
    % Continent---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="cont plot" format="ascii">\n');
        for i=1:length(ContPlot_vec)
            fprintf(fid,'        %g \n',single(ContPlot_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="cont plot" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(ContPlot_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteContinent
    % Continent---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Continent" format="ascii">\n');
        for i=1:length(Continent_vec)
            fprintf(fid,'        %g \n',single(Continent_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Continent" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Continent_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteTopography
    % COMPOSITION3---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Topography" format="ascii">\n');
        for i=1:length(Topography_vec)
            fprintf(fid,'        %g \n',single(Topography_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Topography" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Topography_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteTopoSG
    % COMPOSITION3---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Topography self-grav" format="ascii">\n');
        for i=1:length(TopoSG_vec)
            fprintf(fid,'        %g \n',single(TopoSG_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Topography self-grav" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(TopoSG_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteGeoid
    % COMPOSITION3---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Geoid" format="ascii">\n');
        for i=1:length(Geoid_vec)
            fprintf(fid,'        %g \n',single(Geoid_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Geoid" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Geoid_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteHeatflux
    % COMPOSITION3---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="heat flux" format="ascii">\n');
        for i=1:length(HF_vec)
            fprintf(fid,'        %g \n',single(HF_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="heat flux" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(HF_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteStress
    % STRESS---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Stress" format="ascii">\n');
        for i=1:length(Stress_vec)
            fprintf(fid,'        %g \n',single(Stress_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Stress" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Stress_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end


if WriteEdot
    % STRAIN RATE---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Strain rate" format="ascii">\n');
        for i=1:length(Strainrate_vec)
            fprintf(fid,'        %g \n',single(Strainrate_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Strain rate" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Strainrate_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteDamage
    % Damage---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="DamA" format="ascii">\n');
        for i=1:length(Damage_vec)
            fprintf(fid,'        %g \n',single(Damage_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Damage" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Damage_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteVelocity
    
    % PRESSURE---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="P" format="ascii">\n');
        for i=1:length(Pressure_vec)
            fprintf(fid,'        %g \n',single(Pressure_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="P" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Pressure_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------

    % VELOCITY---------------  : VELOCITY IS A 3-component vector
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">\n');
        for i=1:length(T)
            fprintf(fid,'   %g %g %g \n',single(Velocity_vec(i,:)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Velocity_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
    
end

if WriteResidue
    
    % PRESSURE---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Res_Mag" format="ascii">\n');
        for i=1:length(ResidueP_vec)
            fprintf(fid,'        %g \n',single(ResidueP_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Res_Mag" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(ResidueP_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------

    % VELOCITY---------------  : VELOCITY IS A 3-component vector
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Residue" NumberOfComponents="3" format="ascii">\n');
        for i=1:length(T)
            fprintf(fid,'   %g %g %g \n',single(Residue_vec(i,:)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Residue" NumberOfComponents="3" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Residue_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
    
end

if WriteAge
    % Age---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Age" format="ascii">\n');
        for i=1:length(Age_vec)
            fprintf(fid,'        %g \n',single(Age_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Age" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Age_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteDiv
    % Age---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Divergence" format="ascii">\n');
        for i=1:length(Div_vec)
            fprintf(fid,'        %g \n',single(Div_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Divergence" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Div_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteVor
    % Age---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Vorticity" format="ascii">\n');
        for i=1:length(Vor_vec)
            fprintf(fid,'        %g \n',single(Vor_vec(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Vorticity" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Vor_vec(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

fprintf(fid,'    </PointData>\n');

%--------------------------------------------------------------------------
% Add coordinates of structured grid 
%--------------------------------------------------------------------------
fprintf(fid,'    <Points>\n');

% ASCII
if ASCII
    fprintf(fid,'      <DataArray type="Float32" Name="Array" NumberOfComponents="3" format="ascii">\n');
    for i=1:size(Points,1)
        fprintf(fid,'         %g %g %g \n',[Points(i,:)]);
    end
else
     fprintf(fid,'      <DataArray type="Float32" Name="Array" NumberOfComponents="3" format="appended" offset="%i" >\n',int32(Offset));
     Offset = Offset + length(Points(:))*sizeof_Float32 + 1*sizeof_UInt32;
end
fprintf(fid,'      </DataArray>\n');
fprintf(fid,'    </Points>\n');


%--------------------------------------------------------------------------
% Add CELLS data
%--------------------------------------------------------------------------
fprintf(fid,'    <Cells>\n');

% Connectivity -----------
if ASCII
    % ASCII:
    fprintf(fid,'      <DataArray type="Int32" Name="connectivity" format="ascii">\n');
    for i=1:size(ElementNumbers,1)
        fprintf(fid,'        %i %i %i %i %i %i \n',int32(ElementNumbers(i,:)));
    end
    
else
   fprintf(fid,'      <DataArray type="Int32" Name="connectivity" format="appended" offset="%i">\n',int32(Offset));
   Offset = Offset + length(ElementNumbers(:))*sizeof_UInt32 + 1*sizeof_UInt32;
    
end
fprintf(fid,'      </DataArray>\n');

% Offsets -----------
offsets = cumsum(ones(size(ElementNumbers,1),1)*6);
if ASCII
    % ASCII:
    fprintf(fid,'  <DataArray type="Int32" Name="offsets" format="ascii">\n');
    for i=1:size(ElementNumbers,1)
        fprintf(fid,'        %i\n',int32(offsets(i)));
    end

else
    fprintf(fid,'      <DataArray type="Int32" Name="offsets" format="appended" offset="%i">\n',int32(Offset));
    Offset = Offset + length(offsets(:))*sizeof_UInt32 + 1*sizeof_UInt32;
end
fprintf(fid,'      </DataArray>\n');

% types -----------
types = ones(size(ElementNumbers,1),1)*13;
if ASCII
    % ASCII:
    fprintf(fid,'      <DataArray type="UInt8" Name="types" format="ascii">\n');
    for i=1:size(ElementNumbers,1)
        fprintf(fid,'        %i\n',uint8(13));
    end

else
    fprintf(fid,'      <DataArray type="UInt8" Name="types" format="appended" offset="%i">\n',int32(Offset));
    Offset = Offset + length(types(:))*sizeof_UInt8 + 1*sizeof_UInt32;
end
fprintf(fid,'      </DataArray>\n');


fprintf(fid,'    </Cells>\n');
% -----------------------


fprintf(fid,'    </Piece>\n');
fprintf(fid,'  </UnstructuredGrid>\n');


%--------------------------------------------------------------------------

if ~ASCII
    % Append binary data in raw format: the order in which data arrays are
    % added should be the same as how they are defined above
    fprintf(fid,'  <AppendedData encoding="raw"> \n');
    fprintf(fid,'_');

    % Add Temperature in binary format
    fwrite(fid,uint32(length(Temperature(:))*sizeof_Float32),'uint32');
    fwrite(fid,single(Temperature(:)).'          ,   'float32');
    
    if WriteResidualT
        % Add Viscosity in binary format
        fwrite(fid,uint32(length(resT_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(resT_vec(:)).'     ,   'float32');
    end
    
    if WriteViscosity
        % Add Viscosity in binary format
        fwrite(fid,uint32(length(Viscosity_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Viscosity_vec(:)).'     ,   'float32');
    end

    if WriteComposition
        % Add Composition in binary format
        fwrite(fid,uint32(length(Composition_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Composition_vec(:)).'   ,   'float32');
    end
    
    if WriteContRoot
        % Add Composition1 in binary format
        fwrite(fid,uint32(length(Composition1_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Composition1_vec(:)).'   ,   'float32');
    end
    
    if WriteContCrust
        % Add Composition2 in binary format
        fwrite(fid,uint32(length(Composition2_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Composition2_vec(:)).'   ,   'float32');
    end
    
    if WriteContBelt
        % Add Composition in binary format
        fwrite(fid,uint32(length(Composition3_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Composition3_vec(:)).'   ,   'float32');
    end
        
    if WriteContSuture
        % Add Composition in binary format
        fwrite(fid,uint32(length(Composition4_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Composition4_vec(:)).'   ,   'float32');
    end
    
    if WriteContPlot
        % Add Continent in binary format
        fwrite(fid,uint32(length(ContPlot_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(ContPlot_vec(:)).'   ,   'float32');
    end
    
    if WriteContinent
        % Add Continent in binary format
        fwrite(fid,uint32(length(Continent_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Continent_vec(:)).'   ,   'float32');
    end 
    
    if WriteTopography
        % Add Composition in binary format
        fwrite(fid,uint32(length(Topography_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Topography_vec(:)).'   ,   'float32');
    end
    
    if WriteTopoSG
        % Add Composition in binary format
        fwrite(fid,uint32(length(TopoSG_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(TopoSG_vec(:)).'   ,   'float32');
    end
    
    if WriteGeoid
        % Add Composition in binary format
        fwrite(fid,uint32(length(Geoid_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Geoid_vec(:)).'   ,   'float32');
    end
    
   if WriteHeatflux
        % Add Composition in binary format
        fwrite(fid,uint32(length(HF_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(HF_vec(:)).'   ,   'float32');
    end
    
    if WriteStress
        % Add Stress in binary format
        fwrite(fid,uint32(length(Stress_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Stress_vec(:)).'  ,   'float32');
    end
    
    if WriteEdot
        % Add Strain rate in binary format
        fwrite(fid,uint32(length(Strainrate_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Strainrate_vec(:)).'  ,   'float32');
    end
    
    if WriteDamage
        % Add Strain rate in binary format
        fwrite(fid,uint32(length(Damage_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Damage_vec(:)).'  ,   'float32');
    end
    
    if WriteVelocity
        % Add Pressure in binary format
        fwrite(fid,uint32(length(Pressure_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Pressure_vec(:)).'      ,   'float32');

        % Add Velocity in binary format
        fwrite(fid,uint32(length(Velocity_vec(:))*sizeof_Float32),'uint32');
        fwrite(fid,single(Velocity_vec).' ,   'float32');
    end
    
    if WriteResidue
        % Add Pressure in binary format
        fwrite(fid,uint32(length(ResidueP_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(ResidueP_vec(:)).'      ,   'float32');

        % Add Velocity in binary format
        fwrite(fid,uint32(length(Residue_vec(:))*sizeof_Float32),'uint32');
        fwrite(fid,single(Residue_vec).' ,   'float32');
    end
    
    if WriteAge
        % Add Age in binary format
        fwrite(fid,uint32(length(Age_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Age_vec(:)).'     ,   'float32');
    end

    if WriteDiv
        % Add Age in binary format
        fwrite(fid,uint32(length(Div_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Div_vec(:)).'     ,   'float32');
    end
    
    if WriteVor
        % Add Age in binary format
        fwrite(fid,uint32(length(Vor_vec)*sizeof_Float32),'uint32');
        fwrite(fid,single(Vor_vec(:)).'     ,   'float32');
    end    
    
    % Add Coordinates in binary format
    fwrite(fid,uint32(length(Points(:))*sizeof_Float32),'uint32');
    fwrite(fid,single(Points).' ,   'float32');
    
    % Add Element connectivity in binary format
    fwrite(fid,uint32(length(ElementNumbers(:))*sizeof_UInt32),'uint32');
    fwrite(fid,uint32(ElementNumbers).' ,   'int32');

    % Add offsets array in binary format
    fwrite(fid,uint32(length(offsets(:))*sizeof_UInt32),'uint32');
    fwrite(fid,uint32(offsets).' ,   'int32');
    
    % Add types array in binary format
    fwrite(fid,uint32(length(types(:))*sizeof_UInt32),'uint32');
    fwrite(fid,uint8(types).' ,   'int8');
    
    fprintf(fid,'  </AppendedData> \n');
end
fprintf(fid,'</VTKFile>\n');

fclose(fid);

if ~ASCII
    disp(['Created Binary YinYang XML-VTK output file ',fname_vtk])
else
    disp(['Created ASCII  YinYang XML-VTK output file ',fname_vtk])
end

cd(start_dir);                                  % Change back to original directory

