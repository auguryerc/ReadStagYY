% Write STAG_3D format from matlab in XML-VTK format, using ASCII or BINARY file
% format
%
%
%
% Boris Kaus, 21 Feb. 2008


ASCII           = logical(0);                       % ASCII or BINARY?

% Change to correct directory
start_dir       = pwd;
cd(directory);


% Generate input data in correct format
Points          = single([X_3D(:) Y_3D(:) Z_3D(:)]);
T_3D            = single(T_3D);


% Transform the arrays in vector-format, and change them to the appropriate
% precision. This example file assumes that data comes in single precision
%
% Double precision requires changing in ths code below (Float32 into Float64)
T               = single(T_3D(:));                                                  % Temperature
Points          = [single(X_3D(:))  single(Y_3D(:))     single(Z_3D(:))     ];      % Coordinates

% Define which arrays to write
if WriteVelocity
    P               = single(P_3D(:));                                                  % Pressure
    Velocity        = [single(VX_3D(:)) single(VY_3D(:))    single(VZ_3D(:))    ];      % Velocity                                       % Velocity @ every point
end

if WriteViscosity
    Eta               = single(ETA_3D(:));                                                % Viscosity
end

if WriteComposition
     Composition      = single(C_3D(:));                                                 % Composition
end

if WriteContRoot
     Composition1     = single(C1_3D(:));                                                % Continental Root
end

if WriteContCrust
     Composition2     = single(C2_3D(:));                                                % Composition
end

if WriteContBelt
     Composition3     = single(C3_3D(:));                                                % Composition
end

if WriteContSuture
     Composition4     = single(C4_3D(:));                                                % Composition
end

if WriteContPlot
     ContPlot         = single(CPL_3D(:));                                               % Composition
end

if WriteContinent
     Continent        = single(NRC_3D(:));                                               % Composition
end

if WriteTopography
    Topography        = single(TP_3D(:));
end

if WriteTopoSG
    TopoSG            = single(TPSG_3D(:));
end

if WriteGeoid
    Geoid             = single(G_3D(:));
end

if WriteStress
    Stress            = single(S_3D(:));                                                  % Stress
end

if WriteEdot
    Strainrate           = single(E_3D(:));                                               % Strainrate
end

if WriteDamage
    Damage           = single(D_3D(:));                                               % Strainrate
end

if WriteAge
    Age           = single(AGE_3D(:));                                               % Strainrate
end

%==========================================================================
% Definitions and initialization
sizeof_Float32  =   4;      
sizeof_Float64  =   4;     
sizeof_UInt32   =   4;     
Offset          =   0;      % Initial offset

%--------------------------------------------------------------------------
% Write the header for a structured grid:
%--------------------------------------------------------------------------
fname_vtk       = [fname,'_',num2str(1000000+fname_number),'.vts'];
fid             = fopen(fname_vtk,'w','b');           % note the 'b': not doing BigEndian does not work with MATLAB!
fprintf(fid,'<?xml version="1.0"?> \n');
fprintf(fid,'<VTKFile type="StructuredGrid" version="0.1" byte_order="BigEndian" >\n');
fprintf(fid,'  <StructuredGrid  WholeExtent="%i %i %i %i %i %i">\n', [0 size(X_3D,1)-1 0 size(X_3D,2)-1 0 size(X_3D,3)-1]);
fprintf(fid,'  <Piece Extent="%i %i %i %i %i %i">\n', [0 size(X_3D,1)-1 0 size(X_3D,2)-1 0 size(X_3D,3)-1]);

%--------------------------------------------------------------------------
% Add point-wise data
%--------------------------------------------------------------------------
fprintf(fid,'    <PointData Scalars="T" Vectors="Velocity"  >\n');

% TEMPERATURE -----------
if ASCII
    % ASCII:
    fprintf(fid,'      <DataArray type="Float32" Name="T" format="ascii">\n');
    for i=1:length(T)
        fprintf(fid,'        %g \n',single(T(i)));
    end
else
    % BINARY:
    fprintf(fid,'      <DataArray type="Float32" Name="T" format="appended" offset="%i">\n', int32(Offset));
    Offset = Offset + length(T(:))*sizeof_Float32 + 1*sizeof_UInt32;
end
fprintf(fid,'      </DataArray>\n');
% -----------------------

if WriteViscosity
    % VISCOSITY---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Viscosity" format="ascii">\n');
        for i=1:length(Eta)
            fprintf(fid,'        %g \n',single(Eta(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Viscosity" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Eta(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteComposition
    % COMPOSITION---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Composition" format="ascii">\n');
        for i=1:length(Composition)
            fprintf(fid,'        %g \n',single(Composition(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Composition" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Composition(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end
    
if WriteContRoot
    % COMPOSITION1-----------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Root" format="ascii">\n');
        for i=1:length(Composition1)
           fprintf(fid,'        %g \n',single(Composition1(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Root" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Composition1(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteContCrust
    % COMPOSITION1-----------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Crust" format="ascii">\n');
        for i=1:length(Composition2)
           fprintf(fid,'        %g \n',single(Composition2(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Crust" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Composition2(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end


if WriteContBelt
    % COMPOSITION1-----------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Belt" format="ascii">\n');
        for i=1:length(Composition3)
           fprintf(fid,'        %g \n',single(Composition3(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Belt" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Composition3(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteContSuture
    % COMPOSITION1-----------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Suture" format="ascii">\n');
        for i=1:length(Composition4)
           fprintf(fid,'        %g \n',single(Composition4(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Cont Suture" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Composition4(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteContPlot
    % COMPOSITION1-----------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="cont plot" format="ascii">\n');
        for i=1:length(ContPlot)
           fprintf(fid,'        %g \n',single(ContPlot(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="cont plot" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(ContPlot(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteContinent
    % COMPOSITION1-----------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Continent" format="ascii">\n');
        for i=1:length(Continent)
           fprintf(fid,'        %g \n',single(Continent(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Continent" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Continent(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteTopography
   % Topography---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Topography" format="ascii">\n');
        for i=1:length(Topography)
            fprintf(fid,'        %g \n',single(Topography(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Topography" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Topography(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteTopoSG
   % Topography self-grav ---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Topo self-grav" format="ascii">\n');
        for i=1:length(TopoSG)
            fprintf(fid,'        %g \n',single(TopoSG(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Topo self-grav" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(TopoSG(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteGeoid
   % Geoid---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Geoid" format="ascii">\n');
        for i=1:length(Geoid)
            fprintf(fid,'        %g \n',single(Geoid(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Geoid" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Geoid(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end



if WriteStress
   % STRESS---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Stress" format="ascii">\n');
        for i=1:length(Stress)
            fprintf(fid,'        %g \n',single(Stress(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Stress" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Stress(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteEdot
   % STRAIN RATE---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Strain rate" format="ascii">\n');
        for i=1:length(Strainrate)
            fprintf(fid,'        %g \n',single(Strainrate(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Strainrate" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Strainrate(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteDamage
   % Damage---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Damage" format="ascii">\n');
        for i=1:length(Damage)
            fprintf(fid,'        %g \n',single(Damage(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Damage" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Damage(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end

if WriteAge
   % Damage---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Age" format="ascii">\n');
        for i=1:length(Age)
            fprintf(fid,'        %g \n',single(Age(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Age" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Age(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
end


if WriteVelocity
    
    % PRESSURE---------------
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="P" format="ascii">\n');
        for i=1:length(P)
            fprintf(fid,'        %g \n',single(P(i)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="P" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(P(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------

    % VELOCITY---------------  : VELOCITY IS A 3-component vector
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">\n');
        for i=1:length(T)
            fprintf(fid,'   %g %g %g \n',single(Velocity(i,:)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="appended" offset="%i">\n',int32(Offset));
        Offset = Offset + length(Velocity(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
    % -----------------------
    
end


fprintf(fid,'    </PointData>\n');
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

fprintf(fid,'    <Celldata>\n');
fprintf(fid,'    </Celldata>\n');


%--------------------------------------------------------------------------
% Add coordinates of structured grid 
%--------------------------------------------------------------------------
fprintf(fid,'    <Points>\n');

% ASCII
if ASCII
    fprintf(fid,'      <DataArray type="Float32" Name="Array" NumberOfComponents="3" format="ascii">\n');
    for i=1:size(Points,1)
        fprintf(fid,' %g %g %g \n',Points(i,:));
    end
else
    fprintf(fid,'      <DataArray type="Float32" Name="Array" NumberOfComponents="3" format="appended" offset="%i" >\n',int32(Offset));
    
end
fprintf(fid,'      </DataArray>\n');
fprintf(fid,'    </Points>\n');
%--------------------------------------------------------------------------

fprintf(fid,'  </Piece> \n');
fprintf(fid,'  </StructuredGrid> \n');


if ~ASCII
    % Append binary data in raw format: the order in which data arrays are
    % added should be the same as how they are defined above
    fprintf(fid,'  <AppendedData encoding="raw"> \n');
    fprintf(fid,'_');

    % Add Temperature in binary format
    fwrite(fid,uint32(length(T)*sizeof_Float32),'uint32');
    fwrite(fid,single(T).'      ,   'float32');
    
    if WriteViscosity
        % Add Viscosity in binary format
        fwrite(fid,uint32(length(Eta)*sizeof_Float32),'uint32');
        fwrite(fid,single(Eta).'      ,   'float32');
    end
    
    if WriteComposition
        % Add Composition in binary format
        fwrite(fid,uint32(length(Composition)*sizeof_Float32),'uint32');
        fwrite(fid,single(Composition).'      ,   'float32');
    end
    
    if WriteContRoot
        % Add Composition in binary format
        fwrite(fid,uint32(length(Composition1)*sizeof_Float32),'uint32');
        fwrite(fid,single(Composition1).'      ,   'float32');
    end
    
    if WriteContCrust
        % Add Composition in binary format
        fwrite(fid,uint32(length(Composition2)*sizeof_Float32),'uint32');
        fwrite(fid,single(Composition2).'      ,   'float32');
    end
    
    if WriteContBelt
        % Add Composition in binary format
        fwrite(fid,uint32(length(Composition3)*sizeof_Float32),'uint32');
        fwrite(fid,single(Composition3).'      ,   'float32');
    end
    
    if WriteContSuture
        % Add Composition in binary format
        fwrite(fid,uint32(length(Composition4)*sizeof_Float32),'uint32');
        fwrite(fid,single(Composition4).'      ,   'float32');
    end
    
    if WriteContPlot
        % Add Composition in binary format
        fwrite(fid,uint32(length(ContPlot)*sizeof_Float32),'uint32');
        fwrite(fid,single(ContPlot).'      ,   'float32');
    end
    
    if WriteContinent
        % Add Composition in binary format
        fwrite(fid,uint32(length(Continent)*sizeof_Float32),'uint32');
        fwrite(fid,single(Continent).'      ,   'float32');
    end
    
    if WriteTopography
        % Add Topography in binary format
        fwrite(fid,uint32(length(Topography)*sizeof_Float32),'uint32');
        fwrite(fid,single(Topography).'      ,   'float32');
    end
    
    if WriteTopoSG
        % Add TopoSG in binary format
        fwrite(fid,uint32(length(TopoSG)*sizeof_Float32),'uint32');
        fwrite(fid,single(TopoSG).'      ,   'float32');
    end
    
    if WriteGeoid
        % Add Geoid in binary format
        fwrite(fid,uint32(length(Geoid)*sizeof_Float32),'uint32');
        fwrite(fid,single(Geoid).'      ,   'float32');
    end
    
    if WriteStress
        % Add Stress in binary format
        fwrite(fid,uint32(length(Stress)*sizeof_Float32),'uint32');
        fwrite(fid,single(Stress).'      ,   'float32');
    end
    
    if WriteEdot
        % Add Strainrate in binary format
        fwrite(fid,uint32(length(Strainrate)*sizeof_Float32),'uint32');
        fwrite(fid,single(Strainrate).'      ,   'float32');
    end
    
    if WriteDamage
        % Add Strainrate in binary format
        fwrite(fid,uint32(length(Damage)*sizeof_Float32),'uint32');
        fwrite(fid,single(Damage).'      ,   'float32');
    end
    
    if WriteAge
        % Add Age in binary format
        fwrite(fid,uint32(length(Age)*sizeof_Float32),'uint32');
        fwrite(fid,single(Age).'      ,   'float32');
    end
    
    if WriteVelocity
        % Add Pressure in binary format
        fwrite(fid,uint32(length(P)*sizeof_Float32),'uint32');
        fwrite(fid,single(P).'      ,   'float32');

        % Add Velocity in binary format
        fwrite(fid,uint32(length(Velocity(:))*sizeof_Float32),'uint32');
        fwrite(fid,single(Velocity).' ,   'float32');
    end

    % Add Coordinates in binary format
    fwrite(fid,uint32(length(Points(:))*sizeof_Float32),'uint32');
    fwrite(fid,single(Points).' ,   'float32');

    fprintf(fid,'  </AppendedData> \n');
end


fprintf(fid,'</VTKFile>\n');
fclose(fid);

cd(start_dir)

if ~ASCII
%    disp(['Created Binary XML-VTK output file ',fname_vtk]);
else
%    disp(['Created ASCII XML-VTK output file ',fname_vtk]);
end