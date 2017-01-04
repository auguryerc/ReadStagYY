function [fname_number,fname_vtk]=create_vtu(path,name,X_yin,Y_yin,Z_yin,X_yang,Y_yang,Z_yang,Fields,time)
start_dir       = pwd;
cd(path);%.inv_par

%==========================================================================
% 1) Take the surface of the 2 grids, patch together and triangulate
%==========================================================================
X_yin_s                 =   X_yin(:,:,end);
Y_yin_s                 =   Y_yin(:,:,end);
Z_yin_s                 =   Z_yin(:,:,end);

X_yang_s                 =   X_yang(:,:,end);
Y_yang_s                 =   Y_yang(:,:,end);
Z_yang_s                 =   Z_yang(:,:,end);

R_yin_s                  =   sqrt(X_yin_s.^2 + Y_yin_s.^2 + Z_yin_s.^2);
Th_yin_s              =   atan2(sqrt(X_yin_s.^2 + Y_yin_s.^2),Z_yin_s);
Ph_yin_s                =   atan2(Y_yin_s,X_yin_s);


R_yang_s                  =   sqrt(X_yang_s.^2 + Y_yang_s.^2 + Z_yang_s.^2);
Th_yang_s              =   atan2(sqrt(X_yang_s.^2 + Y_yang_s.^2),Z_yang_s);
Ph_yang_s                =   atan2(Y_yang_s,X_yang_s);


% Cut off the corners from grid #1, which seems to do #2 as well (PJT):
theta12             = acos(sin(Th_yin_s).*sin(Ph_yin_s));  % theta coords of grid 1 in coord system of grid 2
ind_corner          = find( (theta12>pi/4 & Ph_yin_s>pi/2) | (theta12<3*pi/4 & Ph_yin_s<-pi/2 ) );


% Form indices of remaining
ind                           = 1:numel(Ph_yin_s);
ind([ind_corner])=[];
%ind([ind_corner1 ind_corner2])=[];

% Create an array with unique r, theta & phi values
numberYin = ones(size(R_yin_s));
numberYin(find(R_yin_s==R_yin_s))           =   find(R_yin_s==R_yin_s);
numberYang = ones(size(R_yang_s));
numberYang(find(R_yang_s==R_yang_s))          =   find(R_yang_s==R_yang_s)+ max(numberYin(:));

%Closed Yin grid, w/o corners and completely closed
R_s                   =   [R_yin_s(ind),      ]';
Th_s               =   [Th_yin_s(ind)   ]';
Ph_s                 =   [Ph_yin_s(ind),    ]';
NumYin              =   [numberYin(ind)]';


% x,y,z coordinates of complete grid:
X_s                   =   [X_yin_s(ind) -X_yin_s(ind)];
Y_s                   =   [Y_yin_s(ind)  Z_yin_s(ind)];
Z_s                   =   [Z_yin_s(ind)  Y_yin_s(ind)];
NumYang             =   [NumYin + max(numberYin(:)) ];
Number              =   [NumYin(:);    NumYang(:)  ]';

%   remove redundant coordinates
tri                 =   convhulln([X_s(:), Y_s(:),Z_s(:)]);       % simple way to grid it
triYingYang         =   Number(tri);

X_s                   =   [X_yin_s(:); X_yang_s(:)];
Y_s                   =   [Y_yin_s(:); Y_yang_s(:)];
Z_s                   =   [Z_yin_s(:); Z_yang_s(:)];

%==========================================================================
% triYingYang now contains the numbers of all triangles
%==========================================================================

%==========================================================================
% 2) Create a 3D grid with tetrahedron elements
%==========================================================================

% Number all gridpoints we have
NUMBER_1                 = ones(size(X_yin));
NUMBER_2                 = ones(size(X_yang));
NUMBER_1(find(NUMBER_1)) = find(NUMBER_1);
NUMBER_2(find(NUMBER_2)) = find(NUMBER_2) + max(NUMBER_1(:));


% Make a loop over all levels
ElementNumbers      = [];
for iz=1:size(X_yang,3)-1
    
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
Points(NUMBER_1(:),1)       = X_yin(:);
Points(NUMBER_2(:),1)       = X_yang(:);
Points(NUMBER_1(:),2)       = Y_yin(:);
Points(NUMBER_2(:),2)       = Y_yang(:);
Points(NUMBER_1(:),3)       = Z_yin(:);
Points(NUMBER_2(:),3)       = Z_yang(:);

Radius                      = sqrt(sum(Points.^2,2));

Field_Names=fieldnames(Fields);
for k=1:length(Field_Names)
    Field_yin=Fields.(Field_Names{k})(:,1:length(X_yin(1,:,1)),:);
    Field_yang=Fields.(Field_Names{k})(:,length(X_yin(1,:,1))+1:2*length(X_yin(1,:,1)),:);
    eval(strcat(Field_Names{k},'=zeros(max(NUMBER_2(:)),1);'));
    eval(strcat(Field_Names{k},'(NUMBER_1(:))=Field_yin(:);'));
    eval(strcat(Field_Names{k},'(NUMBER_2(:))=Field_yang(:);'));
end


%==========================================================================
% 3) Write VTK file (unstructured mesh)
%==========================================================================
ElementNumbers  = ElementNumbers-1;      % VTK is zero-based


%==========================================================================
% Definitions and innitialization
sizeof_Float32  =   4;
sizeof_Float64  =   4;
sizeof_UInt32   =   4;
sizeof_UInt8    =   1;

Offset          =   0;      % Initial offset

%--------------------------------------------------------------------------
% Write the header for a structured grid:
%--------------------------------------------------------------------------
fname_number    =num2str(1000000+time)
fname_vtk       = [name,'_',fname_number,'.vtu'];
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
for k=1:length(Field_Names)
    fprintf(fid,strcat('      <DataArray type="Float32" Name="',Field_Names{k},'" format="appended" offset="%i">\n'), int32(Offset));
    eval(strcat('Offset = Offset + length(',Field_Names{k},'(:))*sizeof_Float32 + 1*sizeof_UInt32;'));
    fprintf(fid,'      </DataArray>\n');
end

% -----------------------

fprintf(fid,'    </PointData>\n');

%--------------------------------------------------------------------------
% Add coordinates of structured grid
%--------------------------------------------------------------------------
fprintf(fid,'    <Points>\n');

fprintf(fid,'      <DataArray type="Float32" Name="Array" NumberOfComponents="3" format="appended" offset="%i" >\n',int32(Offset));
Offset = Offset + length(Points(:))*sizeof_Float32 + 1*sizeof_UInt32;
fprintf(fid,'      </DataArray>\n');
fprintf(fid,'    </Points>\n');


%--------------------------------------------------------------------------
% Add CELLS data
%--------------------------------------------------------------------------
fprintf(fid,'    <Cells>\n');

% Connectivity -----------

fprintf(fid,'      <DataArray type="Int32" Name="connectivity" format="appended" offset="%i">\n',int32(Offset));
Offset = Offset + length(ElementNumbers(:))*sizeof_UInt32 + 1*sizeof_UInt32;

fprintf(fid,'      </DataArray>\n');

% Offsets -----------
offsets = cumsum(ones(size(ElementNumbers,1),1)*6);
fprintf(fid,'      <DataArray type="Int32" Name="offsets" format="appended" offset="%i">\n',int32(Offset));
Offset = Offset + length(offsets(:))*sizeof_UInt32 + 1*sizeof_UInt32;
fprintf(fid,'      </DataArray>\n');

% types -----------
types = ones(size(ElementNumbers,1),1)*13;
fprintf(fid,'      <DataArray type="UInt8" Name="types" format="appended" offset="%i">\n',int32(Offset));
Offset = Offset + length(types(:))*sizeof_UInt8 + 1*sizeof_UInt32;
fprintf(fid,'      </DataArray>\n');


fprintf(fid,'    </Cells>\n');
% -----------------------


fprintf(fid,'    </Piece>\n');
fprintf(fid,'  </UnstructuredGrid>\n');


%--------------------------------------------------------------------------

% Append binary data in raw format: the order in which data arrays are
% added should be the same as how they are defined above
fprintf(fid,'  <AppendedData encoding="raw"> \n');
fprintf(fid,'_');

% Add Temperature in binary format
for k=1:length(Field_Names)
    eval(strcat('Field=',Field_Names{k},';'));
    fwrite(fid,uint32(length(Field(:))*sizeof_Float32),'uint32');
    fwrite(fid,single(Field(:)).'          ,   'float32');
    
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

fprintf(fid,'</VTKFile>\n');

fclose(fid);

disp(['Created Binary YinYang XML-VTK output file ',fname_vtk])
cd(start_dir);                                  % Change back to original directory

