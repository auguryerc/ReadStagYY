function Create_PVD_file(FileNames,fname,directory)
%
% Creates a Paraview data file that contains all datasets and times
%
% Boris Kaus, Feb. 26, 2008


% Change to correct directory
start_dir       = pwd;
cd(directory);

fname_pvd       = [fname,'.pvd'];
fid             = fopen(fname_pvd,'w','b');           % note the 'b': not doing BigEndian does not work with MATLAB!
fprintf(fid,'<?xml version="1.0"?> \n');
fprintf(fid,'<VTKFile type="Collection" version="0.1" byte_order="BigEndian" >\n');
fprintf(fid,'<Collection>\n');

% Write all filenames & times
for i=1:length(FileNames);
    Data = FileNames{i};
    fprintf(fid,'  <DataSet timestep="%f" part="001" file="%s"  name="Asmb:FRAME"/>\n',single(Data{1}),Data{3});
end

fprintf(fid,'</Collection>\n');
fprintf(fid,'</VTKFile>\n');
fclose(fid);

% Change back to original directory
cd(start_dir);