%
% usg_5x4x3.m
%   To generate an unstructured grid for PFLOTRAN
%
% Date: 18-July, 2011
% Author: Gautam Bisht (bishtg@ornl.gov)
%

function usg_5x4x3()

nx = 5;
ny = 4;
nz = 3;
v_x = [0 10 21 33 46 60];
v_y = [0 13 25 36 46];
v_z = [0 15 35 60];

% Compute Vertex Id
count = 0;
for kk = 1 : nz+1
    for jj = 1 : ny+1
        for ii = 1 : nx+1
            count = count + 1;
            vert_id(ii,jj,kk) = count; % 1-based
            vert_1d(count,:) = [v_x(ii) v_y(jj) v_z(kk)];
        end
    end
end

% Compute vertices forming a cell
count = 0;
for kk = 1 : nz
    for jj = 1 : ny
        for ii = 1 : nx
            count = count + 1;
            cell_id(ii,jj,kk) = count; % 1-based
            cell_connections(count,1) = 8;
            cell_connections(count,2) = vert_id(ii  , jj  , kk  );
            cell_connections(count,3) = vert_id(ii+1, jj  , kk  );
            cell_connections(count,4) = vert_id(ii+1, jj+1, kk  );
            cell_connections(count,5) = vert_id(ii  , jj+1, kk  );
            cell_connections(count,6) = vert_id(ii  , jj  , kk+1);
            cell_connections(count,7) = vert_id(ii+1, jj  , kk+1);
            cell_connections(count,8) = vert_id(ii+1, jj+1, kk+1);
            cell_connections(count,9) = vert_id(ii  , jj+1, kk+1);
        end
    end
end

% Top regions
count = 0;
kk = nz;
for jj = 1 : ny
    for ii = 1 : nx
        count = count + 1;
        top_cell_id(count,1) = cell_id(ii,jj,kk);
        top_sideset(count,1) = 4;
        top_sideset(count,2) = cell_connections( cell_id(ii,jj,kk) , 6);
        top_sideset(count,3) = cell_connections( cell_id(ii,jj,kk) , 7);
        top_sideset(count,4) = cell_connections( cell_id(ii,jj,kk) , 8);
        top_sideset(count,5) = cell_connections( cell_id(ii,jj,kk) , 9);
    end
end

% Bottom regions
count = 0;
kk = 1;
for jj = 1 : ny
    for ii = 1 : nx
        count = count + 1;
        bot_cell_id(count,1) = cell_id(ii,jj,kk);
        bot_sideset(count,1) = 4;
        bot_sideset(count,2) = cell_connections( cell_id(ii,jj,kk) , 2);
        bot_sideset(count,3) = cell_connections( cell_id(ii,jj,kk) , 5);
        bot_sideset(count,4) = cell_connections( cell_id(ii,jj,kk) , 4);
        bot_sideset(count,5) = cell_connections( cell_id(ii,jj,kk) , 3);
    end
end

% West regions
count = 0;
ii = 1;
for kk = 1 : nz
    for jj = 1 : ny
        count = count + 1;
        west_cell_id(count,1) = cell_id(ii,jj,kk);
        west_sideset(count,1) = 4;
        west_sideset(count,2) = cell_connections( cell_id(ii,jj,kk) , 5);
        west_sideset(count,3) = cell_connections( cell_id(ii,jj,kk) , 2);
        west_sideset(count,4) = cell_connections( cell_id(ii,jj,kk) , 6);
        west_sideset(count,5) = cell_connections( cell_id(ii,jj,kk) , 9);
    end
end

% East regions
count = 0;
ii = nx;
for kk = 1 : nz
    for jj = 1 : ny
        count = count + 1;
        east_cell_id(count,1) = cell_id(ii,jj,kk);
        east_sideset(count,1) = 4;
        east_sideset(count,2) = cell_connections( cell_id(ii,jj,kk) , 3);
        east_sideset(count,3) = cell_connections( cell_id(ii,jj,kk) , 4);
        east_sideset(count,4) = cell_connections( cell_id(ii,jj,kk) , 8);
        east_sideset(count,5) = cell_connections( cell_id(ii,jj,kk) , 7);
    end
end

% South regions
count = 0;
jj = 1;
for kk = 1 : nz
    for ii = 1 : nx
        count = count + 1;
        south_cell_id(count,1) = cell_id(ii,jj,kk);
        south_sideset(count,1) = 4;
        south_sideset(count,2) = cell_connections( cell_id(ii,jj,kk) , 2);
        south_sideset(count,3) = cell_connections( cell_id(ii,jj,kk) , 3);
        south_sideset(count,4) = cell_connections( cell_id(ii,jj,kk) , 7);
        south_sideset(count,5) = cell_connections( cell_id(ii,jj,kk) , 6);
    end
end

% North regions
count = 0;
jj = ny;
for kk = 1 : nz
    for ii = 1 : nx
        count = count + 1;
        north_cell_id(count,1) = cell_id(ii,jj,kk);
        north_sideset(count,1) = 4;
        north_sideset(count,2) = cell_connections( cell_id(ii,jj,kk) , 4);
        north_sideset(count,3) = cell_connections( cell_id(ii,jj,kk) , 5);
        north_sideset(count,4) = cell_connections( cell_id(ii,jj,kk) , 9);
        north_sideset(count,5) = cell_connections( cell_id(ii,jj,kk) , 8);
    end
end

tmp(:,1) = uint64(cell_connections(:,1));
tmp(:,2) = uint64(cell_connections(:,2));
tmp(:,3) = uint64(cell_connections(:,3));
tmp(:,4) = uint64(cell_connections(:,4));
tmp(:,5) = uint64(cell_connections(:,5));
tmp(:,6) = uint64(cell_connections(:,6));
tmp(:,7) = uint64(cell_connections(:,7));
tmp(:,8) = uint64(cell_connections(:,8));
tmp(:,9) = uint64(cell_connections(:,9));

bot_sideset_1(:,1) = uint64(bot_sideset(:,1));
bot_sideset_1(:,2) = uint64(bot_sideset(:,2));
bot_sideset_1(:,3) = uint64(bot_sideset(:,3));
bot_sideset_1(:,4) = uint64(bot_sideset(:,4));
bot_sideset_1(:,5) = uint64(bot_sideset(:,5));

top_sideset_1(:,1) = uint64(top_sideset(:,1));
top_sideset_1(:,2) = uint64(top_sideset(:,2));
top_sideset_1(:,3) = uint64(top_sideset(:,3));
top_sideset_1(:,4) = uint64(top_sideset(:,4));
top_sideset_1(:,5) = uint64(top_sideset(:,5));

east_sideset_1(:,1) = uint64(east_sideset(:,1));
east_sideset_1(:,2) = uint64(east_sideset(:,2));
east_sideset_1(:,3) = uint64(east_sideset(:,3));
east_sideset_1(:,4) = uint64(east_sideset(:,4));
east_sideset_1(:,5) = uint64(east_sideset(:,5));

north_sideset_1(:,1)= uint64(north_sideset(:,1));
north_sideset_1(:,2) = uint64(north_sideset(:,2));
north_sideset_1(:,3) = uint64(north_sideset(:,3));
north_sideset_1(:,4) = uint64(north_sideset(:,4));
north_sideset_1(:,5) = uint64(north_sideset(:,5));

south_sideset_1(:,1) = uint64(south_sideset(:,1));
south_sideset_1(:,2) = uint64(south_sideset(:,2));
south_sideset_1(:,3) = uint64(south_sideset(:,3));
south_sideset_1(:,4) = uint64(south_sideset(:,4));
south_sideset_1(:,5) = uint64(south_sideset(:,5));

west_sideset_1(:,1) = uint64(west_sideset(:,1));
west_sideset_1(:,2) = uint64(west_sideset(:,2));
west_sideset_1(:,3) = uint64(west_sideset(:,3));
west_sideset_1(:,4) = uint64(west_sideset(:,4));
west_sideset_1(:,5) = uint64(west_sideset(:,5));

south_faceid = [south_cell_id ones(size(south_cell_id))*1];
east_faceid  = [east_cell_id  ones(size(east_cell_id) )*2];
north_faceid = [north_cell_id ones(size(north_cell_id))*3];
west_faceid  = [west_cell_id  ones(size(west_cell_id) )*4];
bot_faceid   = [bot_cell_id   ones(size(bot_cell_id)  )*5];
top_faceid   = [top_cell_id   ones(size(top_cell_id)  )*6];

all = [1:60]';

%
% Write ASCII files
%

% 1) Write the mesh
fid = fopen('usg_5x4x3.dat','w');
fprintf(fid,'%d %d\n',size(cell_connections,1),size(vert_1d,1))
for ii = 1 : size(cell_connections, 1)
    for jj = 1 : size(cell_connections, 2)
        fprintf(fid,'%d ',cell_connections(ii,jj));
    end
    fprintf(fid,'\n');
end

for ii = 1 : size(vert_1d, 1)
    for jj = 1 : size(vert_1d, 2)
        fprintf(fid,'%f ',vert_1d(ii,jj));
    end
    fprintf(fid,'\n');
end

fclose(fid);

% 2) Write regions:
%
% % % 2.1) Bottom
% % fid = fopen('Bottom.dat','w');
% % for ii = 1 : length(bot_cell_id)
% %     fprintf(fid,'%d\n',bot_cell_id(ii));
% % end
% % fclose(fid)
% % 
% % fid = fopen('Bottom_faceid.dat','w');
% % for ii = 1 : size(bot_faceid,1)
% %     fprintf(fid,'%d %d\n',bot_faceid(ii,1), bot_faceid(ii,2));
% % end
% % fclose(fid)
% % 
% % fid = fopen('Bottom_sidesets.dat','w');
% % for ii = 1 : size(bot_sideset,1)
% %     for jj = 1 : size(bot_sideset,2)
% %         fprintf(fid,'%d ',bot_sideset(ii,jj));
% %     end
% %     fprintf(fid,'\n');
% % end
% % fclose(fid)
% % 
% % % 2.2) Top
% % fid = fopen('Top.dat','w');
% % for ii = 1 : length(top_cell_id)
% %     fprintf(fid,'%d\n',top_cell_id(ii));
% % end
% % fclose(fid)
% % 
% % fid = fopen('Top_faceid.dat','w');
% % for ii = 1 : size(top_faceid,1)
% %     fprintf(fid,'%d %d\n',top_faceid(ii,1), top_faceid(ii,2));
% % end
% % fclose(fid)
% % 
% % fid = fopen('Top_sidesets.dat','w');
% % for ii = 1 : size(top_sideset,1)
% %     for jj = 1 : size(top_sideset,2)
% %         fprintf(fid,'%d ',top_sideset(ii,jj));
% %     end
% %     fprintf(fid,'\n');
% % end
% % fclose(fid)
% % 
% % % 2.3) South
% % fid = fopen('South.dat','w');
% % for ii = 1 : length(south_cell_id)
% %     fprintf(fid,'%d\n',south_cell_id(ii));
% % end
% % fclose(fid)
% % 
% % fid = fopen('South_faceid.dat','w');
% % for ii = 1 : size(south_faceid,1)
% %     fprintf(fid,'%d %d\n',south_faceid(ii,1), south_faceid(ii,2));
% % end
% % fclose(fid)
% % 
% % fid = fopen('South_sidesets.dat','w');
% % for ii = 1 : size(south_sideset,1)
% %     for jj = 1 : size(south_sideset,2)
% %         fprintf(fid,'%d ',south_sideset(ii,jj));
% %     end
% %     fprintf(fid,'\n');
% % end
% % fclose(fid)
% % 
% % % 2.4) North
% % fid = fopen('North.dat','w');
% % for ii = 1 : length(north_cell_id)
% %     fprintf(fid,'%d\n',north_cell_id(ii));
% % end
% % fclose(fid)
% % 
% % fid = fopen('North_faceid.dat','w');
% % for ii = 1 : size(north_faceid,1)
% %     fprintf(fid,'%d %d\n',north_faceid(ii,1), north_faceid(ii,2));
% % end
% % fclose(fid)
% % 
% % fid = fopen('North_sidesets.dat','w');
% % for ii = 1 : size(north_sideset,1)
% %     for jj = 1 : size(north_sideset,2)
% %         fprintf(fid,'%d ',north_sideset(ii,jj));
% %     end
% %     fprintf(fid,'\n');
% % end
% % fclose(fid)
% % 
% % % 2.5) West
% % fid = fopen('West.dat','w');
% % for ii = 1 : length(west_cell_id)
% %     fprintf(fid,'%d\n',west_cell_id(ii));
% % end
% % fclose(fid)
% % 
% % fid = fopen('West_faceid.dat','w');
% % for ii = 1 : size(west_faceid,1)
% %     fprintf(fid,'%d %d\n',west_faceid(ii,1), west_faceid(ii,2));
% % end
% % fclose(fid)
% % 
% % fid = fopen('West_sidesets.dat','w');
% % for ii = 1 : size(west_sideset,1)
% %     for jj = 1 : size(west_sideset,2)
% %         fprintf(fid,'%d ',west_sideset(ii,jj));
% %     end
% %     fprintf(fid,'\n');
% % end
% % fclose(fid)
% % 
% % % 2.6) East
% % fid = fopen('East.dat','w');
% % for ii = 1 : length(east_cell_id)
% %     fprintf(fid,'%d\n',east_cell_id(ii));
% % end
% % fclose(fid)
% % 
% % fid = fopen('East_faceid.dat','w');
% % for ii = 1 : size(east_faceid,1)
% %     fprintf(fid,'%d %d\n',east_faceid(ii,1), east_faceid(ii,2));
% % end
% % fclose(fid)
% % 
% % fid = fopen('East_sidesets.dat','w');
% % for ii = 1 : size(east_sideset,1)
% %     for jj = 1 : size(east_sideset,2)
% %         fprintf(fid,'%d ',east_sideset(ii,jj));
% %     end
% %     fprintf(fid,'\n');
% % end
% % fclose(fid)
% % 
% % % 2.7) All
% % fid = fopen('All.dat','w');
% % for ii = 1 : length(all)
% %     fprintf(fid,'%d\n',all(ii))
% % end
% % fclose(fid);

%
% Write HDF5 file
%
hdf5write('usg_5x4x3.h5', '/Domain/Cells', tmp', ...
    '/Domain/Vertices',vert_1d', ...
    '/Regions/Bottom/Cell Ids',uint64(bot_cell_id), ...
    '/Regions/Top/Cell Ids',uint64(top_cell_id), ...
    '/Regions/South/Cell Ids',uint64(south_cell_id), ...
    '/Regions/North/Cell Ids',uint64(north_cell_id), ...
    '/Regions/West/Cell Ids',uint64(west_cell_id), ...
    '/Regions/East/Cell Ids',uint64(east_cell_id), ...
    '/Regions/Bottom_sidesets/Vertex Ids',uint64(bot_sideset_1'), ...
    '/Regions/Top_sidesets/Vertex Ids',uint64(top_sideset_1'), ...
    '/Regions/South_sidesets/Vertex Ids',uint64(south_sideset_1'), ...
    '/Regions/North_sidesets/Vertex Ids',uint64(north_sideset_1'), ...
    '/Regions/West_sidesets/Vertex Ids',uint64(west_sideset_1'), ...
    '/Regions/East_sidesets/Vertex Ids',uint64(east_sideset_1'), ...
    '/Regions/Bottom_faceid/Cell Ids',uint64(bot_cell_id'), ...
    '/Regions/Bottom_faceid/Face Ids',uint64(bot_faceid'), ...
    '/Regions/Top_faceid/Cell Ids',uint64(top_cell_id'), ...
    '/Regions/Top_faceid/Face Ids',uint64(top_faceid(:,2)'), ...
    '/Regions/South_faceid/Cell Ids',uint64(south_cell_id'), ...
    '/Regions/South_faceid/Face Ids',uint64(south_faceid(:,2)'), ...
    '/Regions/North_faceid/Cell Ids',uint64(north_cell_id'), ...
    '/Regions/North_faceid/Face Ids',uint64(north_faceid(:,2)'), ...
    '/Regions/West_faceid/Cell Ids',uint64(west_cell_id'), ...
    '/Regions/West_faceid/Face Ids',uint64(west_faceid(:,2)'), ...
    '/Regions/East_faceid/Cell Ids',uint64(east_cell_id'), ...
    '/Regions/East_faceid/Face Ids',uint64(east_faceid(:,2)'), ...
    '/Regions/All/Cell Ids',uint64(all)');

