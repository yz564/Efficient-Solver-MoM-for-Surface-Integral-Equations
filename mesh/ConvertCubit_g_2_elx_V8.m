function ConvertCubit_g_2_elx_V8(nstr, RotYDeg)
%%% datestr(now) = 13-Jul-2018 16:38:55
%%% search for *.g files & convert them into *.elx files automatically
% ConvertCubit_g_2_elx_V8('*RotY30*', 30)
% ConvertCubit_g_2_elx_V8('*RotY60*', 60)
clc; %clear;  
if(nargin < 1 )
    nstr = '*';
    RotYDeg = 0;
end
if(nargin < 2 ) 
    RotYDeg = 0;
end
global MatRotY
beta = RotYDeg*pi/180;
% MatRotY = roty( RotYDeg );
cosb = cos(beta);
sinb = sin(beta);
MatRotY = [ cosb 0 sinb; 0 1 0; -sinb 0 cosb];
%%
gtype = '.g'; etype = '.elx';
gv = dir([nstr, gtype]);
elxv = dir([nstr, etype]);
%% remove affixes
M = size(gv, 1);
for m = 1:M
    gname = gv(m).name; L = length(gname) - length(gtype);
    gnam = gname(1:L); gv(m).name = gnam;
end
N = size( elxv, 1 );
%% convert *.g to *.ncg
disp('Convert *.g to *.elx: '); 
for m = 1:M
    gname = gv(m).name; gexist = false;
    for n = 1:N
        gexist = strcmp([gname etype], elxv(n).name);
        if(gexist)
            break;
        end
    end
    if(~gexist)
        ncg2elxV3([gname gtype], [gname etype], false);
        fprintf(1, '\t%s\t%s\r', [gname, etype, ': '], 'converted successfully!');
    else
        fprintf(1, '\t%s\t%s\r', [gname, etype, ': '], 'NO need to be converted!');
    end
end

end

function [ConnectAttrib, NodXYZ, nelem_nnod_natt, elem_type_vec] = ncg2elxV3(ncg_file, file_out, pf)
%% read *.ncg (mesh file generated from Cubit output), and write *.elx files
% % for each subdomain: *.elx --- node indices of each elem & node coordinates
% % attrib* --- material ID
% % datestr(now) = 21-Jun-2016 20:15:01, by Breeze
% % Works for Tetra10 & Shell4, resulting from ncg2tetra10v2.m & ncg2quad4v2.m 
% clear; clc;
% file_out = 'Case1Quad4';
% ncg_file = 'Case1Quad4_v1.ncg';
% file_out = 'Case1Tetra10Subd1';
% ncg_file = 'Case1Tetra4_v7_Subd1.ncg';
% 24-Nov-2017 14:45:26: V3, Tri3 added
if(nargin< 3)
    pf = true;
end
global MatRotY
% clc; clear; fclose('all');
% ncg_file = 'SquarehLW5in_Tri1040.g';

finfo = ncinfo(ncg_file);
ncid = netcdf.open(ncg_file, 'NOWRITE');

% disp(finfo);
% FID = fopen(ncg_file);
%% check file content: num_nod_per_el* = 10
dim_num = size(finfo.Dimensions, 2);
for d = 1:dim_num
    dimid = d - 1;
    [tch, dimlen] = netcdf.inqDim(ncid,dimid);
    if( 1 == strcmp(tch, 'num_nodes')  )
    	num_node = dimlen;
    elseif( 1 == strcmp(tch, 'num_elem')  )
    	num_elem = dimlen;
    elseif( 1 == strcmp(tch, 'num_el_blk')  )
    	num_blk = dimlen;
    end
end

% % (:, [NO. of elem per block, NO. of node per elem, NO. of attributes])
nelem_nnod_natt = zeros(num_blk, 3);
for cnt = 1:num_blk
    blk_el_str = ['num_el_in_blk', int2str(cnt)];
    blk_nn_str = ['num_nod_per_el', int2str(cnt)];
    blk_at_str = ['num_att_in_blk', int2str(cnt)];
    
    for d = 1:dim_num
        dimid = d - 1;
        [tch, dimlen] = netcdf.inqDim(ncid,dimid);
        if( 1 == strcmp(tch, blk_el_str)  )
            nelem = dimlen; 
        elseif( 1 == strcmp(tch, blk_nn_str)  )
            nnod = dimlen; 
        elseif( 1 == strcmp(tch, blk_at_str)  )
            natt = dimlen; 
        end
    end
    
    nelem_nnod_natt(cnt, :) = [nelem, nnod, natt];
end

if( (sum(nelem_nnod_natt(:,1)) ~= num_elem)||(any(nelem_nnod_natt(:,3) ~= 1) ) )
    disp('ERROR: in subdomain element numbers or attributes !!!');
end
%% read attrib*, connect* & coord
var_num = size(finfo.Variables, 2);
att_connect_eltype = cell(num_blk, 3); %% attribute, connect ID, elem type
elem_type_vec = cell(num_blk, 1);
for cnt = 1:num_blk
    att_connect_eltype{cnt, 1} = ['attrib', int2str(cnt)];
    att_connect_eltype{cnt, 2} = ['connect', int2str(cnt)];
    for v = 1:var_num
        vid = v - 1;
        tch = netcdf.inqVar(ncid,vid);
        if( 1 == strcmp(tch, att_connect_eltype{cnt, 2})  )
            typestr = finfo.Variables(v).Attributes.Name;
            if( 1 == strcmp( typestr, 'elem_type') )
                typstr = finfo.Variables(v).Attributes.Value;
            else
                disp('ERROR: in Variables !!!');
            end
            break;
        end
    end 
    
    % % % 13-Jul-2018 16:09:09
    if(strcmp( typstr, 'SHELL4'))
        att_connect_eltype{cnt, 3} = 'QUAD4';
    elseif(strcmp( typstr, 'TETRA10'))
        att_connect_eltype{cnt, 3} = typstr;
    elseif(strcmp( typstr, 'SHELL9'))
        att_connect_eltype{cnt, 3} = 'QUAD9';
    elseif(strcmp( typstr, 'TRI3'))
        att_connect_eltype{cnt, 3} = typstr;
    else
        disp('ERROR: Element Type needs to be modified !!!');
        return;
    end
    elem_type_vec{cnt} = att_connect_eltype{cnt, 3};
end


ConnectAttrib = cell(num_blk, 2);
for cnt = 1:num_blk
    % % read attrib*
    cstr = att_connect_eltype{cnt, 1};
    for v = 1:var_num
        vid = v - 1;
        tch = netcdf.inqVar(ncid,vid);
        if( 1 == strcmp(tch, cstr)  )
            att_vec = netcdf.getVar(ncid, vid);
            break;
        end
    end
    
    nelem = nelem_nnod_natt(cnt, 1);
    natt = nelem_nnod_natt(cnt, 3);
    
    ConnectAttrib{cnt, 2} = transpose( reshape(att_vec, [natt nelem]) );
    
    % % read connect*
    cstr = att_connect_eltype{cnt, 2}; 
    for v = 1:var_num
        vid = v - 1;
        tch = netcdf.inqVar(ncid,vid);
        if( 1 == strcmp(tch, cstr)  )
            con_vec = netcdf.getVar(ncid, vid);
            break;
        end
    end
    nnod = nelem_nnod_natt(cnt, 2); %%% NO. of node per elem
    
    ConnectAttrib{cnt, 1} = transpose( reshape(con_vec, [nnod nelem]) );
end

% % read coordinates of node
NodXYZ = zeros(num_node, 3);
for v = 1:var_num
    vid = v - 1;
    tch = netcdf.inqVar(ncid,vid);
    if( 1 == strcmp(tch, 'coord')  )
        xyz_vec = netcdf.getVar(ncid, vid);
        NodXYZ = reshape(xyz_vec, [num_node 3]);
        break;
    elseif( 1 == strcmp(tch, 'coordx')  )
        xyz_vec = netcdf.getVar(ncid, vid);
        NodXYZ(:, 1) = xyz_vec ; 
    elseif( 1 == strcmp(tch, 'coordy')  )
        xyz_vec = netcdf.getVar(ncid, vid);
        NodXYZ(:, 2) = xyz_vec ; 
    elseif( 1 == strcmp(tch, 'coordz')  )
        xyz_vec = netcdf.getVar(ncid, vid);
        NodXYZ(:, 3) = xyz_vec ;
    end
end
%%
avgxyz = mean( mean( abs(NodXYZ) ) );

nodf0 = abs( NodXYZ ) < avgxyz*1.0E-12;
NodXYZ(nodf0) = 0;
%% format of Prefix+Sid.elx :
%   NO. of Elem, NO. of Nodes per Elem, Elem Type, NO. of Block
%   #Block ID,  NO. of Elem, attrib*
%   Node Indices per Elem
%   #Num_Nodes:  NO. of Nodes
%   NodeX  NodeY   NodeZ
L = length(file_out); affix = '.elx';
if(strcmp(file_out((L-3):L), affix))
    file_out = file_out(1:(L-4));
end
fwid = fopen([file_out, affix], 'w');
fprintf(fwid, '#Date: %s, by Breeze\r\n',  datestr(now));
fprintf(fwid, '%4d%8d%8d \t#num_blk, num_elem, num_node\r\n',...
        num_blk, num_elem, num_node);
fprintf(fwid, '%4d%12s \t#num_nod_per_elem, elem_type\r\n', ...
    nelem_nnod_natt(1,2), att_connect_eltype{1,3} );

wfmt = '';
for m = 1:nelem_nnod_natt(1,2)
    wfmt = [wfmt, '%0d ']; %#ok<AGROW>
end

for cnt = 1:num_blk
    elem_typ_str = att_connect_eltype{1,3};  L = length(elem_typ_str);
    while ((elem_typ_str(L)>='0') && (elem_typ_str(L)<='9'))
        elem_typ_str = elem_typ_str(1:(L-1)); L = length(elem_typ_str);
    end
    if( strcmp(elem_typ_str, 'SHELL'))        
        if( length( ConnectAttrib{cnt, 2}(1, :)) == 2 )
            % %       material_inside*1000 + material_outside
            shell_attrib = ConnectAttrib{cnt, 2}(1, 1)*1000 + ConnectAttrib{cnt, 2}(1, 2);
        else
            shell_attrib = ConnectAttrib{cnt, 2}(1);
        end
        fprintf(fwid, '%8d%4d%8d \t#num_elem_blk, blk_id, atrib_blk\r\n', ...
          nelem_nnod_natt(cnt,1), cnt, shell_attrib );
    else
        fprintf(fwid, '%8d%4d%8d \t#num_elem_blk, blk_id, atrib_blk\r\n', ...
         nelem_nnod_natt(cnt,1), cnt, ConnectAttrib{cnt, 2}(1) );
    end
    for m = 1:nelem_nnod_natt(cnt,1)
        fprintf(fwid, [wfmt, '\r\n'], ConnectAttrib{cnt, 1}(m,:) );
    end
end
%%% print node coordinates
% fprintf(fwid, '#num_node = %d\r\n', num_node);

p=NodXYZ';
t=con_vec;
save two_plate_v1  p ...
            t 
            
for m = 1:num_node
    v = NodXYZ(m,:);
    w = MatRotY*(v');
    fprintf(fwid, '%.16G %.16G %.16G \n', w);
end
%%
netcdf.close(ncid);
fclose(fwid); fclose('all');
if(pf)
    disp(['Mesh converted successfully: ', file_out, affix , ' !!!']);
end
end