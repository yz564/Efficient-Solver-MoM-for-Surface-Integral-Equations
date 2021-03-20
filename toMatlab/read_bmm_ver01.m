%Function for reading binary matrix market format file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filename is the path of bmm file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%By Wenhao Xu, Xi'an Jiaotong University & Duke University
function A=read_bmm_ver01(filename)

%Initialization
fid=fopen(filename,'rb');
is_sparse=fread(fid,1,'bool');
is_real=fread(fid,1,'bool');
row_num=fread(fid,1,'int32');
col_num=fread(fid,1,'int32');

%Read data
if is_sparse
    nonzero_num=fread(fid,1,'int32');
    row_ind=fread(fid,nonzero_num,'int32')+1;
    col_ind=fread(fid,nonzero_num,'int32')+1;
    if is_real
        val=fread(fid,nonzero_num,'double');
    else
        real_val=fread(fid,nonzero_num,'double');
        imag_val=fread(fid,nonzero_num,'double');
        val=complex(real_val,imag_val);
    end
    A=sparse(row_ind,col_ind,val,row_num,col_num);
else
    if is_real
        A=fread(fid,[row_num,col_num],'double');
    else
        real_A=fread(fid,[row_num,col_num],'double');
        imag_A=fread(fid,[row_num,col_num],'double');
        A=complex(real_A,imag_A);
    end
end

%Close file
fclose(fid);

end