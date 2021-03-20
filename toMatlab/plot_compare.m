Cpp=read_bmm_ver01("ScatterFields_two_plate_v3.bmm");
Ex=Cpp(:,1);
Ey=Cpp(:,2);
Ez=Cpp(:,3);

Hx=Cpp(:,4);
Hy=Cpp(:,5);
Hz=Cpp(:,6);

Ex0=reshape(Ex,51,51).';
Ey0=reshape(Ey,51,51).';
Ez0=reshape(Ez,51,51).';
Hx0=reshape(Hx,51,51).';
Hy0=reshape(Hy,51,51).';
Hz0=reshape(Hz,51,51).';

figure(1)
subplot(2,3,1)
imagesc(abs(Ex0))
colormap jet
colorbar

subplot(2,3,2)
imagesc(abs(Ey0))
colormap jet
colorbar

subplot(2,3,3)
imagesc(abs(Ez0))
colormap jet
colorbar



subplot(2,3,4)
imagesc(abs(Hx0))
colormap jet
colorbar
subplot(2,3,5)
imagesc(abs(Hy0))
colormap jet
colorbar
subplot(2,3,6)
imagesc(abs(Hz0))
colormap jet
colorbar



path='./two_plate_v3_NearField1.hfe'
fid=fopen(path,'r');
if fid<0
 disp('打开文件失败!');
 return;
else
FormatString=repmat('%f ',1,9);
a =cell2mat(textscan(fid,FormatString,'HeaderLines',17)); %从174行开始读取7373*18的矩阵数据
end
fclose(fid);


Hx_s_feko=a(:,4)+1i*a(:,5); Hx_1_image0=reshape(Hx_s_feko,51,51);
Hy_s_feko=a(:,6)+1i*a(:,7); Hy_1_image0=reshape(Hy_s_feko,51,51);
Hz_s_feko=a(:,8)+1i*a(:,9); Hz_1_image0=reshape(Hz_s_feko,51,51);


figure(2)
subplot(2,3,1)
imagesc(abs(Hx_1_image0))
colormap jet
colorbar

subplot(2,3,2)
imagesc(abs(Hy_1_image0))
colormap jet
colorbar

subplot(2,3,3)
imagesc(abs(Hz_1_image0))
colormap jet
colorbar


subplot(2,3,4)
imagesc(abs(Hx0))
colormap jet
colorbar
subplot(2,3,5)
imagesc(abs(Hy0))
colormap jet
colorbar
subplot(2,3,6)
imagesc(abs(Hz0))
colormap jet
colorbar


errHx=norm(reshape((Hx_1_image0),51*51,1)-reshape((Hx0),51*51,1))/norm(reshape(Hx_1_image0,51*51,1))
errHy=norm(reshape((Hy_1_image0),51*51,1)-reshape((Hy0),51*51,1))/norm(reshape(Hy_1_image0,51*51,1))
errHz=norm(reshape((Hz_1_image0),51*51,1)-reshape((Hz0),51*51,1))/norm(reshape(Hz_1_image0,51*51,1))