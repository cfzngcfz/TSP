%% ����������·������
% ���룺
% D     ��������֮��ľ���
% Chrom ����Ĺ켣
function len=ObjF(Chrom,D)
len=0;
temp=[Chrom Chrom(1)];
for j=1:length(temp)-1
	len=D(temp(j),temp(j+1))+len;
end
end