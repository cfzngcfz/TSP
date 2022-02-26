%% 计算各个体的路径长度
% 输入：
% D     两两城市之间的距离
% Chrom 个体的轨迹
function len=ObjF(Chrom,D)
len=0;
temp=[Chrom Chrom(1)];
for j=1:length(temp)-1
	len=D(temp(j),temp(j+1))+len;
end
end