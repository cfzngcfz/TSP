%% 进化逆转函数
%输入
%SelCh 被选择的个体
%D     个城市的距离矩阵
%输出
%SelCh  进化逆转后的个体
function SelCh=TSP_Reverse(SelCh,D)
[row,col]=size(SelCh);
for i=1:size(SelCh,1)
    ObjV(i,1)=ObjF(SelCh(i,:),D);      %计算初始种群个体的目标函数值
end
% ObjV=ObjF(D,SelCh);  %计算路径长度
SelCh1=SelCh;
for i=1:row
    r1=randsrc(1,1,[1:col]);
    r2=randsrc(1,1,[1:col]);
    mininverse=min([r1 r2]);
    maxinverse=max([r1 r2]);
    SelCh1(i,mininverse:maxinverse)=SelCh1(i,maxinverse:-1:mininverse);
end
for i=1:size(SelCh,1)
    ObjV1(i,1)=ObjF(SelCh1(i,:),D);      %计算初始种群个体的目标函数值
end
% ObjV1=ObjF(D,SelCh1);  %计算路径长度
index=ObjV1<ObjV;
SelCh(index,:)=SelCh1(index,:);