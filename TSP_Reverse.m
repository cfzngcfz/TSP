%% ������ת����
%����
%SelCh ��ѡ��ĸ���
%D     �����еľ������
%���
%SelCh  ������ת��ĸ���
function SelCh=TSP_Reverse(SelCh,D)
[row,col]=size(SelCh);
for i=1:size(SelCh,1)
    ObjV(i,1)=ObjF(SelCh(i,:),D);      %�����ʼ��Ⱥ�����Ŀ�꺯��ֵ
end
% ObjV=ObjF(D,SelCh);  %����·������
SelCh1=SelCh;
for i=1:row
    r1=randsrc(1,1,[1:col]);
    r2=randsrc(1,1,[1:col]);
    mininverse=min([r1 r2]);
    maxinverse=max([r1 r2]);
    SelCh1(i,mininverse:maxinverse)=SelCh1(i,maxinverse:-1:mininverse);
end
for i=1:size(SelCh,1)
    ObjV1(i,1)=ObjF(SelCh1(i,:),D);      %�����ʼ��Ⱥ�����Ŀ�꺯��ֵ
end
% ObjV1=ObjF(D,SelCh1);  %����·������
index=ObjV1<ObjV;
SelCh(index,:)=SelCh1(index,:);