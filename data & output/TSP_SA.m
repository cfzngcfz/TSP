clc;
load cities3795_data.mat;

%��ʼ������
T_start=1000;%��ʼ�¶�
T_end=1e-3;%��ֹ�¶�
L=200;%���¶��µĵ�������(����)
q=0.9;%��������

%������м��໥����
n=size(cities,1);
D=zeros(n,n);
for i=1:n
    for j=1:n
        if i~=j
            D(i,j)=sqrt(sum((cities(i,:)-cities(j,:)).^2));
        else
            D(i,j)=1e-4;
        end
    end
end

t=ceil(double(solve(['1000*(0.9)^x=',num2str(T_end)])));%�����������
count=0;%��������
Obj=zeros(t,1);%Ŀ��ֵ�����ʼ��
track=zeros(t,n);%ÿ��������·�߾����ʼ��
S1=randperm(n);%�������һ����ʼ·��

%����
while T_start>T_end
    count=count+1;%���µ�������
    temp=zeros(L,n+1);
    for k=1:L
        S2=TSP_SA_NewAnswer(S1);%�����½�
        [S1,R]=TSP_SA_Metropolis(S1,S2,D,T_start);%Metropolis�����ж��Ƿ�����½�
        temp(k,:)=[S1,R];%��¼��һ·�߼���·��
    end
    
    %��¼ÿ�ε������̵�����·��
    [d0,index]=min(temp(:,end));%�ҳ���ǰ�¶��µ�����·��
    if count==1 || d0<Obj(count-1)
        Obj(count)=d0;%�����ǰ�¶�������·��С����һ·�̣����¼��ǰ·��
    else
        Obj(count)=Obj(count-1);%�����ǰ�¶�������·�̴�����һ·�̣����¼��һ·��
    end
    track(count,:)=temp(index,1:end-1);%��¼��ǰ�¶ȵ�����·��
    T_start=0.9*T_start;%����
    disp(['��ǰ��������:',num2str(count)]);
end

%������
%disp('���Ž�:')
S=track(end,:);

%���·��
S=[S,S(1)];
N=length(S);
p=num2str(S(1));
for i=2:N
    p=[p,'->',num2str(S(i))];
end
disp(['���Ž�:',p]);

%p=OutputPath(S);
disp(['�ܾ���:',num2str(TSP_SA_PathLength(D,S))]);
xlswrite('result.xlsx',(S)');

