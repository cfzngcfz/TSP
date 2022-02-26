clc;
load cities1000_data.mat;

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

%��ʼ������
m=5;%��������
alpha=1;%��Ϣ����Ҫ�̶�����
beta=5;%����������Ҫ�̶�����
rho=0.1;%��Ϣ�ػӷ�����
Q=1;%��ϵ��
Eta=1./D;%��������
Tau=ones(n,n);%��Ϣ�ؾ���
Talbe=zeros(m,n);%·����¼��
iter=1;%����������ֵ
iter_max=50;%����������
Route_best=zeros(iter_max,n);%�������·��
Length_best=zeros(iter_max,1);%�������·���ĳ���
Length_ave=zeros(iter_max,1);%����·����ƽ������

%����Ѱ�����·��
while iter<=iter_max
    %��������������ϵ�������
    start=zeros(m,1);
    for i=1:m
        temp=randperm(n);
        start(i)=temp(1);
    end
    Table(:,1)=start;
    %������ռ�
    cities_index=1:n;
    %������ϵ�·��ѡ��
    for i=1:m
        %�������·��ѡ��
        for j=2:n
            tabu=Table(i,1:(j-1));%�ѷ��ʵĳ��м���(���ɱ�)
            allow_index=~ismember(cities_index,tabu);
            allow=cities_index(allow_index);%�����ʳ��м���
            P=allow;
            %������м�ת�Ƹ���
            for k=1:length(allow)
                P(k)=Tau(tabu(end),allow(k))^alpha*Eta(tabu(end),allow(k))^beta;
            end
            P=P/sum(P);
            %���̶ķ�ѡ����һ�����ʳ���
            Pc=cumsum(P);
            target_index=find(Pc>=rand);
            target=allow(target_index(1));
            Table(i,j)=target;
        end
    end
    %����������ϵ�·������
    Length=zeros(m,1);
    for i=1:m
        Route=Table(i,:);
        for j=1:(n-1)
            Length(i)=Length(i)+D(Route(j),Route(j+1));
        end
        Length(i)=Length(i)+D(Route(n),Route(1));
    end
    %�������·�����뼰ƽ������
    if iter==1
        [min_Length,min_index]=min(Length);
        Length_best(iter)=min_Length;
        Length_ave(iter)=mean(Length);
        Route_best(iter,:)=Table(min_index,:);
    else
        [min_Length,min_index]=min(Length);
        Length_best(iter)=min(Length_best(iter-1),min_Length);
        Length_ave(iter)=mean(Length);
        if Length_best(iter)==min_Length
            Route_best(iter,:)=Table(min_index,:);
        else
            Route_best(iter,:)=Route_best(iter-1,:);
        end
    end
    %������Ϣ��
    Delta_Tau=zeros(n,n);
    %������ϼ���
    for i=1:m
        %������м���
        for j=1:(n-1)
            Delta_Tau(Table(i,j),Table(i,j+1))=Delta_Tau(Table(i,j),Table(i,j+1))+Q/Length(i);
        end
        Delta_Tau(Table(i,n),Table(i,1))=Delta_Tau(Table(i,n),Table(i,1))+Q/Length(i);
    end
    Tau=(1-rho)*Tau+Delta_Tau;
    %����������1�����·����¼��
    iter=iter+1;
    Table=zeros(m,n);
    disp(['��ǰ��������:',num2str(iter)]);
end

%���
[Shortest_Length,index]=min(Length_best);
Shortest_Route=Route_best(index,:);
disp(['��̾���:',num2str(Shortest_Length)]);
disp(['���·��:',num2str([Shortest_Route Shortest_Route(1)])]);
xlswrite('result.xlsx',(Shortest_Route)');