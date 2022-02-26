%% ��22�� ��Ⱥ�㷨���Ż����㡪������������(TSP)�Ż�

%% ��ջ�������
clear
clc

% ��������
location=[16.47,96.10
    16.47,94.44
    18.09,92.54
    22.39,93.37
    25.23,97.24
    22.00,96.05
    20.47,97.02
    17.20,96.29
    16.30,97.38
    14.05,98.12
    18.53,98.38
    21.52,95.59
    19.41,97.13
    20.09,92.55];%����������λ��

% ������м��໥����
Dist=zeros(size(location,1),size(location,1));
for i=1:size(location,1)
    for j=i+1:size(location,1)
        Dist(i,j)=((location(i,1)-location(j,1))^2+(location(i,2)-location(j,2))^2)^0.5;
        Dist(j,i)=Dist(i,j);
    end
end
SizePop = 50;                              % ��������
MaxGen = 200;                      % ���������� 
StopGen=50;             %���Ž����ٱ��ִ���
%% ��ʼ������

alpha = 1;                                      % ��Ϣ����Ҫ�̶�����
beta = 5;                                       % ����������Ҫ�̶�����
rho = 0.1;                                      % ��Ϣ�ػӷ�����
Q = 1;                                          % ��ϵ��
Eta = 1./Dist;                                  % ��������
Tau = ones(size(location,1));                   % ��Ϣ�ؾ���



% Route_best = zeros(MaxGen,size(location,1));      % �������·��       
% Length_best = zeros(MaxGen,1);     % �������·���ĳ���  
% Length_ave = zeros(MaxGen,1);      % ����·����ƽ������  





    %�㷨����
	gen = 1;                            % ����������ֵ
    Table = zeros(SizePop,size(location,1));       	% ·����¼��
    start=zeros(SizePop,1);
	for i=1:SizePop
        start(i,1)=randperm(size(location,1),1);
	end
	Table(:,1)=start; 
	% ������ռ�
	location_index=1:size(location,1);
	% �������·��ѡ��
	for i=1:SizePop
        % �������·��ѡ��
        for j=2:size(location,1)
            tabu=Table(i,1:(j - 1));           % �ѷ��ʵĳ��м���(���ɱ�)
            allow=location_index(~ismember(location_index,tabu));  % �����ʵĳ��м���
            P=ones(1,length(allow));
            % ������м�ת�Ƹ���
            for k=1:length(allow)
                P(k)=Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;%ת�Ƹ��ʵķ���
            end
            P=P/sum(P);%ת�Ƹ���
            % ���̶ķ�ѡ����һ�����ʳ���
            Pc=cumsum(P);     
            Table(i,j)=allow(find(Pc>=rand,1));
        end
	end
	% ����������ϵ�·������
	ObjV=zeros(SizePop,1);
	for i=1:SizePop
        Route=Table(i,:);
        for j=1:(size(location,1) - 1)
            ObjV(i)=ObjV(i)+Dist(Route(j),Route(j + 1));
        end
        ObjV(i)=ObjV(i)+Dist(Route(size(location,1)),Route(1));
	end
	% ������Ϣ��
	Delta_Tau=zeros(size(location,1));
	% ������ϼ���
	for i=1:SizePop
        % ������м���
        for j=1:size(location,1)-1
            Delta_Tau(Table(i,j),Table(i,j+1))=Delta_Tau(Table(i,j),Table(i,j+1))+Q/ObjV(i);
        end
        Delta_Tau(Table(i,size(location,1)),Table(i,1))=Delta_Tau(Table(i,size(location,1)),Table(i,1))+Q/ObjV(i);
	end
	Tau=(1-rho)*Tau+Delta_Tau;
    [bestObjV,bestindex]=min(ObjV);
    trace(gen,1)=bestObjV;  
    trace(gen,2)=mean(ObjV);
    bestsol=Table(bestindex,:);
    % ������ʼ
    gen0=0;                             %��¼���ֲ���Ĵ���   
    while gen0<StopGen && gen<=MaxGen
        gen=gen+1;
        disp(gen)
        Table=zeros(SizePop,size(location,1));
        % ��������������ϵ�������
        start=zeros(SizePop,1);
        for i=1:SizePop
            start(i,1)=randperm(size(location,1),1);
        end
        Table(:,1)=start; 
        % ������ռ�
        location_index=1:size(location,1);
        % �������·��ѡ��
        for i=1:SizePop
            % �������·��ѡ��
            for j=2:size(location,1)
                tabu=Table(i,1:(j - 1));           % �ѷ��ʵĳ��м���(���ɱ�)
                allow=location_index(~ismember(location_index,tabu));  % �����ʵĳ��м���
                P=ones(1,length(allow));
                % ������м�ת�Ƹ���
                for k=1:length(allow)
                    P(k)=Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;%ת�Ƹ��ʵķ���
                end
                P=P/sum(P);%ת�Ƹ���
                % ���̶ķ�ѡ����һ�����ʳ���
                Pc=cumsum(P);     
                Table(i,j)=allow(find(Pc>=rand,1));
            end
        end
        % ����������ϵ�·������
        ObjV=zeros(SizePop,1);
        for i=1:SizePop
            Route=Table(i,:);
            for j=1:(size(location,1) - 1)
                ObjV(i)=ObjV(i)+Dist(Route(j),Route(j + 1));
            end
            ObjV(i)=ObjV(i)+Dist(Route(size(location,1)),Route(1));
        end
        % ������Ϣ��
        Delta_Tau=zeros(size(location,1));
        % ������ϼ���
        for i=1:SizePop
            % ������м���
            for j=1:size(location,1)-1
                Delta_Tau(Table(i,j),Table(i,j+1))=Delta_Tau(Table(i,j),Table(i,j+1))+Q/ObjV(i);
            end
            Delta_Tau(Table(i,size(location,1)),Table(i,1))=Delta_Tau(Table(i,size(location,1)),Table(i,1))+Q/ObjV(i);
        end
        Tau=(1-rho)*Tau+Delta_Tau;
        [newbestObjV,newbestindex]=min(ObjV);
        if newbestObjV<bestObjV
            bestObjV=newbestObjV;
            bestsol=Table(newbestindex,:);
            gen0=0;
        else
            gen0=gen0+1;                %����ֵ���ִ�����1
        end
        trace(gen,1)=bestObjV;         	%��¼ÿ��������ֵ
        trace(gen,2)=sum(ObjV)/SizePop; %��¼ÿһ��������ƽ����Ӧ��
    end

% %% �����ʾ
% [Shortest_Length,index] = min(Length_best);
% Shortest_Route = Route_best(index,:);
% disp(['��̾���:' num2str(Shortest_Length)]);
% disp(['���·��:' num2str([Shortest_Route Shortest_Route(1)])]);
% 
% %% ��ͼ
% figure(1)
% plot([location(Shortest_Route,1);location(Shortest_Route(1),1)],...
%      [location(Shortest_Route,2);location(Shortest_Route(1),2)],'o-');
% grid on
% for i = 1:size(location,1)
%     text(location(i,1),location(i,2),['   ' num2str(i)]);
% end
% text(location(Shortest_Route(1),1),location(Shortest_Route(1),2),'       ���');
% text(location(Shortest_Route(end),1),location(Shortest_Route(end),2),'       �յ�');
% xlabel('����λ�ú�����')
% ylabel('����λ��������')
% title(['��Ⱥ�㷨�Ż�·��(��̾���:' num2str(Shortest_Length) ')'])
% figure(2)
% plot(1:MaxGen,Length_best,'b',1:MaxGen,Length_ave,'r:')
% legend('��̾���','ƽ������')
% xlabel('��������')
% ylabel('����')
% title('������̾�����ƽ������Ա�')