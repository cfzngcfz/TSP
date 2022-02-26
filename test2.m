%% 第22章 蚁群算法的优化计算——旅行商问题(TSP)优化

%% 清空环境变量
clear
clc

% 导入数据
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
    20.09,92.55];%个城市坐标位置

% 计算城市间相互距离
Dist=zeros(size(location,1),size(location,1));
for i=1:size(location,1)
    for j=i+1:size(location,1)
        Dist(i,j)=((location(i,1)-location(j,1))^2+(location(i,2)-location(j,2))^2)^0.5;
        Dist(j,i)=Dist(i,j);
    end
end
SizePop = 50;                              % 蚂蚁数量
MaxGen = 200;                      % 最大迭代次数 
StopGen=50;             %最优解最少保持代数
%% 初始化参数

alpha = 1;                                      % 信息素重要程度因子
beta = 5;                                       % 启发函数重要程度因子
rho = 0.1;                                      % 信息素挥发因子
Q = 1;                                          % 常系数
Eta = 1./Dist;                                  % 启发函数
Tau = ones(size(location,1));                   % 信息素矩阵



% Route_best = zeros(MaxGen,size(location,1));      % 各代最佳路径       
% Length_best = zeros(MaxGen,1);     % 各代最佳路径的长度  
% Length_ave = zeros(MaxGen,1);      % 各代路径的平均长度  





    %算法主体
	gen = 1;                            % 迭代次数初值
    Table = zeros(SizePop,size(location,1));       	% 路径记录表
    start=zeros(SizePop,1);
	for i=1:SizePop
        start(i,1)=randperm(size(location,1),1);
	end
	Table(:,1)=start; 
	% 构建解空间
	location_index=1:size(location,1);
	% 逐个蚂蚁路径选择
	for i=1:SizePop
        % 逐个城市路径选择
        for j=2:size(location,1)
            tabu=Table(i,1:(j - 1));           % 已访问的城市集合(禁忌表)
            allow=location_index(~ismember(location_index,tabu));  % 待访问的城市集合
            P=ones(1,length(allow));
            % 计算城市间转移概率
            for k=1:length(allow)
                P(k)=Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;%转移概率的分子
            end
            P=P/sum(P);%转移概率
            % 轮盘赌法选择下一个访问城市
            Pc=cumsum(P);     
            Table(i,j)=allow(find(Pc>=rand,1));
        end
	end
	% 计算各个蚂蚁的路径距离
	ObjV=zeros(SizePop,1);
	for i=1:SizePop
        Route=Table(i,:);
        for j=1:(size(location,1) - 1)
            ObjV(i)=ObjV(i)+Dist(Route(j),Route(j + 1));
        end
        ObjV(i)=ObjV(i)+Dist(Route(size(location,1)),Route(1));
	end
	% 更新信息素
	Delta_Tau=zeros(size(location,1));
	% 逐个蚂蚁计算
	for i=1:SizePop
        % 逐个城市计算
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
    % 迭代开始
    gen0=0;                             %记录保持不变的代数   
    while gen0<StopGen && gen<=MaxGen
        gen=gen+1;
        disp(gen)
        Table=zeros(SizePop,size(location,1));
        % 随机产生各个蚂蚁的起点城市
        start=zeros(SizePop,1);
        for i=1:SizePop
            start(i,1)=randperm(size(location,1),1);
        end
        Table(:,1)=start; 
        % 构建解空间
        location_index=1:size(location,1);
        % 逐个蚂蚁路径选择
        for i=1:SizePop
            % 逐个城市路径选择
            for j=2:size(location,1)
                tabu=Table(i,1:(j - 1));           % 已访问的城市集合(禁忌表)
                allow=location_index(~ismember(location_index,tabu));  % 待访问的城市集合
                P=ones(1,length(allow));
                % 计算城市间转移概率
                for k=1:length(allow)
                    P(k)=Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;%转移概率的分子
                end
                P=P/sum(P);%转移概率
                % 轮盘赌法选择下一个访问城市
                Pc=cumsum(P);     
                Table(i,j)=allow(find(Pc>=rand,1));
            end
        end
        % 计算各个蚂蚁的路径距离
        ObjV=zeros(SizePop,1);
        for i=1:SizePop
            Route=Table(i,:);
            for j=1:(size(location,1) - 1)
                ObjV(i)=ObjV(i)+Dist(Route(j),Route(j + 1));
            end
            ObjV(i)=ObjV(i)+Dist(Route(size(location,1)),Route(1));
        end
        % 更新信息素
        Delta_Tau=zeros(size(location,1));
        % 逐个蚂蚁计算
        for i=1:SizePop
            % 逐个城市计算
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
            gen0=gen0+1;                %最优值保持次数加1
        end
        trace(gen,1)=bestObjV;         	%记录每代的最优值
        trace(gen,2)=sum(ObjV)/SizePop; %记录每一代进化中平均适应度
    end

% %% 结果显示
% [Shortest_Length,index] = min(Length_best);
% Shortest_Route = Route_best(index,:);
% disp(['最短距离:' num2str(Shortest_Length)]);
% disp(['最短路径:' num2str([Shortest_Route Shortest_Route(1)])]);
% 
% %% 绘图
% figure(1)
% plot([location(Shortest_Route,1);location(Shortest_Route(1),1)],...
%      [location(Shortest_Route,2);location(Shortest_Route(1),2)],'o-');
% grid on
% for i = 1:size(location,1)
%     text(location(i,1),location(i,2),['   ' num2str(i)]);
% end
% text(location(Shortest_Route(1),1),location(Shortest_Route(1),2),'       起点');
% text(location(Shortest_Route(end),1),location(Shortest_Route(end),2),'       终点');
% xlabel('城市位置横坐标')
% ylabel('城市位置纵坐标')
% title(['蚁群算法优化路径(最短距离:' num2str(Shortest_Length) ')'])
% figure(2)
% plot(1:MaxGen,Length_best,'b',1:MaxGen,Length_ave,'r:')
% legend('最短距离','平均距离')
% xlabel('迭代次数')
% ylabel('距离')
% title('各代最短距离与平均距离对比')