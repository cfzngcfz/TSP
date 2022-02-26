clc,clear
close all
%输入数据
load('data1.mat')
% load('data2.mat')

Dist=zeros(size(location,1),size(location,1));
for i=1:size(location,1)
    for j=i+1:size(location,1)
        Dist(i,j)=((location(i,1)-location(j,1))^2+(location(i,2)-location(j,2))^2)^0.5;
        Dist(j,i)=Dist(i,j);
    end
end
%算法通用参数
SizePop=100;             %种群规模
MaxGen=500;            	%迭代次数上限
StopGen=50;             %最优解最少保持代数
%%
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('本程序采用不同的人工智能算法求解TSP问题')
disp('人工智能算法选择：')
fprintf(' 1.遗传算法    2.混合粒子群算法    3.模拟退火算法    4.蚁群算法 \n')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
number_ANN=input('请输入人工智能算法编号=');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if number_ANN==1            %遗传算法

    %算法参数
    GGAP=0.9;               %选择概率（可选参数，默认为1）
    Pc=0.9;                 %交叉概率
    Pm=0.05;                %变异概率
    %算法主体
    for i=1:SizePop
        Chrom(i,:)=randperm(size(Dist,1)); %初始化种群
        ObjV(i,1)=ObjF(Chrom(i,:),Dist); 	%计算初始种群个体的目标函数值
    end
    [bestObjV,bestindex]=min(ObjV);
    bestsol=Chrom(bestindex,:);
    %进化开始
    gen=1;                              %初始遗传代数
    gen0=0;                             %记录保持不变的代数   
    while gen0<StopGen && gen<=MaxGen
        disp(gen)
        FitnV=ranking(ObjV);          	%分配适应度值(Assign fitness values)
        %ranking的输入越小，输出越大
        SelCh=TSP_Select(Chrom,FitnV,GGAP);%选择
        SelCh=TSP_Recombin(SelCh,Pc);   %交叉
        SelCh=TSP_Mutate(SelCh,Pm);     %变异
        SelCh=TSP_Reverse(SelCh,Dist);     %逆转
        Chrom=TSP_Reins(Chrom,SelCh,ObjV);%重插入子代的新种群
         for i=1:SizePop
            ObjV(i,1)=ObjF(Chrom(i,:),Dist);%计算初始种群个体的目标函数值
        end
        [newbestObjV,newbestindex]=min(ObjV);
        [worestObjV,worestindex]=max(ObjV);
        if newbestObjV<bestObjV
            bestObjV=newbestObjV;
            bestsol=Chrom(newbestindex,:);
            gen0=0;
        else
            gen0=gen0+1;                %最优值保持次数加1
        end
        Chrom(worestindex,:)=bestsol;
        ObjV(worestindex,:)=bestObjV;
        trace(gen,1)=bestObjV;         	%记录每代的最优值
        trace(gen,2)=sum(ObjV)/SizePop; %记录每一代进化中平均适应度
        gen=gen+1;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif number_ANN==2      %混合粒子群算法
    
    num=size(location,1);                 %计算位置数量
    for i=1:SizePop
        individual(i,:)=randperm(num);    %生成初始解
    end
    individual_ObjV=ones(SizePop,1)*Inf;%记录个体最优值
    bestObjV=Inf;                       %记录群体最优值
    gen=1;                              %初始遗传代数
    gen0=0;                             %记录保持不变的代数   
    while gen0<StopGen && gen<=MaxGen
        disp(gen)
        newindividual=individual;
        %计算适应度值
        for i=1:SizePop
            ObjV(i,:)=ObjF(individual(i,:),Dist);
        end
        %更新当前最优和历史最优
        for i=1:SizePop
            if ObjV(i)<individual_ObjV(i)
                individual_ObjV(i)=ObjV(i);
                individual_sol(i,:)=individual(i,:);
            end
        end
        [newbestObjV,index]=min(individual_ObjV);
        if newbestObjV<bestObjV
            bestObjV=newbestObjV;
            gen0=0;
            bestsol=individual_sol(index,:);
        else 
            gen0=gen0+1;
        end
        trace(gen,1)=bestObjV;
        trace(gen,2)=mean(individual_ObjV);
        %交叉和变异操作
        for i=1:SizePop
            %与个体最优进行交叉
            temp1=randperm(num);
            c1=temp1(1);
            c2=temp1(2);
            cros=individual_sol(i,min(c1,c2):max(c1,c2));
            %删除与交叉区域相同元素
            temp2=newindividual(i,:);
            for j=1:size(cros,2)
                temp2(temp2==cros(j))=[];
            end
            newindividual(i,:)=[temp2 cros];
            %新路径长度变短则接受
            newObjV=ObjF(newindividual(i,:),Dist);
            if newObjV<ObjV(i)
                individual(i,:)=newindividual(i,:);
                ObjV(i)=newObjV;
            end
            % 与全体最优进行交叉
            temp3=randperm(num);
            c1=temp3(1);
            c2=temp3(2);
            cros=bestsol(min(c1,c2):max(c1,c2));
            %删除与交叉区域相同元素
            temp4=newindividual(i,:);
            for j=1:size(cros,2)
                temp4(temp4==cros(j))=[];
            end
            newindividual(i,:)=[temp4 cros];
            %新路径长度变短则接受
            newObjV=ObjF(newindividual(i,:),Dist);
            if newObjV<ObjV(i)
                individual(i,:)=newindividual(i,:);
                ObjV(i)=newObjV;
            end
            %变异操作
            temp5=randperm(num);
            c1=temp5(1);
            c2=temp5(2); 
            temp6=newindividual(i,c1);
            newindividual(i,c1)=newindividual(i,c2);
            newindividual(i,c2)=temp6;
            %新路径长度变短则接受
            newObjV=ObjF(newindividual(i,:),Dist);
            if newObjV<ObjV(i)
                individual(i,:)=newindividual(i,:);
                ObjV(i)=newObjV;
            end
        end
        gen=gen+1;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif number_ANN==3      %模拟退火算法
    
    %模拟退火算法参数
    T0=1000;   %初始温度
    Tend=1e-6;  %终止温度
    q=0.97;    %降温速率
    num=size(location,1);    %城市的个数
    %生成初始解
    S1=randperm(num);
    bestObjV=ObjF(S1,Dist);                    	%记录最优值
    bestsol=S1;                                 %记录最优解
    %开始迭代
    gen=1;                                      %记录迭代次数
    gen2=0;                                     %记录温度下降过程中目标函数值不变的次数  
    while T0>Tend && gen2<StopGen*2
        disp(gen)                               %输出当前迭代次数
        temp=[];
        gen3=1;
        S2=TSP_NewAnswer(S1);
        [S1,R]=TSP_Metropolis(S1,S2,Dist,T0);  %Metropolis 抽样算法
        temp(gen3,:)=[S1 R];
        gen4=0;                                 %记录固定温度下目标函数值不变的次数
        while  gen3<=MaxGen && gen4<StopGen*5
            gen3=gen3+1;
            S2=TSP_NewAnswer(S1);%产生新解
            [S1,R]=TSP_Metropolis(S1,S2,Dist,T0);%Metropolis法则判断是否接受新解
            temp(gen3,:)=[S1 R];          %记录下一路线的及其路程
            if temp(gen3,end)<temp(gen3-1,end)
                gen4=0;
            else
                gen4=gen4+1;
            end
        end
        %记录每次迭代过程的最优路线
        [newbestObjV,index]=min(temp(:,end)); 	%找出当前温度下最优值
        if newbestObjV<bestObjV
            bestObjV=newbestObjV;             	%如果当前温度下最优值小于上一温度下最优值，则记录当前最优值
            bestsol=temp(index,1:end-1);
            gen2=0;
        else
            gen2=gen2+1;
        end
        trace(gen,1)=bestObjV;                	%记录当前温度的最优值
        T0=q*T0;                               	%降温
        gen=gen+1;                              %更新迭代次数
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif number_ANN==4      %蚁群算法
    
    %初始化参数
    alpha = 1;                                  %信息素重要程度因子
    beta = 5;                                  	%启发函数重要程度因子
    rho = 0.1;                                 	%信息素挥发因子
    Q = 1;                                     	%常系数
    Eta = 1./Dist;                             	%启发函数
    Tau = ones(size(location,1));             	%信息素矩阵
    %算法主体
	gen=1;                            % 迭代次数初值
    disp(gen)
    Table=zeros(SizePop,size(location,1));       	% 路径记录表
	for i=1:SizePop
        Table(i,1)=randperm(size(location,1),1);
	end
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
        ObjV(i)=ObjF(Table(i,:),Dist);
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
        for i=1:SizePop
            Table(i,1)=randperm(size(location,1),1);
        end
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
            ObjV(i)=ObjF(Table(i,:),Dist);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    disp('！！！人工智能算法编号超出范围！！！')
end
%%
%结果绘图
%图1：优化过程图
if size(trace,2)==1
    plot(1:size(trace,1),trace(:,1));
    legend('各代最优值');
elseif size(trace,2)==2
    plot(1:size(trace,1),trace(:,1),'r-',1:size(trace,1),trace(:,2),'b--');
    legend('各代最优值','各代平均值');
end
title(['目标函数值优化曲线  ' '终止代数＝' num2str(size(trace,1))],'fontsize',12);
xlabel('进化代数','fontsize',12);
ylabel('目标函数值','fontsize',12);
xlim([1,size(trace,1)])
%图2：路线规划图
figure;
hold on
plot(location(:,1),location(:,2),'r.', 'MarkerSize', 20)
for i=1:size(location,1)
    text(location(i,1)+0.1,location(i,2)+0.1,num2str(i),'color',[1,0,0]);
end
plot(location(bestsol(1),1),location(bestsol(1),2),'rp','MarkerSize',20)
R=[bestsol bestsol(1)]; %一个随机解(个体)
A=location(R,:);
for i=2:size(A,1)
    [arrowx,arrowy] = dsxy2figxy(gca,A(i-1:i,1),A(i-1:i,2));%坐标转换
    annotation('textarrow',arrowx,arrowy,'HeadWidth',8,'color',[0,0,1]);
end
hold off
xlabel('横坐标')
ylabel('纵坐标')
title(['优化路径(最短总距离:' num2str(bestObjV) ')'])
box on
grid on
%% 
%结果输出
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('TSP最优路线：')
R=[bestsol,bestsol(1)];
p=num2str(R(1));
for i=2:length(R)
    p=[p,'->',num2str(R(i))];
end
disp(p)
disp(['最短总距离：',num2str(bestObjV)]);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')