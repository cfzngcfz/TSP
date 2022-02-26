clc;
load cities3795_data.mat;

%初始化参数
T_start=1000;%初始温度
T_end=1e-3;%终止温度
L=200;%各温度下的迭代次数(链长)
q=0.9;%降温速率

%计算城市间相互距离
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

t=ceil(double(solve(['1000*(0.9)^x=',num2str(T_end)])));%计算迭代次数
count=0;%迭代计数
Obj=zeros(t,1);%目标值矩阵初始化
track=zeros(t,n);%每代的最优路线矩阵初始化
S1=randperm(n);%随机产生一个初始路线

%迭代
while T_start>T_end
    count=count+1;%更新迭代次数
    temp=zeros(L,n+1);
    for k=1:L
        S2=TSP_SA_NewAnswer(S1);%产生新解
        [S1,R]=TSP_SA_Metropolis(S1,S2,D,T_start);%Metropolis法则判断是否接受新解
        temp(k,:)=[S1,R];%记录下一路线及其路程
    end
    
    %记录每次迭代过程的最优路线
    [d0,index]=min(temp(:,end));%找出当前温度下的最优路线
    if count==1 || d0<Obj(count-1)
        Obj(count)=d0;%如果当前温度下最优路程小于上一路程，则记录当前路程
    else
        Obj(count)=Obj(count-1);%如果当前温度下最优路程大于上一路程，则记录上一路程
    end
    track(count,:)=temp(index,1:end-1);%记录当前温度的最优路线
    T_start=0.9*T_start;%降温
    disp(['当前迭代次数:',num2str(count)]);
end

%输出结果
%disp('最优解:')
S=track(end,:);

%输出路径
S=[S,S(1)];
N=length(S);
p=num2str(S(1));
for i=2:N
    p=[p,'->',num2str(S(i))];
end
disp(['最优解:',p]);

%p=OutputPath(S);
disp(['总距离:',num2str(TSP_SA_PathLength(D,S))]);
xlswrite('result.xlsx',(S)');

