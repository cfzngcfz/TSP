function [S,R]=TSP_SA_Metropolis(S1,S2,D,T)
R1=TSP_SA_PathLength(D,S1);%计算路线长度
R2=TSP_SA_PathLength(D,S2);
N=length(S1);%得到城市个数
dC=R2-R1;%计算能量之差
if dC<0
    S=S2;
    R=R2;
elseif exp(-dC/T)>=rand%以一定概率接受新路线
    S=S2;
    R=R2;
else
    S=S1;
    R=R1;
end