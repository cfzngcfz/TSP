function [S,R]=TSP_SA_Metropolis(S1,S2,D,T)
R1=TSP_SA_PathLength(D,S1);%����·�߳���
R2=TSP_SA_PathLength(D,S2);
N=length(S1);%�õ����и���
dC=R2-R1;%��������֮��
if dC<0
    S=S2;
    R=R2;
elseif exp(-dC/T)>=rand%��һ�����ʽ�����·��
    S=S2;
    R=R2;
else
    S=S1;
    R=R1;
end