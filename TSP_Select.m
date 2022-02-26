%% 选择操作
%输入
%Chrom 种群
%FitnV 适应度值
%GGAP：选择概率
%输出
%SelCh  被选择的个体
function SelCh=TSP_Select(Chrom,FitnV,GGAP)
NSel=ceil(size(Chrom,1)*GGAP/2)*2;
ChrIx=TSP_Sus(FitnV,NSel);
SelCh=Chrom(ChrIx,:);