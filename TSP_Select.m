%% ѡ�����
%����
%Chrom ��Ⱥ
%FitnV ��Ӧ��ֵ
%GGAP��ѡ�����
%���
%SelCh  ��ѡ��ĸ���
function SelCh=TSP_Select(Chrom,FitnV,GGAP)
NSel=ceil(size(Chrom,1)*GGAP/2)*2;
ChrIx=TSP_Sus(FitnV,NSel);
SelCh=Chrom(ChrIx,:);