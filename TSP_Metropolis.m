function [S,R]=TSP_Metropolis(S1,S2,Dist,T)
% ����
% S1��  ��ǰ��
% S2:   �½�
% D:    ��������������е�֮��ľ��룩
% T:    ��ǰ�¶�
% ���
% S��   ��һ����ǰ��
% R��   ��һ����ǰ���·�߾���
R1=ObjF(S1,Dist);           %����·�߳���
R2=ObjF(S2,Dist);           %����·�߳���
if R2<R1                    %����������� ������·��
    S=S2;
    R=R2;
elseif exp((R1-R2)/T)>=rand	%��exp(-dC/T)���ʽ�����·��
    S=S2;
    R=R2;
else                        %��������·��
    S=S1;
    R=R1;
end
end