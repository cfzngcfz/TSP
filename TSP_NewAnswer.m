function S2=TSP_NewAnswer(S1)
% ����
% S1:��ǰ��
% ���
% S2���½�
S2=S1;
a=randperm(length(S1));%�����������λ�� ��������
temp=S2(a(1));
S2(a(1))=S2(a(2));
S2(a(2))=temp;         %�õ�һ����·��
end