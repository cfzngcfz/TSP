function S2=TSP_SA_NewAnswer(S1)
N=length(S1);
S2=S1;
a=round(rand(1,2)*(N-1)+1);%���������������λ��
W=S2(a(1));
S2(a(1))=S2(a(2));
S2(a(2))=W;
