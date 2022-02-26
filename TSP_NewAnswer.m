function S2=TSP_NewAnswer(S1)
% 输入
% S1:当前解
% 输出
% S2：新解
S2=S1;
a=randperm(length(S1));%产生两个随机位置 用来交换
temp=S2(a(1));
S2(a(1))=S2(a(2));
S2(a(2))=temp;         %得到一个新路线
end