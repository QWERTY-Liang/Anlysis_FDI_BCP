%data为输入数据；m为模板维数，通常取2,3,4；r为模板匹配相似容忍度，通常取0.1-0.25倍的信号标准差
function output = TL6_FuzzyEn(data,m,r)
N = length(data);
% m = 2;
n = 2;
% r = 0.25;
r = r*std(data);
%% 模版构成
for i = 1:N-m+1
    XX(i,:) = data(i:i+m-1);%构造m维模板，模板总数为N-m+1,每个模板包含m个数
    u0(i,1:m) = sum(XX(i,:))/m;%计算第i个模板m个点的均值
end
X = XX - u0;%用每个模板减去自身的均值，得到去均值后的模板
%% 计算距离
for i = 1:N-m+1
    data2 = X(i,:);%选择第i个模板
    data22 = X;
    data22(i,:) = [];%将m维模板中的第i个模板去掉
    d2 = [ones(N-m,1) * data2] - data22;%计算第i个模板到其他N-m个模板的对应点的距离，即消除了自身模板匹配
    for k = 1:N-m
        D2(k) = exp(-(max(abs(d2(k,:)))/r)^n);%用模糊函数计算两个模板之间的相似值，n通常取2
    end
    c2(i) = sum(D2)/(N-m);%用所有的模板相似值的和除以模板匹配的数目即得相似模板的比例
    clear data22
end
%% 计算平均相似率
y2 = sum(c2)/(N-m+1);


%% 模版构成
m = m + 1;
for i = 1:N-m+1
    YY(i,:) = data(i:i+m-1);
    v0(i,1:m) = sum(YY(i,:))/m;
end
Y = YY - v0;
%% 计算距离
for i = 1:N-m+1
    data3 = Y(i,:);
    data33 = Y;
    data33(i,:) = [];
    d3 = [ones(N-m,1) * data3] - data33;
    for k = 1:N-m
        D3(k) = exp(-(max(abs(d3(k,:)))/r)^n);
    end
    c3(i) = sum(D3)/(N-m);
    clear data33
end
y3 = sum(c3)/(N-m+1);
%% 计算模糊熵近似熵
output = -log(y3/y2);%y2为m维模板的平均相似率，y3为m+1维的平均相似率