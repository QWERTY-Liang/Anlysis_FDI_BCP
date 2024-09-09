%测试 结束时10 到应答时的时间差（基本为232）
for i=1:6016
    if size(EEG.epoch(i).eventlatency,2)==6
    test(i)=EEG.epoch(i).eventlatency{1,6}-EEG.epoch(i).eventlatency{1,5};
    elseif size(EEG.epoch(i).eventlatency,2)==5
    test(i)= 0;
    end

end

histogram(test)

%测试 图片时101-104 到开始时间差
for i=1:6016

    test(i)=EEG.epoch(i).eventlatency{1,4}-EEG.epoch(i).eventlatency{1,2};


end

histogram(test)

%%
