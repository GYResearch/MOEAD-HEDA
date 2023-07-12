function HD = ILH(Pop,Global,HD_pre)
    PopDec = Pop.decs;
    PopObj = Pop.objs;
    W = 10;
    epsilon = 0.01;
    gamma = Global.evaluated/Global.evaluation; 
    num = size(PopDec,1);
    n1 = Global.LD;
    HD = zeros(n1,W);
    A = zeros(n1,W);
    h = zeros(1,num);
    %∑«÷ß≈‰≈≈–Ú
    [Front,~] = NDSort(PopObj,Pop.cons,num);
    for k=1:num
        h(1,k) = (num-Front(k)+1)/num;
    end
    indexes = find(Front <= 3);
    for n=1:n1
        for w = 1:W
            temp = find(PopDec(indexes,n)==w);
            T1 = h(1,temp);
            F = setdiff([1:num],temp);
            T2 = h(1,F);
            A(n,w) = sum(T1*1) + sum(T2*epsilon);
            HD(n,w) = (1-gamma)*HD_pre(n,w)+gamma*A(n,w)/num;
        end
    end
end