function [HC,p] = FCHM(PopDec,Global)
    num = size(PopDec,1);
    n1 = Global.LC;

    W=10;
    [cls,~]=kmeans(PopDec,W);
    epsilon = 0.05;
    bound = zeros(n1,2);
    bound(:,1) = Global.lower(1:n1);
    bound(:,2) = Global.upper(1:n1);
    HC = zeros(n1,W+2);

    PopDec1 = sort(PopDec)'; 
    p1 = PopDec1(:,[1,2,end-1,end]); 
    p = zeros(n1,W+3);
    A = zeros(n1,W+2);
    for n=1:n1 
        p(n,1) = bound(n,1);
        p(n,W+3) = bound(n,2);
        p(n,2) = max(p1(n,1)-0.5*(p1(n,2)-p1(n,1)),bound(n,1));
        p(n,W+2) = min(p1(n,4)+0.5*(p1(n,4)-p1(n,2)),bound(n,2));

        for i = 1:W-1
            bin_width =0.5*(max(PopDec(cls==i,n))+min(PopDec(cls==i+1,n))); 
            p(n,i+2) = p(n,i+1)+bin_width*i;
        end
        for w = 1:W+2
            if w>=2 && w<=W+1
                A(n,w) = length(find(PopDec(:,n)>=p(n,w)&PopDec(:,n)<=p(n,w+1)));
            end
            if (w==1 || w==W+2) && p(n,w+1) > p(n,w)
                A(n,w) = epsilon;
            end
            if (w==1 || w==W+2) && p(n,w+1) == p(n,w)
                A(n,w) = 0;
            end 
            HC(n,w) = A(n,w)/num;
        end
    end
end