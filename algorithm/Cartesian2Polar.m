function polar = Cartesian2Polar(cartesian)
%Cartesian2Polar 

%   �˴���ʾ��ϸ˵��
    polar=cartesian;
    for i=1:size(cartesian,1)
        polar(i,:)=Cartesian_polar(cartesian(i,:),[]);
    end

end

