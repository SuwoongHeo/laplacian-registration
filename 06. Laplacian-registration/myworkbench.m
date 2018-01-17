landmark = zeros(17,1);
for i = 1:5
    a = find(VT(:,1)==landmark3(i).Position(1,1) & VT(:,2)==landmark3(i).Position(1,2) & VT(:,3)==landmark3(i).Position(1,3));
    landmark = [landmark ; a(1,1)];
end
for i = 1:8
    a = find(VS(:,1)==FW_land2(i).Position(1,1) & VS(:,2)==FW_land2(i).Position(1,2) & VS(:,3)==FW_land2(i).Position(1,3));
    landmark(i+9,1) = a(1,1);
end

a = [0 1 0];
b = [1 0 0];
theta = atan2d(norm(cross(a,b)),dot(a,b));
a = find(atan2d(norm(cross(new_NS(i,:),new_NT)),dot(new_NS(i,:),new_NT)) < 45);

sqrt(sum((VS(marker(1,1),:)-VT(marker(1,2),:)).*(VS(marker(1,1),:)-VT(marker(1,2),:))));
i = marker(1,1);
j = marker(1,2);
test = atan2d(norm(cross(new_NS(i,:),new_NT(j,:))),dot(new_NS(i,:),new_NT(j,:)));