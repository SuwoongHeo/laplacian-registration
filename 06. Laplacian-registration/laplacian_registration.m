function x = laplacian_registration(VS, FS, NS, VT, FT, NT, marker, wl, wc, dist, theta)
% VS : Source vertex
% VT : Target vertex
% FS : Source face index
% FT : Target face index
% NS : Source normal
% NT : Target normal
% marker : [source landmark ; target landmark]
% wl : laplacian weight
% wc : correspondence weight
% dist : distance constraint of each vertex in corresponding process
% theta : angle constraint of normal in corresponding process

% normalize vertex
tmean = 0;
tstd = sqrt(2);
VS = normPts(VS, tmean, tstd);
VT = normPts(VT, tmean, tstd);
figure, dispMesh(VS, FS, [.8 .8 .8],.3);
hold on, dispMesh(VT, FT, [0 0 .8], .5);


% find coordinate by PCA
mean_VS = sum(VS) / length(VS);
mean_VT = sum(VT) / length(VT);

VS_p = VS - repmat(mean_VS, [length(VS),1]);
VT_p = VT - repmat(mean_VT, [length(VT),1]);

Cov_VS = VS_p'*VS_p;
Cov_VT = VT_p'*VT_p;

[vector_VS value_VS] = eig(Cov_VS);
[vector_VT value_VT] = eig(Cov_VT);

% vector_VS = [vector_VS(:,1),vector_VS(:,2),vector_VS(:,3)];
vector_VT = [vector_VT(:,2),vector_VT(:,3),vector_VT(:,1)];

VS = (vector_VS'*VS')';
VT = (vector_VT'*VT')';

figure, dispMesh(VS, FS, [.8 .8 .8],.3);
figure, dispMesh(VT, FT, [.8 .8 .8],.3);

figure, dispMesh(VS, FS, [.8 .8 .8],.3);
hold on, quiver3(mean_VS(1,1), mean_VS(1,2), mean_VS(1,3), vector_VS(1,1), vector_VS(2,1), vector_VS(3,1));
quiver3(mean_VS(1,1), mean_VS(1,2), mean_VS(1,3), vector_VS(1,2), vector_VS(2,2), vector_VS(3,2));
quiver3(mean_VS(1,1), mean_VS(1,2), mean_VS(1,3), vector_VS(1,3), vector_VS(2,3), vector_VS(3,3));

figure, dispMesh(VT, FT, [.8 .8 .8],.3);
hold on, quiver3(mean_VT(1,1), mean_VT(1,2), mean_VT(1,3), vector_VT(1,1), vector_VT(2,1), vector_VT(3,1));
quiver3(mean_VT(1,1), mean_VT(1,2), mean_VT(1,3), vector_VT(1,2), vector_VT(2,2), vector_VT(3,2));
quiver3(mean_VT(1,1), mean_VT(1,2), mean_VT(1,3), vector_VT(1,3), vector_VT(2,3), vector_VT(3,3));




[R, t, s, res] = similarity_fitting(VT(marker(:,2),:), VS(marker(:,1),:));
VT = (VT*(diag(s)*R)' + repmat(t, length(VT), 1));

figure, dispMesh(VS, FS, [.8 .8 .8],.3);
hold on, dispMesh(VT, FT, [0 0 .8], .5);
scatter3(VS(marker(:,1),1), VS(marker(:,1),2), VS(marker(:,1),3));
scatter3(VT(marker(:,2),1), VT(marker(:,2),2), VT(marker(:,2),3));
% add corresponding point
new_NS = zeros(length(VS),3);
new_NT = zeros(length(VT),3);
for i = 1 : length(FS)
    new_NS(FS(i,1),:) = NS(3*i-2,:);
    new_NS(FS(i,2),:) = NS(3*i-1,:);
    new_NS(FS(i,3),:) = NS(3*i,:);
end
for i = 1 : length(FT)
    new_NT(FT(i,1),:) = NT(3*i-2,:);
    new_NT(FT(i,2),:) = NT(3*i-1,:);
    new_NT(FT(i,3),:) = NT(3*i,:);
end

marker2 = [];
for i = 1 : length(VS)
    dist_set = sum((VT-repmat(VS(i,:), [length(VT),1])).*(VT-repmat(VS(i,:), [length(VT),1])),2);
    min_dist = min(dist_set);
    min_index = find(dist_set==min_dist);
    if (min_dist < dist && atan2d(norm(cross(new_NS(i,:),new_NT(min_index,:))),dot(new_NS(i,:),new_NT(min_index,:))) < theta)
        marker2 = [marker2;[i min_index]];
    end
end      

L = cotmatrix(VS, FS);
% C = sparse(length(marker),length(VS));
L = wl * L;

% normalize matrix
IJV = [(1:length(marker2))' marker2(:,1) ones(length(marker2),1)];
C = sparse(IJV(:,1), IJV(:,2), IJV(:,3), length(marker2), length(L));
C = wc*C;

% for i = 1:length(marker)
%     C(i,marker(i,1)) = 1;
% end
M = [L;C];

b = L*VS;
b = [b; VT(marker2(:,2),:)];
% for i = 1:length(marker)
%     b = [b;VT(marker(i,2),:)];
% end
x = M\b;

% normalize vertex
tmean = 0;
tstd = sqrt(2);
x = normPts(x, tmean, tstd);

figure, dispMesh(x, FS, [.8 .8 .8], 1);
hold on, dispMesh(VT, FT, [0 0 .8], .5);
end