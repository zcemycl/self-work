clear all; clc
% K-means
% 3 clusters sample

% cluster 1: mean = [2,3] sigma = [1 1.5; 1.5 3]
m1 = [2,3]; s1 = [1 1.5; 1.5 3];
R1 = mvnrnd(m1,s1,500);
% cluster 2: mean = [7,7] sigma = [1 0; 0 1]
m2 = [7,7]; s2 = [1 0; 0 1];
R2 = mvnrnd(m2,s2,250);
% cluster 2: mean = [5,2] sigma = [1 0.5; 1.5 1]
m3 = [5,2]; s3 = [2 1; 1 2];
R3 = mvnrnd(m3,s3,700);

subplot(2,2,1)
hold on 
plot(R1(:,1),R1(:,2),'g+');
plot(R2(:,1),R2(:,2),'ko');
plot(R3(:,1),R3(:,2),'rx');
hold off

% whole dataset = R1,R2,R3
R = [R1;R2;R3];
% figure(1);
% plot(R(:,1),R(:,2),'b.');
Ra = R(:,1); Rb = R(:,2);

% draw initial K points for grouping
G = datasample(R,3);
hold on 
plot(G(:,1),G(:,2),'k*');

Cost = [];

for iter = 1:20
ld = zeros(length(R),2); %label and distance
for i = 1:length(R)
    d1 = norm(R(i,:)-G(1,:));
    d2 = norm(R(i,:)-G(2,:));
    d3 = norm(R(i,:)-G(3,:));
    D = [d1,d2,d3];
    
    %minimum distance
    mindis = min(D);
    %label
    label = find(D == mindis);
    
    ld(i,:) = [mindis,label];  
end
i1 = find(ld(:,2)==1);
tmp1 = [Ra(i1),Rb(i1)];
i2 = find(ld(:,2)==2);
tmp2 = [Ra(i2),Rb(i2)];
i3 = find(ld(:,2)==3);
tmp3 = [Ra(i3),Rb(i3)];
G(1,:) = [sum(tmp1(:,1)),sum(tmp1(:,2))]/length(tmp1);
G(2,:) = [sum(tmp2(:,1)),sum(tmp2(:,2))]/length(tmp2);
G(3,:) = [sum(tmp3(:,1)),sum(tmp3(:,2))]/length(tmp3);


subplot(2,2,3)
hold on 
plot(tmp1(:,1),tmp1(:,2),'g+')
plot(tmp2(:,1),tmp2(:,2),'k+')
plot(tmp3(:,1),tmp3(:,2),'r+')
plot(G(:,1),G(:,2),'bo')
drawnow;
%hold off

Cost(iter) = sum(ld(:,1)); %
subplot(2,2,4)
%hold on
plot(Cost)
drawnow;
%hold off

end
