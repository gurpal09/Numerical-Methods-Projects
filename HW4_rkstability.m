theta = 0:0.001:2*pi;
a = zeros(4,length(theta));

for k = 1:length(theta)
    c = [1./24 1./6 0.5 1 1-exp(i*theta(k))];
    a(:,k) = roots(c);
end

hold on
plot(a(1,:),'co:')
plot(a(2,:),'co:')
plot(a(3,:),'co:')
plot(a(4,:),'co:')
hold off