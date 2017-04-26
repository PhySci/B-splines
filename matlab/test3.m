n = 100;
x = linspace(0,10,n);
y = cos(x) %x.^3+1;
w = ones(n,1);

sp1 = spaps1(x,y,-1,w,2)
sp2 = spaps1(x,y,-0.5,w,2)
sp3 = spaps1(x,y,-0.1,w,2)

subplot(311)
plot(x,y,'ro',x,fnval(x,sp1))
legend('p = 1')

subplot(312)
plot(x,y,'ro',x,fnval(x,sp2))
legend('p = 0.5')

subplot(313)
plot(x,y,'ro',x,fnval(x,sp3))
legend('p = 0.1')


