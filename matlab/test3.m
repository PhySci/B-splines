n = 10;
x = linspace(0,5,10);
y = x.^3+1;
w = ones(n,1);

spaps1(x,y,-1,w,2)
plot(x,y,'ro')
xlim([-0.1, 5.1])
ylim([-1 130])


