syms x;     %y=f(x)中的自变量x，表示飞机的水平坐标
syms A;     %f(x)的三次项系数
syms B;     %f(x)的二次项系数
syms C;     %f(x)的一次项系数
syms D;     %f(x)的常数项
syms l;     %飞机距降落点的水平距离
syms h;     %飞机的初始降落高度
syms u;     %飞机的水平飞行速度
syms g;     %重力加速度,实验中取9.8m/(s*s)
syms T;     %飞机的滑行时间
syms a;     %飞机滑行过程中的加速度,取正值,-a表示减速
syms d;     %飞机跑道的长度,即飞机滑行距离
syms t;     %飞机飞行时间
format long;


%--------------------  实验一(1)  --------------------%
%--------------  求解飞机的降落曲线f(x)  --------------%
R = solve(A*(-l)^3+B*(-l)^2+C*(-l)+D==h,D==0,...
    3*A*(-l)^2+2*B*(-l)+C==0,C==0,A,B,C,D);
r1 = double(subs(R.A,{l,h},{80,10}));
r2 = double(subs(R.B,{l,h},{80,10}));
r3 = double(subs(R.C,{l,h},{80,10}));
r4 = double(subs(R.D,{l,h},{80,10}));
F = [r1 r2 r3 r4]


%--------------------  实验一(2)  --------------------%
%------------  求水平距离l所允许的最小值  -------------%
l_min = solve(subs(6*h*u^2/l^2+12*h*u^2*x/l^3==g/10,x,0),l);
l_min = vpa(subs(l_min(1),{h,u,g},{10,540/3600,0.0098}))


%--------------------  实验一(3)  --------------------%
%----------------  模拟飞机降落全过程  ----------------%
%一些准备工作
Q = solve(u*T-a*T^2/2==d,u-a*T==0,a,T);
q1 = double(subs(Q.a,{d,u},{3.6,540}));
q2 = double(subs(Q.T,{d,u},{3.6,540}));
G = 0;
T1 = double(80/540);
T2 = T1+q2;

%飞机降落过程中x和y关于t的函数
F_x_t = [540 -80];
F_y_t = sym2poly(subs(3*h*x^2/l^2+2*h*x^3/l^3,...
    {h,x,l},{10,540*t-80,80}));

%飞机滑行过程中x和y关于t的函数(假定t从0开始取,实际上是从T1开始取)
G_x_t = [ -540^2/(4*3.6) 540 0];
G_y_t = 0;

%开始模拟
h=figure('numbertitle','off','name',' Lab01');
%绘制飞机的降落曲线
xx = linspace(-80,0,5000);
yy = polyval(F,xx);
g1 = plot(xx,yy);
hold on;
%绘制飞机的滑行曲线
xx = linspace(0,3.6,300);
yy = polyval(G,xx);
g2 = plot(xx,yy);
hold on;
%一些显示的调整
set(g1,'Linestyle','-','color','g','Linewidth',1.5);
set(g2,'Linestyle','-','color','m','Linewidth',1.5);
xlabel('x');
ylabel('y');
axis([-90,10,0,10.05]); 

x0 = -80;
y0 = 10;
p1 = plot(x0,y0,'color','r','marker','.','markersize',15); 
t = 0;
dt = double(80/540000);
while 1
    if ~ishandle(h),return,end
    if (t<T1)
        set(p1,'xdata',polyval(F_x_t,t),'ydata',polyval(F_y_t,t));
    elseif ( t>=T1 && t<=T2 )
        set(p1,'xdata',polyval(G_x_t,(t-T1)),'ydata',...
            polyval(G_y_t,(t-T1)));
    else
        t = 0;
    end
    t = t + dt;
    drawnow
end
