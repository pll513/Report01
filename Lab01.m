syms x;     %y=f(x)�е��Ա���x����ʾ�ɻ���ˮƽ����
syms A;     %f(x)��������ϵ��
syms B;     %f(x)�Ķ�����ϵ��
syms C;     %f(x)��һ����ϵ��
syms D;     %f(x)�ĳ�����
syms l;     %�ɻ��ཱུ����ˮƽ����
syms h;     %�ɻ��ĳ�ʼ����߶�
syms u;     %�ɻ���ˮƽ�����ٶ�
syms g;     %�������ٶ�,ʵ����ȡ9.8m/(s*s)
syms T;     %�ɻ��Ļ���ʱ��
syms a;     %�ɻ����й����еļ��ٶ�,ȡ��ֵ,-a��ʾ����
syms d;     %�ɻ��ܵ��ĳ���,���ɻ����о���
syms t;     %�ɻ�����ʱ��
format long;


%--------------------  ʵ��һ(1)  --------------------%
%--------------  ���ɻ��Ľ�������f(x)  --------------%
R = solve(A*(-l)^3+B*(-l)^2+C*(-l)+D==h,D==0,...
    3*A*(-l)^2+2*B*(-l)+C==0,C==0,A,B,C,D);
r1 = double(subs(R.A,{l,h},{80,10}));
r2 = double(subs(R.B,{l,h},{80,10}));
r3 = double(subs(R.C,{l,h},{80,10}));
r4 = double(subs(R.D,{l,h},{80,10}));
F = [r1 r2 r3 r4]


%--------------------  ʵ��һ(2)  --------------------%
%------------  ��ˮƽ����l���������Сֵ  -------------%
l_min = solve(subs(6*h*u^2/l^2+12*h*u^2*x/l^3==g/10,x,0),l);
l_min = vpa(subs(l_min(1),{h,u,g},{10,540/3600,0.0098}))


%--------------------  ʵ��һ(3)  --------------------%
%----------------  ģ��ɻ�����ȫ����  ----------------%
%һЩ׼������
Q = solve(u*T-a*T^2/2==d,u-a*T==0,a,T);
q1 = double(subs(Q.a,{d,u},{3.6,540}));
q2 = double(subs(Q.T,{d,u},{3.6,540}));
G = 0;
T1 = double(80/540);
T2 = T1+q2;

%�ɻ����������x��y����t�ĺ���
F_x_t = [540 -80];
F_y_t = sym2poly(subs(3*h*x^2/l^2+2*h*x^3/l^3,...
    {h,x,l},{10,540*t-80,80}));

%�ɻ����й�����x��y����t�ĺ���(�ٶ�t��0��ʼȡ,ʵ�����Ǵ�T1��ʼȡ)
G_x_t = [ -540^2/(4*3.6) 540 0];
G_y_t = 0;

%��ʼģ��
h=figure('numbertitle','off','name',' Lab01');
%���Ʒɻ��Ľ�������
xx = linspace(-80,0,5000);
yy = polyval(F,xx);
g1 = plot(xx,yy);
hold on;
%���Ʒɻ��Ļ�������
xx = linspace(0,3.6,300);
yy = polyval(G,xx);
g2 = plot(xx,yy);
hold on;
%һЩ��ʾ�ĵ���
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
