import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib import animation
from scipy.integrate import solve_ivp


G, M, R = 6.67384e-11, 5.972e+24, 6.3781e+6

#지구
plt.clf()
fig = plt.figure(1, figsize = (10,10))
ax = plt.axes(xlim=(-7e+6, 7e+6), ylim=(-7e+6, 7e+6))
title = ax.set_title('')
ax.plot([0], [0], 'k.', label='center of mass of Earth')
Ex = np.linspace(-R, R, 500)
Ey = []
Ey.append((R*R - Ex*Ex)**(0.5))
Ey.append(-(R*R - Ex*Ex)**(0.5))
j = 0

while j<2:
    ax.plot(Ex, Ey[j], 'b-', label='Earth')
    j = j + 1

#우주정거장 궤도 (근지점 408km, 원지점 410km이므로 409km의 원궤도로 함)
MISS , RISS, vISS = 419725, 409e+3, 7706.6
ISS, = ax.plot([], [], 'r.', label='ISS')
dt = 0.01

def updateISS(a):
    t = dt * a
    x, y = (R + RISS)*np.cos((vISS/(R+RISS)) * t), (R + RISS)*np.sin((vISS/(R+RISS)) * t)
    ISS.set_data([x], [y])
    
animISS = animation.FuncAnimation(fig, updateISS, frames=10000, interval=dt*1000)

xISS = []
yISS = []
e = 0

while e < 10000:         #뒤에 숫자는 frames 값과 같음
    tISS = dt * e
    xISS.append((R + RISS)*np.cos((vISS/(R+RISS)) * tISS))
    yISS.append((R + RISS)*np.sin((vISS/(R+RISS)) * tISS))
    e = e + 1

#우주선 쏘기
m = 50000
question = input('포물선 운동(1), 로켓운동(2) --> ')
question = int(question)

#초기속도에 의한 포물선 운동을 하는 경우
if question ==  1:

    def f1(t, u, G, m, M):
        x, y, vx, vy = u
        A = -G*M/((x)**2+(y)**2)**1.5
        return vx, vy, A*x, A*y
        
    
    v0min = input('구하고자 하는 속도 범위의 최솟값을 넣으세요 : ')
    v0max = input('구하고자 하는 속도 범위의 최댓값을 넣으세요 : ')
    v0min = int(v0min)
    v0max = int(v0max)
    
    def fmin1(v0min, v0max):
        t_max = 100
        n = 0
        num = 0
        rmin = []
        v0 = v0min
        
        while v0 < v0max:
            u0 = (R, 0, v0*10, v0)              #v0에 어떠한 배수 곱해서 x축 속도와 y축 속도 크기 조정
            sol = solve_ivp(f1,(0,t_max), u0, args=(G,m,M), t_eval=np.arange(0, t_max, dt), rtol=1e-9)      
            while n < 100:            #요거는 frames/100와 같음
                r = []
                r.append((sol.y[0][n*100]-xISS[n*100])**2 + (sol.y[1][n*100]-yISS[n*100])**2)
                n = n + 1
            rmin.append(min(r))
            num = num + 1
            v0 = v0 + 1
            
        rminmin = min(rmin)
        i = 0
        
        while i < num:
            if rmin[i] == rminmin:
                rv0 = v0min + i
            i = i + 1
        
        print(rv0, '일 때 우주정거장과 가장 근접함')
        
        ru0 = (R, 0, rv0*10, rv0)           #rv0에 어떠한 배수 곱해서 x축 속도와 y축 속도 크기 조정 (위와 같은 값)
        return solve_ivp(f1,(0,t_max), ru0, args=(G,m,M), t_eval=np.arange(0, t_max, dt), rtol=1e-9)
    
    rsol = fmin1(v0min, v0max)

#로켓 운동 하는 경우
elif question == 2:
    vrel = 40    #분출속도 (나로호 궤도 참조)
    mdt = 5    #초당 손실 질량 
    wmin = input('구하고자 하는 각도의 최솟값을 넣으세요(단위: 도) : ')       #w는 지표면과 가속도가 이루는 각도   
    wmin = int(wmin)
    wmax = input('구하고자 하는 각도의 최댓값을 넣으세요(단위: 도) : ')
    wmax = int(wmax)
    
    def f2(t, u, G, m, M, w):
        x, y, vx, vy = u
        A = -G*M/((x)**2+(y)**2)**1.5
        B = vrel*(mdt/m)*np.sin(w*(math.pi/180))
        C = vrel*(mdt/m)*np.cos(w*(math.pi/180))
        D = 1/((x)**2+(y)**2)**0.5
        return vx, vy, B*D*x-C*D*y+A*x, B*D*y+C*D*x+A*y
    
    def fmin2(wmin, wmax):
        t_max = 100
        n = 0
        num = 0
        rmin = []
        w = wmin
        
        while w < wmax:
            u0 = (R, 0, 0, 0)
            sol = solve_ivp(f2,(0,t_max), u0, args=(G,m,M,w), t_eval=np.arange(0, t_max, dt), rtol=1e-9)      
            while n < 100:            #요거는 frames/100와 같음
                r = []
                r.append((sol.y[0][n*100]-xISS[n*100])**2 + (sol.y[1][n*100]-yISS[n*100])**2)
                n = n + 1
            rmin.append(min(r))
            num = num + 1
            w = w + 1
            
        rminmin = min(rmin)
        i = 0
        
        while i < num:
            if rmin[i] == rminmin:
                rw = wmin + i
            i = i + 1
        w = rw
        print(w, '(도)일 때 우주정거장과 가장 근접함')
        
        return solve_ivp(f2,(0,t_max), u0, args=(G,m,M,w), t_eval=np.arange(0, t_max, dt), rtol=1e-9)

    rsol = fmin2(wmin,wmax)
    
spaceship, = ax.plot([],[], 'g.', label = 'spaceship')
spaceship_trace, = ax.plot([], [], 'k:')
x_list, y_list = [], []

def updates(i):
    x, y = rsol.y[0][i],rsol.y[1][i]
    x_list.append(x); y_list.append(y)
    spaceship.set_data([x], [y])
    spaceship_trace.set_data(x_list, y_list)
    title.set_text(f'Time = {rsol.t[i]:7.2f}')
anims = animation.FuncAnimation(fig, updates, frames=len(rsol.t), interval=dt*1000)
plt.show()








