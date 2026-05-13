# %%
import numpy as np
import matplotlib.pyplot as plt
import json 
#%%
G_val = 6.67430 * 10**(-8)
pi = 3.14159
to_cm = 1.49598 * 10**13    # 1 AU
to_g = 1.98854 * 10**33     # 1 Msun
to_cmpers = (G_val/to_cm*to_g)**0.5     # cm/s
to_kmpers = to_cmpers/10**5             # km/s
to_s = 1/(G_val/to_cm**3*to_g)**0.5       # s
c = 2.99792*10**10/to_cmpers


#%%

def abs_vector(x):
    return (x[0]**2+x[1]**2)**0.5

def Newtonian_accel(relM, R, V):  # 1기준. 2: m2, m1, -R, -V
    R_abs = abs_vector(R)
    rd = R@V/R_abs 
    an = -1*relM/R_abs**3 * R 
    an_d = -relM/R_abs**3 * V + relM/R_abs**4 *3*rd*R 
    
    return an, an_d

def PN_accel(myM, relM, R, V):
    
    R_abs = abs_vector(R) 
    aa = (myM+relM)/R_abs**3*(-R)
    V_abs = abs_vector(V)
    rd = R@V/R_abs
    M = myM + relM
    eta = myM*relM/M**2
    
    A = 8/5*eta*M/R_abs*rd*(17/3*M/R_abs + 3*V_abs**2)
    B = -8/5*eta*M/R_abs*(3*M/R_abs + V_abs**2)
    
    aPN = relM*(A*R/R_abs + B*V)/R_abs**2/c**5
    
    rdd = (V_abs**2 + R@aa - rd**2)/R_abs
    vd = V@aa/V_abs
    Ad = 8/5*eta*M*( 17/3*M*( rdd/R_abs**2 - 2*rd**2/R_abs**3 ) + 3*( rdd*V_abs**2/R_abs - rd**2*V_abs**2/R_abs**2 + rd/R_abs*2*V_abs*vd ) )
    Bd = 8/5*eta*M/R_abs**2 *rd *(3*M/R_abs + V_abs**2) - 8/5*eta*M/R_abs*(-3*M/R_abs**2*rd + 2*V@aa)
    
    aPN_d = relM* ( -3*A*rd*R/R_abs**4 + (Ad*R + A*V - 2*B*rd*V)/R_abs**3 + (Bd*V + B*aa)/R_abs**2 )/c**5

    return aPN, aPN_d

def update(m1, m2, dt, pos1, pos2, vel1, vel2):
    # dt 후의 위치, 속도 벡터 반환
    
    # 위치벡터, 상대속도
    R = pos1 - pos2
    V = vel1 - vel2
    a1, a1d = Newtonian_accel(m2, R, V)
    a2, a2d = Newtonian_accel(m1, -R, -V)
    
    # 임시 위치 속도 만들기
    rp1 = pos1 + vel1*dt + a1*dt**2/2 + a1d*dt**3/6 
    vp1 = vel1 + a1*dt + a1d*dt**2/2
    rp2 = pos2 + vel2*dt + a2*dt**2/2 + a2d*dt**3/6
    vp2 = vel2 + a2*dt + a2d*dt**2/2
    
    # predicted a 만들기 
    Rp = rp1 - rp2
    Vp = vp1 - vp2
    ap1, ap1d = Newtonian_accel(m2, Rp, Vp)
    ap2, ap2d = Newtonian_accel(m1, -Rp, -Vp)
    
    # 고차 미분 이용하여 다음 위치, 속도 추정
    a1ddd = 12*(a1 - ap1)/dt**3 + 6*(a1d + ap1d)/dt**2
    a1dd = -6*(a1 - ap1)/dt**2 - 2*(2*a1d + ap1d)/dt
    a2ddd = 12*(a2 - ap2)/dt**3 + 6*(a2d + ap2d)/dt**2
    a2dd = -6*(a2 - ap2)/dt**2 - 2*(2*a2d + ap2d)/dt
    
    pos1_next = rp1 + a1dd*dt**4/24 + a1ddd*dt**5/120
    pos2_next = rp2 + a2dd*dt**4/24 + a2ddd*dt**5/120
    vel1_next = vp1 + a1dd*dt**3/6 + a1ddd*dt**4/24
    vel2_next = vp2 + a2dd*dt**3/6 + a2ddd*dt**4/24
    
    return pos1_next, pos2_next, vel1_next, vel2_next

def updatePN(m1, m2, dt, pos1, pos2, vel1, vel2):
    # dt 후의 위치, 속도 벡터 반환
    
    # relative pos, vel, accel (1 기준, 2는 방향만 다 반대)
    R = pos1 - pos2
    V = vel1 - vel2
    # newtonian acceleration
    a1n, a1nd = Newtonian_accel(m2, R, V)
    a2n, a2nd = Newtonian_accel(m1, -R, -V)
    # PN aceeleration
    a1PN, a1PNd = PN_accel(m1, m2, R, V)
    a2PN, a2PNd = PN_accel(m2, m1, -R, -V)
    ###### PN 적용된 가속도 
    a1 = a1n + a1PN
    a1d = a1nd + a1PNd
    a2 = a2n + a2PN
    a2d = a2nd + a2PNd
    
    
    # 임시 위치 속도 만들기(동일)
    rp1 = pos1 + vel1*dt + a1*dt**2/2 + a1d*dt**3/6 
    vp1 = vel1 + a1*dt + a1d*dt**2/2
    rp2 = pos2 + vel2*dt + a2*dt**2/2 + a2d*dt**3/6
    vp2 = vel2 + a2*dt + a2d*dt**2/2
    
    
    ###### predicted a 만들기
    Rp = rp1 - rp2
    Vp = vp1 - vp2
    ap1n, ap1nd = Newtonian_accel(m2, Rp, Vp)
    ap2n, ap2nd = Newtonian_accel(m1, -Rp, -Vp)
    ap1PN, ap1PNd = PN_accel(m1, m2, Rp, Vp)
    ap2PN, ap2PNd = PN_accel(m2, m1, -Rp, -Vp)
    ap1 = ap1n + ap1PN
    ap1d = ap1nd + ap1PNd
    ap2 = ap2n + ap2PN 
    ap2d = ap2nd + ap2PNd


    # 고차 미분 이용하여 다음 위치, 속도 추정
    a1ddd = 12*(a1 - ap1)/dt**3 + 6*(a1d + ap1d)/dt**2
    a1dd = -6*(a1 - ap1)/dt**2 - 2*(2*a1d + ap1d)/dt
    a2ddd = 12*(a2 - ap2)/dt**3 + 6*(a2d + ap2d)/dt**2
    a2dd = -6*(a2 - ap2)/dt**2 - 2*(2*a2d + ap2d)/dt
    
    pos1_next = rp1 + a1dd*dt**4/24 + a1ddd*dt**5/120
    pos2_next = rp2 + a2dd*dt**4/24 + a2ddd*dt**5/120
    vel1_next = vp1 + a1dd*dt**3/6 + a1ddd*dt**4/24
    vel2_next = vp2 + a2dd*dt**3/6 + a2ddd*dt**4/24
    
    return pos1_next, pos2_next, vel1_next, vel2_next

def pred(m1, m2, e, a, Nperiod, timestepmax):
    period = (4*pi**2*a**3/(m1+m2))**0.5
    # m1, m2, e, a 로 초기 위치, 속도 벡터를 구하기
    # 질량중심계의 두 질량의 최대 거리 == m2 초점일 때 rmax = a(1+e) = (a1+a2)(1+e)와 같다.
    pos1 = np.array([-m2/(m1+m2)*a*(1+e), 0] )
    pos2 = np.array([m1/(m1+m2)*a*(1+e) ,0])
    plt.axis("equal")
    
    # 상대속도 이용하여 COM 좌표계에서의 속도 벡터
    v_a_rel = ((m1+m2)/a*(1-e)/(1+e))**0.5
    vel1 = np.array([0, m2/(m1+m2)*v_a_rel])  # 초기 up
    vel2 = np.array([0, -m1/(m1+m2)*v_a_rel])   # 초기 down

    # 주기 당 스텝 개수:
    period_step = int(timestepmax/Nperiod)
    
    # plot 위한 datapoint 저장
    x1_vals, y1_vals = np.zeros(timestepmax), np.zeros(timestepmax)
    x2_vals, y2_vals = np.zeros(timestepmax), np.zeros(timestepmax)
    
    # 스텝 개수만큼의 dt list를 빙빙돌기위한 dtlist: 원일점 시작 증가했다가 감소하면서 더해서 1주기되는 리스트 
    dt_list = [(i+1)**(2*e) for i in range(int(period_step/2))]
    dt_list = dt_list[::-1]+dt_list
    if period_step%2 != 0:
        dt_list.append(dt_list[-1])
    dt_list = [period/sum(dt_list)*dt_list[i] for i in range(len(dt_list))]
    
    # numerically 위치, 속도 예측
    step = 0
    time = 0
    P = 0 
    distance = ((pos2[0]-pos1[0])**2  + (pos2[1]-pos1[1])**2)**0.5
    print('STEP\tTIME(s)\t\t\tDISTANCE(AU)')
    for i in range(Nperiod): 
        P += 1
        for j in range(period_step):
            step += 1
            dt = dt_list[j]
            time += dt 
            # 위치, 속도 dt 후로 update
            pos1, pos2, vel1, vel2 = update(m1, m2, dt, pos1, pos2, vel1, vel2)
            distance = ((pos2[0]-pos1[0])**2  + (pos2[1]-pos1[1])**2)**0.5
            print('%d\t%1.11e\t%1.11e'%(step, time*to_s, distance))
            x1_vals[step-1], y1_vals[step-1] = pos1[0], pos1[1]
            x2_vals[step-1], y2_vals[step-1] = pos2[0], pos2[1]
    
    # answer = dict()
    # answer["period"] = P
    # answer["step"] = step
    # answer["time(s)"] = time*to_s
    # answer["distance(AU)"] = distance
    # with open('./output/output'+fnum, 'w', encoding='utf-8') as make_file:
    #     json.dump(answer, make_file, indent="\t")
    
    return x1_vals, y1_vals, x2_vals, y2_vals

def predPN(m1, m2, e, a, timestepmax):
    plt.figure(figsize = (20,20))
    # m1, m2, e, a 로 초기 위치, 속도 벡터를 구하기
    pos1 = np.array([-m2/(m1+m2)*a*(1+e), 0] )
    pos2 = np.array([m1/(m1+m2)*a*(1+e) ,0])
    plt.axis("equal")
    
    # 상대속도 이용하여 COM 좌표계에서의 속도 벡터 
    v_a_rel = ((m1+m2)/a*(1-e)/(1+e))**0.5
    vel1 = np.array([0, m2/(m1+m2)*v_a_rel])  # 초기 up
    vel2 = np.array([0, -m1/(m1+m2)*v_a_rel])   # 초기 down
    
    # plot 위한 datapoint 저장할 리스트 
    x1_vals, y1_vals =  [], []
    x2_vals, y2_vals = [], []

    # PN 위한 슈바르츠실트반지름
    Rsch = 2*(m1+m2)/c**2
    
    # numerically 위치, 속도 예측
    step = 0
    time = 0
    distance = ((pos2[0]-pos1[0])**2  + (pos2[1]-pos1[1])**2)**0.5
    print('STEP\tTIME(s)\t\t\tDISTANCE(AU)')
    while(True): 
        # dt = distance/100000 -> 40만 스텝에서 컴퓨터 멈춤
        # dt = distance/10000 -> 10만 스텝정도
        # dt = distance/50000 -> 멈춤
        # dt = distance/20000 -> 노트북에서 돌아가는 정도..
        dt = distance/75000 # olaf1 8코어에서 돌아감.
  
        step += 1
        time += dt
        if step>timestepmax:
            break
        elif distance < 2*Rsch: 
            # print('MERGE!!')
            # 거리가 슈바르츠실트보다 가까워지면 merge alert 날리고 while문 종료
            break
        
        pos1, pos2, vel1, vel2 = updatePN(m1, m2, dt, pos1, pos2, vel1, vel2)
        distance = ((pos2[0]-pos1[0])**2  + (pos2[1]-pos1[1])**2)**0.5
        
        print('%d\t%1.11e\t%1.11e'%(step, time*to_s, distance))
        
        x1_vals.append(pos1[0])
        y1_vals.append(pos1[1])
        x2_vals.append(pos2[0])
        y2_vals.append(pos2[1])
        
    return x1_vals, y1_vals, x2_vals, y2_vals
        
#%%
def TwoBody(fnum, name, m1, m2, e, a, Nperiod, timestepmax, PNon):
    if PNon == False:
        x1_vals, y1_vals, x2_vals, y2_vals = pred(m1, m2, e, a, Nperiod, timestepmax)
        plt.plot(x1_vals, y1_vals, color = 'grey', marker = 'o',markersize = 6*m1/(m1+m2), label='m1 = '+ str(m1)+' $M_\odot$')
        plt.plot(x2_vals, y2_vals, color = 'b', marker = 'o', markersize = 6*m2/(m1+m2), label='m2 = '+ str(m2)+' $M_\odot$')
    elif PNon == True:
        x1_vals, y1_vals, x2_vals, y2_vals = predPN(m1, m2, e, a, timestepmax)

        plt.plot(x1_vals, y1_vals, color = 'red', linestyle = 'solid', linewidth = m2/(m1+m2)*4,label='m1 = '+ str(m1)+' $M_\odot$')
        plt.plot(x2_vals, y2_vals, color = 'blue', linestyle = 'solid', linewidth = m1/(m1+m2)*4,label='m2 = '+ str(m2)+' $M_\odot$')
        scale = 10**-6
        if fnum == '3.2':
            scale = scale/2
        plt.xlim([-scale, scale])
        plt.ylim([-scale, scale])
    plt.title('Orbit Simulation: '+name)
    plt.xlabel('X position (AU)')
    plt.ylabel('Y position (AU)')
    plt.legend()
    plt.savefig('./FIG/'+name+'.png')
    plt.show()

#%%
# 1.2
m1 = 1  # solar mass
m2 = 1  # solar mass
a = 2   # AU
e = 0
Nperiod = 100 # 100주기
timestepmax = 20000
PNon = False # PN x 

#TwoBody('1.2','Q1-2',m1,m2,e,a,Nperiod,timestepmax,PNon)

#%%
# 1.3
m1 = 1  # solar mass
m2 = 1  # solar mass
a = 2   # AU
e = 0
Nperiod = 100 # 100주기
timestepmax = 40000
PNon = False # PN x 

#TwoBody('1.3','Q1-3',m1,m2,e,a,Nperiod,timestepmax,PNon)


#%%
# 1.4
m1 = 10  # solar mass
m2 = 2  # solar mass
a = 100   # AU
e = 0
Nperiod = 100 # 100주기
timestepmax = 20000
PNon = False # PN x 

#TwoBody('1.4','Q1-4',m1,m2,e,a,Nperiod,timestepmax,PNon)

# %%
# 2.1
m1 = 1  # solar mass
m2 = 1  # solar mass
a = 2   # AU
e = 0.5
Nperiod = 100 # 100주기
timestepmax = 40000
PNon = False # PN x 

#TwoBody('2.1','Q2-1',m1,m2,e,a,Nperiod,timestepmax,PNon)

# %%
# 2.2
m1 = 1  # solar mass
m2 = 2.22*10**14 /to_g  # solar mass
a = 17.8   # AU
e = 0.967
Nperiod = 50
timestepmax = 40000
PNon = False # PN x 

#TwoBody('2.2','Q2-2',m1,m2,e,a,Nperiod,timestepmax,PNon)

# %%
# 3.1
m1 = 10
m2 = 10
a = 2*10**(-5)
e = 0 
Nperiod = None
timestepmax = 1000000
PNon = True
#TwoBody('3.1','Q3-1',m1,m2,e,a,Nperiod,timestepmax,PNon)

# %%
# 3.2
m1 = 10
m2 = 2
a = 2*10**(-5)
e = 0.8
Nperiod = None
timestepmax = 1000000
PNon = True
#TwoBody('3.2','Q3-2',m1,m2,e,a,Nperiod,timestepmax,PNon)
# %%
