# package
import matplotlib.pyplot as plt
import numpy as np
import os
import natsort
import pycbc
from pycbc import psd
from pycbc import frame
from pycbc.waveform import get_td_waveform, get_fd_waveform
from pycbc.detector import add_detector_on_earth, Detector

# cbc information
# list_etmdc1.xlsx 에서 정보 가져오기
# 예시) index 44413 (SNR=43.428922) BBH
t0, tc = 1001688810.88445, 1001688818.2177 
m1, m2 = 123.762727, 95.373487 
s1, s1x, s1y, s1z, s2, s2x, s2y, s2z = 0.2432, 0.119985, 0.027662, -0.209725, 0.1836, -0.180167, 0.008, -0.034423
lambda1, lambda2 = 0, 0
z, dL, ra, decl, psi, inc, phi0 = 1.9425, 15367.8861, 0.631436, -0.356858, 0.430408, 1.857159, 0.312234


# gwf 파일 몇번째에 들어있는 데이터인지 체크하기
fnum = int((t0-1000000000)/2048)
print(fnum)
print(1000000000+2048*fnum)
print(1000000000+2048*(fnum+1))


MDC_PATH="/cvmfs/.../data/" 
E1_file_list=os.listdir(MDC_PATH+"/E1")
E1_file_list=natsort.natsorted(E1_file_list)

E2_file_list=os.listdir(MDC_PATH+"/E2")
E2_file_list=natsort.natsorted(E2_file_list)

E3_file_list=os.listdir(MDC_PATH+"/E3")
E3_file_list=natsort.natsorted(E3_file_list)


# 전체 데이터 읽기
print(MDC_PATH+'/E1/'+E1_file_list[fnum])
data1_read=frame.read_frame(MDC_PATH+'/E1/'+E1_file_list[fnum],'E1:STRAIN')
data2_read=frame.read_frame(MDC_PATH+'/E2/'+E2_file_list[fnum],'E2:STRAIN') 
data3_read=frame.read_frame(MDC_PATH+'/E3/'+E3_file_list[fnum],'E3:STRAIN') 


# data가 길어서 일부만 잘라서 보고싶은 경우 crop 해주기
data1_read = data1_read.crop(t0-1000000000-2048*fnum-12, 1000000000+2048*(fnum+1)-tc-12)
data2_read = data2_read.crop(t0-1000000000-2048*fnum-12, 1000000000+2048*(fnum+1)-tc-12)
data3_read = data3_read.crop(t0-1000000000-2048*fnum-12, 1000000000+2048*(fnum+1)-tc-12)


# time series data --> frequency series data
data1_tilde = data1_read.to_frequencyseries()
data2_tilde = data2_read.to_frequencyseries()
data3_tilde = data3_read.to_frequencyseries()
print(len(data1_tilde))


# PSD(frequency series) 불러와서 데이터랑 길이 맞추기
filename = './PSD_.txt'

delta_f = data1_read.delta_f
length = int(4096 / delta_f) + 1
low_frequency_cutoff = 5
PSD = psd.from_txt(filename, length, delta_f,
                         low_frequency_cutoff, is_asd_file=False)
print(len(PSD))


# template 형성, 길이 맞추기
## argument 이름 정확히 적어주기
hp, hc = pycbc.waveform.get_fd_waveform(approximant='IMRPhenomXPHM',
                             mass1=m1, mass2=m2,
                             delta_f = delta_f, 
                             f_lower=5,
                             spin1x = s1x, spin1y = s1y, spin1z = s1z, 
                             spin2x = s2x, spin2y = s2y, spin2z = s2z,
                             lambda1 = lambda1, lambda2 = lambda2, 
                             redshift=z,distance=dL,ra=ra, 
                             dec=decl,polarization=psi,inclination=inc,coa_phase=phi0)

hp.resize(len(data1_tilde))
hc.resize(len(data1_tilde))



# antenna pattern 
# antenna pattern
E1_lon = None  # [deg] * np.pi/180 → [rad], 비공개
E1_lat = None  # [deg] * np.pi/180 → [rad], 비공개
E1_xangle = None  # 비공개
E1_yangle = None  # 비공개
E1_elevation = None  # 비공개

E2_lon = None  # [deg] * np.pi/180 → [rad], 비공개
E2_lat = None  # [deg] * np.pi/180 → [rad], 비공개
E2_xangle = None  # 비공개
E2_yangle = None  # 비공개
E2_elevation = None  # 비공개

E3_lon = None  # [deg] * np.pi/180 → [rad], 비공개
E3_lat = None  # [deg] * np.pi/180 → [rad], 비공개
E3_xangle = None  # 비공개
E3_yangle = None  # 비공개
E3_elevation = None  # 비공개

# Detector 추가 
add_detector_on_earth("E1",E1_lon,E1_lat,yangle=E1_yangle,xangle=E1_xangle,height=E1_elevation,xlength=10000,ylength=10000)
add_detector_on_earth("E2",E2_lon,E2_lat,yangle=E2_yangle,xangle=E2_xangle,height=E2_elevation,xlength=10000,ylength=10000)
add_detector_on_earth("E3",E3_lon,E3_lat,yangle=E3_yangle,xangle=E3_xangle,height=E3_elevation,xlength=10000,ylength=10000)

d1,d2,d3=Detector("E1"),Detector("E2"),Detector("E3")



# antenna pattern 적용 된 real template 만들기

fp1, fc1=d1.antenna_pattern(ra,decl,psi, tc)
fp2, fc2=d2.antenna_pattern(ra,decl,psi, tc)
fp3, fc3=d3.antenna_pattern(ra,decl,psi, tc)
print("fp1={},fc1={}".format(fp1,fc1))
print("fp2={},fc2={}".format(fp2,fc2))
print("fp3={},fc3={}".format(fp3,fc3))

# real template!
ht1=fp1*hp+fc1*hc
ht2=fp2*hp+fc2*hc
ht3=fp3*hp+fc3*hc


# real template과 data의 snr 계산
snr1 = pycbc.filter.matched_filter(ht1, data1_tilde, psd=PSD, low_frequency_cutoff=5)
snr2 = pycbc.filter.matched_filter(ht2, data2_tilde, psd=PSD, low_frequency_cutoff=5)
snr3 = pycbc.filter.matched_filter(ht3, data3_tilde, psd=PSD, low_frequency_cutoff=5)


# network snr 계산
snr_net = (abs(snr1)**2 + abs(snr2)**2 + abs(snr3)**2)**0.5
peak = snr_net.numpy().argmax()
snrp = snr_net[peak]
time = snr_net.sample_times[peak]
print("SNR_net: {}s, {}".format(time, snrp))
print('t0 = ', t0, ', tc = ', tc)
plt.plot(snr_net.sample_times, snr_net)
plt.ylabel('network snr')
plt.xlabel('time (s)')
plt.show()
