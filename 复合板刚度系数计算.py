# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 18:53:54 2023

@author: Att1ce
"""

import numpy as np

num = int(input('请输入铺层数:'))
t = [0]*num#铺层角度
for i in range(len(t)):
    t[i] = float(input('请依次输入铺层角度:'))
h = float(input('请输入铺层厚度：'))#铺层厚度


print('铺层角度：',t,'每个铺层的厚度：',h)

print('单位为Mpa,N/mm2')

E1 = float(input('E1:'))
E2 = float(input('E2:'))
nu12 = float(input('nu12:'))
G12 = float(input('G12:'))


def calc_stiffnes_matrix(E1,E2,nu12,G12):
   
    #定义各种弹性常数
    nu21 = nu12*E2/E1
    Q_11 = E1/(1-nu12*nu21)
    Q_22 = E2/(1-nu12*nu21)
    Q_12 = nu12*E2/(1-nu12*nu21)
    Q_66 = G12
    print('开始计算刚度矩阵')
    
    Q = np.array([[Q_11,Q_12,0],
               [Q_12,Q_22,0],
               [0,0,Q_66]])
    print('Q:',np.around(Q,decimals=2))
    
    print('开始计算ABD矩阵')  
    #ABD矩阵计算，注意，这是左手系
    z0 = h*num/2
    A = np.zeros((3,3))  
    B = np.zeros((3,3))     
    D = np.zeros((3,3))  
    for i in range(len(t)):
        angle=t[i]*np.pi/180
        #计算坐标变换矩阵
        T_ni = np.array([[(np.cos(angle))**2,(np.sin(angle))**2,-2*np.cos(angle)*np.sin(angle)],
                    [(np.sin(angle))**2,(np.cos(angle))**2,-2*np.cos(angle)*np.sin(angle)],
                    [np.cos(angle)*np.sin(angle),-np.cos(angle)*np.sin(angle),(np.cos(angle)**2-(np.sin(angle))**2)]])
        T_Q = np.dot(T_ni,Q)
        Q_2 = np.dot(T_Q, T_ni.T)#此为坐标变换关系Q_
        
        zi = i * h + z0
        z_i = 0.5 * (zi+zi-h)
        A += Q_2 * h
        B += Q_2 * z_i * h
        D += Q_2 * (z_i**2 + 1/12 * h **2) * h
        
    print('A:',np.around(A,decimals=1))
    print('B:',np.around(B,decimals=1))
    print('D:',np.around(D,decimals=1))
    
calc_stiffnes_matrix(E1,E2,nu12,G12)

           
        
        
        
        
        