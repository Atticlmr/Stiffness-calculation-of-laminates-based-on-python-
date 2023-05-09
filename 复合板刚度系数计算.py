# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 18:53:54 2023

@author: 20051009 李明睿
"""

import numpy as np
import pandas as pd

h = pd.read_csv("data1.csv",usecols=["	铺层厚度/mm"],encoding='gbk')
h = int(np.array(h)[0][0])

num = pd.read_csv("data1.csv",usecols=["铺层层数"],encoding='gbk')
num = int(np.array(num)[0][0])


E1 = pd.read_csv("data1.csv",usecols=["	E_1/Mpa"],encoding='gbk')
E1 = float(np.array(E1)[0][0])

E2 = pd.read_csv("data1.csv",usecols=["	E_2/Mpa"],encoding='gbk')
E2 = float(np.array(E2)[0][0])

nu12 = pd.read_csv("data1.csv",usecols=["	v_12"],encoding='gbk')
nu12 = float(np.array(nu12)[0][0])

G12 = pd.read_csv("data1.csv",usecols=["	G_12"],encoding='gbk')
G12 = float(np.array(G12)[0][0])

#print(num)
t = [0]*num#铺层角度
angle = pd.read_csv("data1.csv",usecols=["	铺层角度（角度制）"],encoding='gbk')

for i in range(len(t)):
    t[i] = float(np.array(angle)[i][0])




print('铺层角度：',t,'每个铺层的厚度：',h)

print('单位为Mpa,N/(mm^2)')




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

           
        
        
        
        
        