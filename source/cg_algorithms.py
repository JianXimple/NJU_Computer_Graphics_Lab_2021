#!/usr/bin/env python
# -*- coding:utf-8 -*-

# 本文件只允许依赖math库
import math
from typing import final

from numpy import place, square


def draw_line(p_list, algorithm):
    """绘制线段

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'，此处的'Naive'仅作为示例，测试时不会出现
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    result = []
    if algorithm == 'Naive':
        if x0 == x1:
            for y in range(y0, y1 + 1):
                result.append((x0, y))
        else:
            if x0 > x1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            k = (y1 - y0) / (x1 - x0)
            for x in range(x0, x1 + 1):
                result.append((x, int(y0 + k * (x - x0))))
    elif algorithm == 'DDA':
        if x0==x1:
            for y in range(min(y0,y1),max(y0,y1)+1):
                result.append((x0,y))
        else:
            k=k = (y1 - y0) / (x1 - x0)
            if abs(k)>1:
                if y0>y1:
                    x0,y0,x1,y1=x1,y1,x0,y0
                x=x0
                for y in range(y0,y1+1):
                    result.append((int(x),y))
                    x+=1/k
            else:
                if x0>x1:
                    x0,y0,x1,y1=x1,y1,x0,y0
                y=y0
                for x in range(x0,x1+1):
                    result.append((x,int(y)))
                    y+=k
    elif algorithm == 'Bresenham':
        if x0==x1:
            if y0>y1:
                y0,y1=y1,y0
            for y in range(y0, y1 + 1):
                result.append((x0, y))
        else:
            if x0 > x1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            k = (y1 - y0) / (x1 - x0)
            def cal(x0,y0,x1,y1,delta_x,delta_y,d):
                p0=2*delta_y-delta_x
                for x in range(x0,x1+1):
                    if p0<0:
                        result.append((x,y0))
                        p0+=2*delta_y
                    else:
                        y0+=d
                        result.append((x,y0))
                        p0=p0+2*delta_y-2*delta_x
                return 0
            if abs(k)<1:
                if k>0:
                    delta_x=x1-x0
                    delta_y=y1-y0
                    cal(x0,y0,x1,y1,delta_x,delta_y,1)
                else:
                    delta_x = x1 - x0
                    delta_y = y1 - y0
                    delta_y = -delta_y
                    cal(x0,y0,x1,y1,delta_x,delta_y,-1)
            else:
                def caly(x0,y0,x1,y1,delta_x,delta_y,d):
                    p0=2*delta_x-delta_y
                    for y in range(y0,y1+1):
                        if p0<0:
                            result.append((x0,y))
                            p0+=2*delta_x
                        else:
                            x0 += d
                            result.append((x0,y))
                            p0=p0+2*delta_x-2*delta_y
                    return 0
                if y0 > y1:
                    x0, y0, x1, y1 = x1, y1, x0, y0
                if k>0:
                    delta_x=x1-x0
                    delta_y=y1-y0
                    caly(x0,y0,x1,y1,delta_x,delta_y,1)
                else:
                    delta_y = y1 - y0
                    delta_x = x1 - x0
                    delta_x = -delta_x
                    caly(x0,y0,x1,y1,delta_x,delta_y,-1)
    return result


def draw_polygon(p_list, algorithm):
    """绘制多边形

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 多边形的顶点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    #print("start polygon")
    result = []
    for i in range(len(p_list)-1):
        line = draw_line([p_list[i], p_list[i+1]], algorithm)
        result += line
    line=draw_line([p_list[0], p_list[-1]], algorithm)
    result +=line
    return result


def draw_ellipse(p_list,algorithm):
    """绘制椭圆（采用中点圆生成算法）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 椭圆的矩形包围框左上角和右下角顶点坐标
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """

    result=[]
    x0,y0=p_list[0]
    x1,y1=p_list[1]
    cx=int((x0+x1)/2)
    cy=int((y0+y1)/2)
    rx=int(abs(x0-x1)/2)
    ry=int(abs(y0-y1)/2)
    #result.append((0,ry))
    x=0
    y=ry
    result.append((x,y))
    result.append((-x,y))
    result.append((x,-y))
    result.append((-x,-y))
    #p1_0=float(pow(ry,2)*pow((x+1),2)+pow(rx,2)*pow((y-1/2),2)-pow(rx,2)*pow(ry,2))
    p1_0=float(pow(ry,2)-pow(rx,2)*ry+pow(rx,2)/4)
    while(pow(rx,2)*y>pow(ry,2)*x):
        result.append((x,y))
        if p1_0<0:
            p1_0+=float(2*pow(ry,2)*x+3*pow(y,2))
            x=x+1
            y=y
        else:
            p1_0+=float(2*pow(ry,2)*x-2*pow(rx,2)*y+2*pow(rx,2)+3*pow(ry,2))
            x=x+1
            y=y-1
    (x,y)=result[-1]
    #p2_0=float(pow(ry,2)*pow((x+1/2),2)+pow(rx,2)*pow((y-1),2)-pow(rx,2)*pow(ry,2))
    p2_0=float(pow(ry,2)*pow(x+0.5,2)+pow(rx,2)*pow(y-1,2)-pow(rx,2)*pow(ry,2))
    while y>=0:
        result.append((x,y))
        if p2_0>0:
            p2_0+=float((-2)*pow(rx,2)*y+3*pow(rx,2))
            x=x
            y=y-1
        else:
            p2_0+=float((-2)*pow(rx,2)*y+3*pow(rx,2)+2*pow(ry,2)*x+2*pow(ry,2))
            x=x+1
            y=y-1
    final=[]
    temp=[]
    for i in result:
        temp.append(i)
        temp.append(((-1)*i[0],i[1]))
        temp.append((i[0],(-1)*i[1]))
        temp.append(((-1)*i[0],(-1)*i[1]))
    for i in temp:
        final.append(((i[0]+cx),(i[1]+cy)))
    return final


def draw_curve(p_list, algorithm):
    """绘制曲线

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 曲线的控制点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'Bezier'和'B-spline'（三次均匀B样条曲线，曲线不必经过首末控制点）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    def cal_bs_point(u, p_list):
        temp = [-u ** 3 + 3 * u ** 2 - 3 * u + 1, 3 * u ** 3 - 6 * u ** 2 + 4, -3 * u ** 3 + 3 * u ** 2 + 3 * u + 1, u ** 3]
        res = 0.0
        for i in range(4):
            res += temp[i] * p_list[i]
        return res / 6
    
    du = 0.001
    result = []
    if algorithm == 'Bezier':
        n = len(p_list) - 1
        result.append(p_list[0])
        u = du
        #取样长度
        while u < 1:
            res = p_list.copy()
            for i in range(n):
                temp = []
                for j in range(len(res) - 1):
                    x1, y1 = res[j]
                    x2, y2 = res[j + 1]
                    temp.append([(1 - u) * x1 + u * x2, (1 - u) * y1 + u * y2])
                res = temp.copy()
            x, y = int(res[0][0] ), int(res[0][1])
            result.append([x, y])
            u += du
        result.append(p_list[-1])
    elif algorithm == 'B-spline':
        u = 0
        while u <= 1:
            for i in range(len(p_list) - 3):
                x_list = [point[0] for point in p_list[i:i + 4]]
                y_list = [point[1] for point in p_list[i:i + 4]]
                result.append([int(cal_bs_point(u, x_list) ), int(cal_bs_point(u, y_list) )])
            u += du
    return result
    


def translate(p_list, dx, dy):
    """平移变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param dx: (int) 水平方向平移量
    :param dy: (int) 垂直方向平移量
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    temp_plist=[]
    for x,y in p_list:
        temp_plist.append((x+dx,y+dy))
    return temp_plist
    


def rotate(p_list, x, y, r):
    """旋转变换（除椭圆外）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 旋转中心x坐标
    :param y: (int) 旋转中心y坐标
    :param r: (int) 顺时针旋转角度（°）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """

    res=[]
    for item in p_list:
        x1=item[0]-x
        y1=item[1]-x
        x2 = x1 * math.cos(r / 180 * math.pi) - y1 * math.sin(r / 180 * math.pi)
        y2 = x1 * math.sin(r / 180 * math.pi) + y1 * math.cos(r / 180 * math.pi)
        res.append((int(x2+x),int(y2+y)))
    return res


def scale(p_list, x, y, s):
    """缩放变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 缩放中心x坐标
    :param y: (int) 缩放中心y坐标
    :param s: (float) 缩放倍数
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    res=[]
    for item in p_list:
        x1=(item[0]-x)*s
        y1=(item[1]-y)*s
        x2=x1+x
        y2=y1+y
        res.append((int(x2),int(y2)))
    return res
    


def clip(p_list, x_min, y_min, x_max, y_max, algorithm):
    """线段裁剪

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param x_min: 裁剪窗口左上角x坐标
    :param y_min: 裁剪窗口左上角y坐标
    :param x_max: 裁剪窗口右下角x坐标
    :param y_max: 裁剪窗口右下角y坐标
    :param algorithm: (string) 使用的裁剪算法，包括'Cohen-Sutherland'和'Liang-Barsky'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1]]) 裁剪后线段的起点和终点坐标
    """
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    if algorithm == 'Cohen-Sutherland':
        while True:
            code0_temp = []
            code1_temp = []
            # if y_max-y0>=0:
            #     code0_temp.append(0)
            # else:
            #     code0_temp.append(1)
            # if y0-y_min>=0:
            #     code0_temp.append(0)
            # else:
            #     code0_temp.append(1)
            # if x_max-x0>=0:
            #     code0_temp.append(0)
            # else:
            #     code0_temp.append(1)
            # if x0-x_min>=0:
            #     code0_temp.append(0)
            # else:
            #     code0_temp.append(1)
            
            # if y_max-y1>=0:
            #     code1_temp.append(0)
            # else:
            #     code1_temp.append(1)
            # if y1-y_min>=0:
            #     code1_temp.append(0)
            # else:
            #     code1_temp.append(1)
            # if x_max-x1>=0:
            #     code1_temp.append(0)
            # else:
            #     code1_temp.append(1)
            # if x1-x_min>=0:
            #     code1_temp.append(0)
            # else:
            #     code1_temp.append(1)

            code0_temp.append(0 if y_max - y0 >= 0 else 1)
            code0_temp.append(0 if y0 - y_min >= 0 else 1)
            code0_temp.append(0 if x_max - x0 >= 0 else 1)
            code0_temp.append(0 if x0 - x_min >= 0 else 1)
            code1_temp.append(0 if y_max - y1 >= 0 else 1)
            code1_temp.append(0 if y1 - y_min >= 0 else 1)
            code1_temp.append(0 if x_max - x1 >= 0 else 1)
            code1_temp.append(0 if x1 - x_min >= 0 else 1)
            code0 = code0_temp[0] * 8 + code0_temp[1] * 4 + code0_temp[2] * 2 + code0_temp[3]
            code1 = code1_temp[0] * 8 + code1_temp[1] * 4 + code1_temp[2] * 2 + code1_temp[3]
            if code0 == 0 and code1 == 0:
                return [[int(x0), int(y0)], [int(x1), int(y1)]]
            elif code0 & code1 != 0:
                return []
            else:
                if code0 == 0:
                    x0, y0, x1, y1 = x1, y1, x0, y0
                    code0, code1 = code1, code0
                if code0 & 8 == 8:
                    u = (y_max - y0) / (y1 - y0)
                    x0 += u * (x1 - x0)
                    y0 = y_max
                elif code0 & 4 == 4:
                    u = (y_min - y0) / (y1 - y0)
                    x0 += u * (x1 - x0)
                    y0 = y_min
                elif code0 & 2 == 2:
                    u = (x_max - x0) / (x1 - x0)
                    y0 += u * (y1 - y0)
                    x0 = x_max
                elif code0 & 1 == 1:
                    u = (x_min - x0) / (x1 - x0)
                    y0 += u * (y1 - y0)
                    x0 = x_min
    elif algorithm == 'Liang-Barsky':
        dx, dy = x1 - x0, y1 - y0
        p = [-dx, dx, -dy, dy]
        q = [x0 - x_min, x_max - x0, y0 - y_min, y_max - y0]
        u1 = 0
        u2 = 1
        for i in range(4):
            if p[i] == 0:
                if q[i] < 0:
                    return []
            elif p[i] < 0:
                u1 = max(u1, q[i] / p[i])
            else:
                u2 = min(u2, q[i] / p[i])
            if u1 > u2:
                return []
        return [[int(x0 + u1 * (x1 - x0) ), int(y0 + u1 * (y1 - y0) )],
                [int(x0 + u2 * (x1 - x0) ), int(y0 + u2 * (y1 - y0) )]]
    
