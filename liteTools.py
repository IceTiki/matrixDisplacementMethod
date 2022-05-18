import math
import random
import os
import numpy as np


class MiscTools:
    numUuid = 0

    @staticmethod
    def geneNumId():
        '''
        分配整数UUID
        '''
        MiscTools.numUuid += 1
        return MiscTools.numUuid

    @staticmethod
    def geneId(len_=16, dict_='0123456789qwertyuiopasdfghjklzxcvbnQWERTYUIOPASDFGHJKLZXCVBNM'):
        '''
        生成随机UUID
        '''
        return str().join(random.choices(dict_, k=len_))

    @staticmethod
    def geneFileFolder(folderDir='./pic', fileName=('a', 'b'), fileType='png'):
        '''
        :params folderDir: 文件夹路径
        :params fileName: 文件名元组
        :params fileType: 文件类型
        :return: 在文件夹里的文件路径元组
        '''
        if not os.path.isdir(folderDir):
            os.makedirs(folderDir)
        fileDir = []
        for i in fileName:
            fileDir.append(os.path.join(folderDir, f'{i}.{fileType}'))
        return tuple(fileDir)


class MathTools:
    @staticmethod
    def uniVector(degrees):
        '''
        输入角度, 输出单位向量
        :param degrees: 角度
        :return: tuple(x,y): 返回单位向量元组
        '''
        radians = math.radians(degrees)
        return (math.cos(radians), math.sin(radians))

    @staticmethod
    def geneSequence(a, b, step):
        '''
        生成等差数列 (范围a~b, 公差为step) 
        :param a: 数列首项
        :param b: 数列末项
        :param step: 数列公差
        :return tuple: 数列元组
        '''
        # 开始生成
        arr = []
        x = a
        while 1:
            arr.append(x)
            x += step
            if step == 0:
                raise Exception('step不应该为0')
            elif step*(x-b) > 0:
                arr.append(b)
                break
            elif step*(x-b) == 0:
                break
        return tuple(arr)

    @staticmethod
    def distanceBetweenTwoPoint(p1: tuple[float, float], p2: tuple[float, float]):
        '''计算两点间的间距'''
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        return (dx**2+dy**2)**(0.5)

    @staticmethod
    def findZeroCross(a: np.matrix):
        '''寻找矩阵中主对角线上的位置, 该位置所处的行列所有元素皆为0'''
        zeroRow = np.where(~a.any(axis=1))[0]
        zeroColumn = np.where(~a.any(axis=0))[1]
        zeroCross = np.intersect1d(zeroRow, zeroColumn)
        return zeroCross
