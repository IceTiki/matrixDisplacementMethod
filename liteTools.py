import math
import random


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
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        return (dx**2+dy**2)**(0.5)
