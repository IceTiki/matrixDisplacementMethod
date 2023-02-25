import math
import random
import os
import numpy as np
import pathlib


class MiscTools:
    numUuid = 0

    @staticmethod
    def gene_int_uuid():
        """
        分配整数UUID
        """
        MiscTools.numUuid += 1
        return MiscTools.numUuid

    @staticmethod
    def gene_word_uuid(
        len_=16, dict_="0123456789qwertyuiopasdfghjklzxcvbnQWERTYUIOPASDFGHJKLZXCVBNM"
    ):
        """
        生成随机UUID
        """
        return str().join(random.choices(dict_, k=len_))

    @staticmethod
    def gene_file_folder(folder_dir="./pic", file_name=("a", "b"), file_type="png"):
        """
        :params folder_dir: 文件夹路径
        :params file_name: 文件名元组
        :params file_type: 文件类型
        :return: 在文件夹里的文件路径元组
        """
        folder = pathlib.Path(folder_dir)
        folder.mkdir(parents=True, exist_ok=True)
        file_dir = tuple(str(folder / f"{i}.{file_type}") for i in file_name)
        return file_dir


class MathTools:
    @staticmethod
    def unit_vector(degrees):
        """
        输入角度, 输出单位向量
        :param degrees: 角度
        :return: tuple(x,y): 返回单位向量元组
        """
        radians = math.radians(degrees)
        return (math.cos(radians), math.sin(radians))

    @staticmethod
    def gene_sequence(a, b, step):
        """
        生成等差数列 (范围a~b, 公差为step)
        :param a: 数列首项
        :param b: 数列末项
        :param step: 数列公差
        :return tuple: 数列元组
        """
        # 开始生成
        arr = []
        x = a
        while 1:
            arr.append(x)
            x += step
            if step == 0:
                raise Exception("step不应该为0")
            elif step * (x - b) > 0:
                arr.append(b)
                break
            elif step * (x - b) == 0:
                break
        return tuple(arr)

    @staticmethod
    def distance_between_2point(
        point1: tuple[float, float], point2: tuple[float, float]
    ):
        """计算两点间的间距"""
        delta_x = point2[0] - point1[0]
        delta_y = point2[1] - point1[1]
        return (delta_x**2 + delta_y**2) ** (0.5)

    @staticmethod
    def find_zero_cross(a: np.matrix):
        """寻找矩阵中主对角线上的位置, 该位置所处的行列所有元素皆为0"""
        zero_row = np.where(~a.any(axis=1))[0]
        zero_column = np.where(~a.any(axis=0))[1]
        zero_cross = np.intersect1d(zero_row, zero_column)
        return zero_cross
