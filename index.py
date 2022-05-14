from asyncio import constants
import numpy as np
import math
import random


class Utils:
    numUuid = 0

    @staticmethod
    def geneNumId():
        Utils.numUuid += 1
        return Utils.numUuid

    @staticmethod
    def geneId(len_=16, dict_='0123456789qwertyuiopasdfghjklzxcvbnQWERTYUIOPASDFGHJKLZXCVBNM'):
        return str().join(random.choices(dict_, k=len_))


class LiteMathTools:
    @staticmethod
    def distanceBetweenTwoPoint(p1: tuple[float, float], p2: tuple[float, float]):
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        return (dx**2+dy**2)**(0.5)


class Node:
    '''节点'''

    def __init__(self, position=(0, 0), lock=(0, 0, 0), load=(0, 0, 0)):
        '''
        :params position: Tuple[float, float]: 节点坐标(x, y)
        :params lock: Tuple[int, int, int]: 节点被约束的分量
        :params load: Tuple[float, float, float]: 节点荷载(μ, ν, θ)
        '''
        # 基本性质
        self.position = position
        self.unlock = tuple((0 if i else 1) for i in lock)
        self.load = load
        self.id = Utils.geneNumId()
        # 接口
        self.element: list[Element] = []
        # 解
        self.solution = {}

    def __str__(self):
        return f'Node {self.id}'


class Element:
    '''单元'''

    def __init__(self, node: tuple[Node, Node], elementEA=1, elementEI=1):
        '''
        :params node: tuple[Node, Node]: 单元连接的节点
        :params elementEA: float: 杆件弹性模量和截面面积的乘积
        :params elementEI: float: 杆件弹性模量和轴惯性矩的乘积
        '''
        # 基本性质
        self.node = node
        self.elementEA = elementEA
        self.elementEI = elementEI
        self.id = Utils.geneNumId()
        # 节点单元联系
        for n in self.node:
            n.element.append(self)
        # 生成性质
        self.length = None
        self.matrix_coordTrans = None
        self.matrix_elementL = None
        self.matrix_elementG = None
        self.update()
        # 解
        self.solution = {}

    def update(self):
        self.length = LiteMathTools.distanceBetweenTwoPoint(
            self.node[0].position, self.node[1].position)
        self.matrix_coordTrans = self.geneMatrix_localToGlobalCoordinateSystem()
        self.matrix_elementL = self.geneElementMatrix_localCoordinateSystem()
        self.matrix_elementG = self.geneElementMatrix_globalCoordinateSystem()

    def geneMatrix_localToGlobalCoordinateSystem(self):
        '''整体转局部坐标转换矩阵'''
        p1 = self.node[0].position
        p2 = self.node[1].position
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        udx = dx/self.length
        udy = dy/self.length
        matrix = [
            [udx, -udy, 0, 0, 0, 0],
            [udy, udx, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, udx, -udy, 0],
            [0, 0, 0, udy, udx, 0],
            [0, 0, 0, 0, 0, 1]
        ]
        return np.matrix(matrix)

    def geneElementMatrix_localCoordinateSystem(self):
        '''单元坐标系下的单元刚度矩阵'''
        ea = self.elementEA
        ei = self.elementEI
        l = self.length
        matrix = [
            [ea/l, 0, 0, -ea/l, 0, 0],
            [0, 12*ei/l**3, 6*ei/l**2, 0, -12*ei/l**3, 6*ei/l**2],
            [0, 6*ei/l**2, 4*ei/l, 0, -6*ei/l**2, 2*ei/l],
            [-ea/l, 0, 0, ea/l, 0, 0],
            [0, -12*ei/l**3, -6*ei/l**2, 0, 12*ei/l**3, -6*ei/l**2],
            [0, 6*ei/l**2, 2*ei/l, 0, -6*ei/l**2, 4*ei/l]
        ]
        matrix = np.matrix(matrix)
        return matrix

    def geneElementMatrix_globalCoordinateSystem(self):
        '''结构坐标系下的单元刚度矩阵'''
        return self.geneElementMatrix_localCoordinateSystem()*self.geneMatrix_localToGlobalCoordinateSystem().T

    def geneLockMatrix(self):
        '''锁定矩阵, 将矩阵/向量中, 被锁定的分量设为0'''
        n1 = self.node[0].unlock
        n2 = self.node[1].unlock
        meaningMatrix = [
            [n1[0], 0, 0, 0, 0, 0],
            [0, n1[1], 0, 0, 0, 0],
            [0, 0, n1[2], 0, 0, 0],
            [0, 0, 0, n2[0], 0, 0],
            [0, 0, 0, 0, n2[1], 0],
            [0, 0, 0, 0, 0, n2[2]]
        ]
        return np.matrix(meaningMatrix)

    def __str__(self):
        return f'Element {self.id}'


class Struction:
    '''整体结构'''

    def __init__(self, firstNode: Node):
        # 生成结构元素列表
        self.id = Utils.geneNumId()
        self.firstNode = firstNode
        self.nodeList: list[Node] = []
        self.elementList: list[Element] = []
        self.findStruction()

    def __str__(self):
        return f'Struction {self.id}'

    def findStruction(self):
        '''从初始节点开始寻找与其连接的节点和单元'''
        nodeList = [self.firstNode]
        elementList = []
        for n in nodeList:
            '''遍历节点列表'''
            for e in n.element:
                '''遍历该节点连接的单元列表'''
                if e not in elementList:
                    elementList.append(e)
                for n2 in e.node:
                    '''遍历该单元连接的节点列表'''
                    if n2 not in nodeList:
                        nodeList.append(n2)
        nodeList.sort(key=lambda x: x.id)
        elementList.sort(key=lambda x: x.id)
        self.nodeList = tuple(nodeList)
        self.elementList = tuple(elementList)
        return self

    def calculate(self):
        '''计算结构内力'''
        matrixSize = len(self.nodeList)*3
        self.matrix_totalStiffness = np.zeros((matrixSize, matrixSize))
        # 组合单钢矩阵
        for element in self.elementList:
            i1 = self.nodeList.index(element.node[0])*3
            i2 = self.nodeList.index(element.node[1])*3
            self.matrix_totalStiffness[i1:i1+3, i1:i1 +
                                       3] += element.matrix_elementG[0:3, 0:3]
            self.matrix_totalStiffness[i1:i1+3, i2:i2 +
                                       3] += element.matrix_elementG[0:3, 3:6]
            self.matrix_totalStiffness[i2:i2+3, i1:i1 +
                                       3] += element.matrix_elementG[3:6, 0:3]
            self.matrix_totalStiffness[i2:i2+3, i2:i2 +
                                       3] += element.matrix_elementG[3:6, 3:6]
        # 处理支座约束: 对应行列设0
        self.matrix_constraint = np.zeros((matrixSize, matrixSize))
        for node in self.nodeList:
            ind = self.nodeList.index(node)*3
            for i, v in enumerate(node.unlock):
                self.matrix_constraint[ind+i][ind+i] = v
        self.matrix_totalStiffness = np.dot(
            self.matrix_constraint, self.matrix_totalStiffness)
        self.matrix_totalStiffness = np.dot(
            self.matrix_totalStiffness, self.matrix_constraint)
        # 处理支座约束: 对应行列的主对角线元素设1
        for node in self.nodeList:
            ind = self.nodeList.index(node)*3
            for i, v in enumerate(node.unlock):
                self.matrix_totalStiffness[ind+i][ind+i] = 0 if v else 1
        # 组合荷载矩阵
        self.loadArray = np.zeros(matrixSize)
        for node in self.nodeList:
            ind = self.nodeList.index(node)*3
            for i, v in enumerate(node.load):
                self.loadArray[ind+i] = v
        self.loadArray = np.dot(self.loadArray, self.matrix_constraint)
        # 解出位移矩阵
        self.matrix_deformation = np.linalg.solve(
            self.matrix_totalStiffness, self.loadArray)
        return self

    def printImage(self):
        self.calculate()
        # 将变形保存到节点中
        for node in self.nodeList:
            ind = self.nodeList.index(node)*3
            node.solution[self.id] = {
                'deformation': self.matrix_deformation[ind: ind+3]}
        # 计算单元杆端应力
        for element in self.elementList:
            n1s = element.node[0].solution[self.id]
            n2s = element.node[1].solution[self.id]
            elementDeformation = np.concatenate((n1s['deformation'], n2s['deformation']))
            elementDeformation = np.dot(element.matrix_coordTrans, elementDeformation)
            elementForce = np.dot(element.matrix_elementL, np.transpose(elementDeformation))
            print(np.round(np.transpose(elementForce),3))

            
        return self

# n1 = Node((-1, -2), (0, 0, 0), (23, 11, -23))
# n2 = Node((3, 1), (1, 1, 0), (23, -124, -23))
# n3 = Node((7, 9), (0, 0, 0), (13, 17, -243))
# n4 = Node((-3, -9), (0, 0, 1), (23, 54, -33))
# e1 = Element((n1, n2), 100, 100)
# e2 = Element((n1, n3), 100, 100)
# e3 = Element((n2, n3), 100, 100)
# e4 = Element((n3, n4), 100, 100)
n1 = Node((0, 0), (1, 1, 0))
n2 = Node((2, 0))
n3 = Node((2, 1), (0, 0, 0), (0, 1, 0))
n4 = Node((3, 1), (0, 1, 0))

e1 = Element((n1, n2))
e2 = Element((n2, n3))
e3 = Element((n3, n4))

c1 = Struction(n1)

print([str(i) for i in c1.nodeList])
print([str(i) for i in c1.elementList])
a = c1.printImage()
print(np.round(a.matrix_totalStiffness, decimals=2))
# for i in a.nodeList:
#     print(i.solution)
# print(np.round(a.matrix_deformation, decimals=4))
# print(a.matrix_constraint)

