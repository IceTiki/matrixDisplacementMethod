import locationCurve
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

    def nodeAllLoad(self):
        '''
        统计单元传来的等效荷载和节点本身的荷载
        :return: load: Tuple[float, float, float]: 节点荷载(μ, ν, θ)
        '''
        load = self.load
        # for element in self.element:
        #     q = element.q
        #     l = element.length
        #     eleLoad = np.array([0, q*l/2, q*l**2/12, 0, q*l/2, -q*l**2/12])
        #     eleLoad = np.dot(element.matrix_coordTrans.T, eleLoad)

        #     load = [a+b for a, b in zip(load, eleLoad)]
        # print([round(i,2) for i in load])
        return load

    def __str__(self):
        return f'Node {self.id}'


class Element:
    '''单元'''

    def __init__(self, node: tuple[Node, Node], elementEA=1, elementEI=1, q=0):
        '''
        :params node: tuple[Node, Node]: 单元连接的节点
        :params elementEA: float: 杆件弹性模量和截面面积的乘积
        :params elementEI: float: 杆件弹性模量和轴惯性矩的乘积
        :params q: float: 均布荷载
        '''
        # 基本性质
        self.node = node
        self.elementEA = elementEA
        self.elementEI = elementEI
        self.q = q
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
        self.matrix_coordTrans = self.geneMatrix_globalToLocalCoordinateSystem()
        self.matrix_elementL = self.geneElementMatrix_localCoordinateSystem()
        self.matrix_elementG = self.geneElementMatrix_globalCoordinateSystem()

    def geneMatrix_globalToLocalCoordinateSystem(self):
        '''整体转局部坐标转换矩阵'''
        p1 = self.node[0].position
        p2 = self.node[1].position
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        udx = dx/self.length
        udy = dy/self.length
        self.rad = math.acos(udx) if udy >= 0 else 2*math.pi - math.acos(udx)
        self.ang = math.degrees(self.rad)
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
        self.matrix_elementL = self.geneElementMatrix_localCoordinateSystem()
        self.matrix_coordTrans = self.geneMatrix_globalToLocalCoordinateSystem()
        self.matrix_elementG = np.dot(
            np.dot(self.matrix_coordTrans.T, self.matrix_elementL), self.matrix_coordTrans)
        return self.matrix_elementG

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
                if not v:
                    self.matrix_totalStiffness[ind+i][ind+i] = 1
        # 组合荷载矩阵
        self.loadArray = np.zeros(matrixSize)
        for node in self.nodeList:
            ind = self.nodeList.index(node)*3
            for i, v in enumerate(node.nodeAllLoad()):
                self.loadArray[ind+i] = v
        self.loadArray = np.dot(self.loadArray, self.matrix_constraint)
        # 解出位移矩阵
        self.matrix_deformation = np.dot(np.linalg.inv(
            self.matrix_totalStiffness), self.loadArray)
        return self

    def printImage(self):
        # 将变形保存到节点中
        for node in self.nodeList:
            ind = self.nodeList.index(node)*3
            node.solution[self.id] = {
                'deformation': self.matrix_deformation[ind: ind+3]}
        # 计算单元杆端应力
        for element in self.elementList:
            n1s = element.node[0].solution[self.id]
            n2s = element.node[1].solution[self.id]
            elementDeformation = np.concatenate(
                (n1s['deformation'], n2s['deformation']))
            elementForce = np.dot(element.matrix_elementG, elementDeformation)
            elementForceInLocal = np.dot(
                element.matrix_coordTrans, np.transpose(elementForce))
            elementForceInLocal = [
                float(elementForceInLocal[i][0]) for i in range(6)]
            element.solution[self.id] = {'force': elementForceInLocal}
        # 画图
        for element in self.elementList:
            f = element.solution[self.id]['force']
            enp = element.node[0].position
            locationCurve.plotBendingMoment(f[2], -f[5], 0,
                                            element.length, enp[0], enp[1], element.ang, 1)
            locationCurve.plotShearingForce(
                f[1], -f[4], element.length, enp[0], enp[1], element.ang, 1)
        locationCurve.show()
        return self


def lp(a, d=''):
    print((f'>>>>>>>>>>{d}<<<<<<<<<<\n' if d else '') + f'{np.round(a, 2)}')

# n1 = Node((0, 0), (1, 1, 1))
# n2 = Node((2, 0))
# n3 = Node((2, 1), load=(0, -1, 0))
# n4 = Node((3, 1), (1, 1, 1))
# e1 = Element((n1, n2),elementEA=99)
# e2 = Element((n2, n3),elementEA=99)
# e3 = Element((n3, n4),elementEA=99)
# c1 = Struction(n1)


n1 = Node((0, 0), (1, 1, 1))
n2 = Node((2, 0), (0, 0, 0))
n3 = Node((2, 1), load=(0, -1, 0))
n4 = Node((3, 1), (1, 1, 1))
e1 = Element((n1, n2), elementEA=10000)
e2 = Element((n2, n3), elementEA=10000)
e3 = Element((n3, n4), elementEA=10000)
c1 = Struction(n1)

c1.calculate()
c1.printImage()
