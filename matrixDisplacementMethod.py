import numpy as np
import math

import constructionPlot
from liteTools import MathTools, MiscTools
'''
本文件定义了三个类, 分别是Node(节点), Element(单元), Struction(结构)。
依赖关系而言, 单元依赖于节点, 结构依赖于单元与节点。
    比如单元有一些参数 (比如单元刚度矩阵) 依赖于节点的参数, 但节点的参数不依赖于单元。(节点的等效荷载算是特例, 所以被做成函数)
'''


class Node:
    '''
    节点: 节点是结构中突变点(截面性质改变/存在集中荷载/分布荷载突变处)
    '''

    def __init__(self, position=(0, 0), lock=(0, 0, 0), load=(0, 0, 0)):
        '''
        :params position: Tuple[float, float]: 节点坐标(x, y)
        :params lock: Tuple[int, int, int]: 节点被约束的分量
        :params load: Tuple[float, float, float]: 节点荷载(μ, ν, θ)
        '''
        # 基本性质
        self.position: tuple[float, float] = position
        self.unlock = tuple((0 if i else 1) for i in lock)
        self.load = load
        self.id = MiscTools.geneNumId()
        # 接口
        self.element: list[Element] = []
        # 解
        self.solution = {}

    def setProperties(self, x=None, y=None, load=None, lock=None):
        '''
        重设节点参数
        :params x: 节点x坐标
        :params y: 节点y坐标
        :params load: 节点荷载
        :params lock: 节点约束
        '''
        # 位置
        position = list(self.position)
        if x != None:
            position[0] = x
        if y != None:
            position[1] = y
        self.position = tuple(position)
        # 荷载
        if load != None:
            self.load = tuple(load)
        # 约束
        if lock != None:
            unlock = tuple((0 if i else 1) for i in lock)
            self.unlock = unlock
        return self

    def geneCalculationLoad(self):
        '''
        统计单元传来的等效荷载, 以及生成总节点荷载
        '''
        load = self.load
        load_equivalent = [0, 0, 0]
        for element in self.element:
            position = element.node.index(self)
            eleLoad = [-i for i in element.eleLoad]
            eleLoad = np.dot(element.matrix_coordTrans.T, eleLoad)
            eleLoad = eleLoad.tolist()[0][position*3:position*3+3]
            load_equivalent = [a+b for a, b in zip(load_equivalent, eleLoad)]
        load_calculate = [a+b for a, b in zip(load_equivalent, load)]
        return load_calculate

    def __str__(self):
        return f'Node {self.id}'


class Element:
    '''
    单元: 单元是两个节点之间的联系
    '''

    def __init__(self, node: tuple[Node, Node], elementEA=10**10, elementEI=1, q=0):
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
        self.id = MiscTools.geneNumId()
        # 节点单元联系
        for n in self.node:
            n.element.append(self)
        # 生成性质
        self.length = None
        self.unitVector = None
        self.rad = None
        self.ang = None
        self.eleLoad = None
        self.matrix_coordTrans = None
        self.matrix_elementL = None
        self.matrix_elementG = None
        self.update()
        # 解
        self.solution = {}

    def update(self):
        self.geneGeometricProperties()
        self.matrix_coordTrans = self.geneMatrix_globalToLocalCoordinateSystem()
        self.matrix_elementL = self.geneElementMatrix_localCoordinateSystem()
        self.matrix_elementG = self.geneElementMatrix_globalCoordinateSystem()
        self.eleLoad = self.geneElementLoad()

    def geneGeometricProperties(self):
        '''生成构件几何性质'''
        p1 = self.node[0].position
        p2 = self.node[1].position
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        l = (dx**2+dy**2)**(0.5)
        udx = dx/l
        udy = dy/l
        self.length = l
        self.unitVector = (udx, udy)
        self.rad = math.acos(udx) if udy >= 0 else 2*math.pi - math.acos(udx)
        self.ang = math.degrees(self.rad)

    def geneMatrix_globalToLocalCoordinateSystem(self):
        '''整体转局部坐标转换矩阵'''
        udx, udy = self.unitVector
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

    def geneElementLoad(self):
        q = self.q
        l = self.length
        self.eleLoad = [0, -q*l/2, -q*l**2/12, 0, -q*l/2, q*l**2/12]
        return self.eleLoad

    def __str__(self):
        return f'Element {self.id}'


class Struction:
    '''
    结构: 由若干有相互联系的节点和单元组成
    '''

    def __init__(self, firstNode: Node):
        '''
        输入结构中的任意节点, 自动寻找与之有联系的所有单元和节点
        '''
        # 生成结构元素列表
        self.id = MiscTools.geneNumId()
        self.firstNode = firstNode
        self.nodeList: list[Node] = []
        self.elementList: list[Element] = []
        self.isCalcultated = False
        self.update()

    def update(self):
        '''
        更新自身
        '''
        self.findStruction()

    def __str__(self):
        return f'Struction {self.id}'

    def findStruction(self):
        '''
        从初始节点开始寻找与其连接的节点和单元, 更新到自身的nodeList和elementList中
        '''
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
        # 更新单元
        for e in self.elementList:
            e.update()
        # 生成总钢矩阵
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
        # 节点荷载向量 (含等效节点荷载)
        self.loadArray = np.zeros(matrixSize)
        for node in self.nodeList:
            ind = self.nodeList.index(node)*3
            for i, v in enumerate(node.geneCalculationLoad()):
                self.loadArray[ind+i] = v
        self.loadArray = np.dot(self.loadArray, self.matrix_constraint)
        # 解出位移矩阵
        self.matrix_deformation = np.dot(np.linalg.inv(
            self.matrix_totalStiffness), self.loadArray)
        self.isCalcultated = True
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
            elementFullForceInLocal = [
                a+b for a, b in zip(elementForceInLocal, element.eleLoad)]
            element.solution[self.id] = {'fullForce': elementFullForceInLocal}
        return self

    def printImage(self, scale=(0.1, 0.1, 0.1), printForce=(True, True, True), outputType=0):
        '''
        根据计算结果绘制内力图
        :params scale: tuple[float, float, float]: 轴力 剪力 弯矩图的放大系数
        :params printForce: tuple[bool, bool, bool]: 是否绘制(轴力 剪力 弯矩)图
        '''
        if not self.isCalcultated:
            self.calculate()
        # 画图
        cplot = constructionPlot.StructionPlot(outputType)
        for element in self.elementList:
            '''逐个单元进行绘图'''
            f = element.solution[self.id]['fullForce']
            enp = element.node[0].position
            if printForce[0]:
                '''绘制轴力图'''
                cplot.plotAxialForce(
                    f[0], element.length, enp[0], enp[1], element.ang, scale=scale[0])
            if printForce[1]:
                '''绘制剪力图'''
                cplot.plotShearingForce(
                    f[1], -f[4], element.length, enp[0], enp[1], element.ang, scale=scale[1])
            if printForce[2]:
                '''绘制弯矩图'''
                cplot.plotBendingMoment(f[2], -f[5], element.q,
                                        element.length, enp[0], enp[1], element.ang, scale=scale[2])
        cplot.show()
        return self
