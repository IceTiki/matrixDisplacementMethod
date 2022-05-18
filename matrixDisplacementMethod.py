import numpy as np
import math

import constructionPlot
from liteTools import MathTools, MiscTools
'''
本文件定义了三个类, 分别是Node(节点), Element(单元), Struction(结构)。
依赖关系而言, 单元依赖于节点, 结构依赖于单元与节点。
    比如单元有一些参数 (比如单元刚度矩阵) 依赖于节点的参数, 但节点的参数不依赖于单元。(节点的等效荷载算是特例, 所以被做成函数)
'''

# 注意，节点接收等效荷载的时候会反向


class Node:
    '''
    节点: 节点是结构中突变点(截面性质改变/存在集中荷载/分布荷载突变处)
    '''

    def __init__(self, position=(0, 0), constraints=(0, 0, 0), load=(0, 0, 0)):
        '''
        :params position: Tuple[float, float]: 节点坐标(x, y)
        :params constraints: Tuple[int, int, int]: 节点被支座约束的分量
        :params load: Tuple[float, float, float]: 节点荷载(μ, ν, θ)
        '''
        # 基本性质
        self.position: tuple[float, float] = position
        self.constraints = constraints
        self.load = load
        self.id = MiscTools.geneNumId()
        # 接口
        self.element: list[Element] = []
        # 解
        self.solution = {}

    def setProperties(self, x=None, y=None, load=None, constraints=None):
        '''
        重设节点参数
        :params x: 节点x坐标
        :params y: 节点y坐标
        :params load: 节点荷载
        :params constraints: 节点支座约束
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
        if constraints != None:
            self.constraints = constraints
        return self

    @property
    def unconstraintComponent(self):
        '''无支座约束的节点分量'''
        return tuple((0 if i else 1) for i in self.constraints)

    @property
    def bindingComponent(self):
        '''受支座或单元约束的非自由分量'''
        junctionConstraints = [0, 0, 0]
        directionVectors = []
        for element in self.element:
            '''遍历节点连接的单元'''
            nodeIndex = element.node.index(self)
            elementConstraints = element.junctions[nodeIndex *
                                                   3: (nodeIndex+1)*3]
            eleV = element.unitVector
            if elementConstraints[0]:
                '''如果有轴向约束'''
                directionVectors.append(eleV)
            if elementConstraints[1]:
                '''如果有法向约束'''
                directionVectors.append((-eleV[1], eleV[0]))
            if elementConstraints[2]:
                '''如果有转角约束'''
                junctionConstraints[2] = 1
        # 考虑支座约束
        if self.constraints[0]:
            directionVectors.append((1, 0))
        if self.constraints[1]:
            directionVectors.append((0, 1))
        if len(directionVectors) < 2:
            raise Exception(f'节点『{self}』缺少来自支座或单元的约束')
        directionVector_1 = directionVectors[0]
        for v in directionVectors[1:]:
            if np.cross(directionVector_1, v) != 0:
                '''有向量与向量一不平行, 则节点的x, y位移被约束'''
                junctionConstraints[0] = 1
                junctionConstraints[1] = 1
                break
        else:
            raise Exception(f'节点『{self}』缺少来自支座或单元的约束')
        return tuple(junctionConstraints)

    @property
    def meaningfulComponent(self):
        '''有意义的变量(没有被支座约束且与构件有连接)'''
        return tuple(1 if a and b else 0 for a, b in zip(self.unconstraintComponent, self.bindingComponent))

    def geneCalculationLoad(self):
        '''
        统计单元传来的等效荷载, 以及生成总节点荷载
        '''
        load = self.load
        load_equivalent = [0, 0, 0]
        for element in self.element:
            position = element.node.index(self)
            nodeEquivalentLoads = [-i for i in element.nodeEquivalentLoads]
            nodeEquivalentLoads = np.dot(
                element.matrix_coordTrans.T, nodeEquivalentLoads)
            nodeEquivalentLoads = nodeEquivalentLoads.tolist()[
                0][position*3:position*3+3]
            load_equivalent = [a+b for a,
                               b in zip(load_equivalent, nodeEquivalentLoads)]
        load_calculate = [a+b for a, b in zip(load_equivalent, load)]
        return load_calculate

    def __str__(self):
        return f'Node {self.id}'


class Element:
    '''
    单元: 单元是两个节点之间的联系
    '''

    def __init__(self, node: tuple[Node, Node], elementEA=10**10, elementEI=1, q=0, junctions=(1, 1, 1, 1, 1, 1)):
        '''
        :params node: tuple[Node, Node]: 单元连接的节点
        :params elementEA: float: 杆件弹性模量和截面面积的乘积
        :params elementEI: float: 杆件弹性模量和轴惯性矩的乘积
        :params q: float: 均布荷载
        :params junctions: 单元与节点之间的连接方式,
            分别代表单元端部与节点(x方向, y方向, 转角)是否绑定, 两个节点总共6个项目
        '''
        # 基本性质
        self.node = node
        self.junctions = junctions
        self.elementEA = elementEA
        self.elementEI = elementEI
        self.q = q
        self.id = MiscTools.geneNumId()
        # 节点单元联系
        for n in self.node:
            n.element.append(self)
        # 解
        self.solution = {}

    def setProperties(self, elementEA=None, elementEI=None, q=None, junctions=None):
        '''
        重设单元参数
        :params elementEA: 单元截面EA
        :params elementEI: 单元截面EI
        :params q: 均布荷载
        :params junctions: 单元与节点之间的连接方式,
            分别代表单元端部与节点(x方向, y方向, 转角)是否绑定, 两个节点总共6个项目
        '''
        if elementEA != None:
            self.elementEA = elementEA
        if elementEI != None:
            self.elementEI = elementEI
        if q != None:
            self.q = q
        if junctions != None:
            self.junctions = junctions

    @ property
    def rad(self):
        '''杆件转角(弧度)'''
        udx, udy = self.unitVector
        return math.acos(udx) if udy >= 0 else 2*math.pi - math.acos(udx)

    @ property
    def ang(self):
        '''杆件转角(弧度)'''
        return math.degrees(self.rad)

    @ property
    def unitVector(self):
        '''单位向量'''
        p1 = self.node[0].position
        p2 = self.node[1].position
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        l = (dx**2+dy**2)**(0.5)
        udx = dx/l
        udy = dy/l
        unitVector = (udx, udy)
        return unitVector

    @ property
    def length(self):
        '''杆件长度'''
        p1 = self.node[0].position
        p2 = self.node[1].position
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        l = (dx**2+dy**2)**(0.5)
        length = l
        return length

    @ property
    def matrix_coordTrans(self):
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

    @ property
    def matrix_elementL(self):
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

    @ property
    def matrix_elementG(self):
        '''结构坐标系下的单元刚度矩阵'''
        matrix_elementG = np.dot(
            np.dot(self.matrix_coordTrans.T, self.matrix_elementL), self.matrix_coordTrans)
        return matrix_elementG

    @ property
    def nodeEquivalentLoads(self):
        '''节点等效荷载'''
        q = self.q
        l = self.length
        nodeEquivalentLoads = [0, -q*l/2, -q*l**2/12, 0, -q*l/2, q*l**2/12]
        return nodeEquivalentLoads

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
        self.isCalcultated = False

    def __str__(self):
        return f'Struction {self.id}'

    @ property
    def nodeList(self):
        '''
        与初始节点有联系的节点列表
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
        nodeList: tuple[Node]
        return tuple(nodeList)

    @ property
    def elementList(self):
        '''
        与初始节点有联系的单元列表
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
        elementList.sort(key=lambda x: x.id)
        elementList: tuple[Element]
        return tuple(elementList)

    @ property
    def allMeaningfulComponents(self):
        '''
        整体刚度矩阵中需要计算的分量
        '''
        unconstraints = ()
        for node in self.nodeList:
            unconstraints += node.meaningfulComponent
        return unconstraints

    @ property
    def matrix_setConstraintsToZero(self):
        '''
        将受约束的结构节点分量设为0的矩阵
        (对向量需要前乘本矩阵)
        (对二维矩阵需要前乘本矩阵再后乘本矩阵)
        (总钢矩阵中, 受约束的节点分量行列被设为0(除了主对角线), 本矩阵用于设0)
        '''
        allUnconstraints = self.allMeaningfulComponents
        # 将无约束分量向量组装到矩阵主对角线上
        matrixSize = len(allUnconstraints)
        matrix_unconstraint = np.zeros((matrixSize, matrixSize))
        for i in range(matrixSize):
            matrix_unconstraint[i, i] = allUnconstraints[i]
        return matrix_unconstraint

    def calculate(self):
        '''计算结构内力'''
        nodeList = self.nodeList
        elementList = self.elementList
        # 生成总钢矩阵
        matrixSize = len(nodeList)*3
        self.matrix_totalStiffness = np.zeros((matrixSize, matrixSize))
        # 组合单钢矩阵
        for element in elementList:
            i1 = nodeList.index(element.node[0])*3
            i2 = nodeList.index(element.node[1])*3
            self.matrix_totalStiffness[i1:i1+3, i1:i1 +
                                       3] += element.matrix_elementG[0:3, 0:3]
            self.matrix_totalStiffness[i1:i1+3, i2:i2 +
                                       3] += element.matrix_elementG[0:3, 3:6]
            self.matrix_totalStiffness[i2:i2+3, i1:i1 +
                                       3] += element.matrix_elementG[3:6, 0:3]
            self.matrix_totalStiffness[i2:i2+3, i2:i2 +
                                       3] += element.matrix_elementG[3:6, 3:6]
        # 处理不需要计算的分量: 对应行列设0
        matrix_setConstraintsToZero = self.matrix_setConstraintsToZero
        self.matrix_totalStiffness = np.dot(np.dot(
            matrix_setConstraintsToZero, self.matrix_totalStiffness), matrix_setConstraintsToZero)
        # 处理不需要计算的分量: 对应行列的主对角线元素设1
        for i, v in enumerate(self.allMeaningfulComponents):
            if not v:
                self.matrix_totalStiffness[i][i] = 1
        # 节点荷载向量 (含等效节点荷载)
        self.loadArray = np.zeros(matrixSize)
        for node in nodeList:
            ind = nodeList.index(node)*3
            for i, v in enumerate(node.geneCalculationLoad()):
                self.loadArray[ind+i] = v
        self.loadArray = np.dot(self.loadArray, matrix_setConstraintsToZero)
        # 解出位移矩阵
        self.matrix_deformation = np.dot(np.linalg.inv(
            self.matrix_totalStiffness), self.loadArray)
        self.isCalcultated = True
        # 将变形保存到节点中
        for node in nodeList:
            ind = nodeList.index(node)*3
            node.solution[self.id] = {
                'deformation': self.matrix_deformation[ind: ind+3]}
        # 计算单元杆端应力
        for element in elementList:
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
                a+b for a, b in zip(elementForceInLocal, element.nodeEquivalentLoads)]
            element.solution[self.id] = {'fullForce': elementFullForceInLocal}
        return self

    def printImage(self, scale=(0.1, 0.1, 0.1), printForce=(True, True, True), outputType=0, figSize=(10, 10), picOutput=('./pic', 'pdf'), decimal=(2, 2, 2)):
        '''
        根据计算结果绘制内力图
        :params scale: tuple[float, float, float]: 轴力 剪力 弯矩图的放大系数
        :params printForce: tuple[bool, bool, bool]: 是否绘制(轴力 剪力 弯矩)图
        :params outputType:
            0: 各内力图独占画布
            1: 各内力图共用画布
        :params figSize: 画布大小(单位: 英寸)
        :params picOutput: 图片保存的位置和类型(如果为None则不保存)
        :params decimal: 数据标记精度(轴力|剪力|弯矩)
        '''
        if not self.isCalcultated:
            self.calculate()
        # 画图
        cplot = constructionPlot.StructionPlot(
            outputType, figSize=figSize, picOutput=picOutput, scale=scale, decimal=decimal)
        for element in self.elementList:
            '''逐个单元进行绘图'''
            f = element.solution[self.id]['fullForce']
            enp = element.node[0].position
            if printForce[0]:
                '''绘制轴力图'''
                cplot.plotAxialForce(
                    f[0], element.length, enp[0], enp[1], element.ang)
            if printForce[1]:
                '''绘制剪力图'''
                cplot.plotShearingForce(
                    f[1], -f[4], element.length, enp[0], enp[1], element.ang)
            if printForce[2]:
                '''绘制弯矩图'''
                cplot.plotBendingMoment(
                    f[2], -f[5], element.q, element.length, enp[0], enp[1], element.ang)
        cplot.show()
        return self
