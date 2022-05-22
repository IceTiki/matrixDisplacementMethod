import numpy as np
import math

import constructionPlot
from liteTools import MathTools, MiscTools


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

    @property
    def calculationLoad(self):
        '''
        统计单元传来的等效荷载, 以及生成总节点荷载
        '''
        load = self.load
        load_equivalent = [0, 0, 0]
        for element in self.element:
            '''遍历单元'''
            position = element.node.index(self)
            nodeEquivalentLoads = element.nodeEquivalentLoads
            # 坐标转换
            nodeEquivalentLoads = np.dot(
                element.matrix_coordTrans.T, nodeEquivalentLoads.T)
            # 将单元对应本节点部分的节点荷载加上
            nodeEquivalentLoads = nodeEquivalentLoads.T.tolist()
            nodeEquivalentLoads = nodeEquivalentLoads[0][position*3:position*3+3]
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

    def __init__(self, node: tuple[Node, Node], elementEA=10**10, elementEI=1, q=(0, 0), junctions=(1, 1, 1, 1, 1, 1)):
        '''
        :params node: tuple[Node, Node]: 单元连接的节点
        :params elementEA: float: 杆件弹性模量和截面面积的乘积
        :params elementEI: float: 杆件弹性模量和轴惯性矩的乘积
        :params q: tuple[float, float]: 均布荷载(轴向, 法向)
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
        :params q: tuple[float, float]: 均布荷载(轴向, 法向)
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

    @property
    def matrix_elementL(self):
        '''单元坐标系下的单元刚度矩阵'''
        ea = self.elementEA
        ei = self.elementEI
        l = self.length
        junctions = self.junctions
        # 轴向
        if junctions[0] and junctions[3]:
            k11 = ea/l
        elif junctions[0] or junctions[3]:
            k11 = 0
        else:
            raise Exception(f'{self}在轴向缺乏约束')
        # 法向和转角
        junctionsYT = tuple(junctions[1:3] + junctions[4:6])
        junctionConditionTY = {
            (0, 0, 0, 0): '可变',
            (0, 0, 0, 1): '可变',
            (0, 0, 1, 0): '可变',
            (0, 0, 1, 1): '静定',
            (0, 1, 0, 0): '可变',
            (0, 1, 0, 1): '可变',
            (0, 1, 1, 0): '静定',
            (0, 1, 1, 1): '定固',
            (1, 0, 0, 0): '可变',
            (1, 0, 0, 1): '静定',
            (1, 0, 1, 0): '静定',
            (1, 0, 1, 1): '铰固',
            (1, 1, 0, 0): '静定',
            (1, 1, 0, 1): '固定',
            (1, 1, 1, 0): '固铰',
            (1, 1, 1, 1): '固固'
        }
        jc = junctionConditionTY[junctionsYT]
        if jc == '可变':
            raise Exception('{self}缺乏约束')
        elif jc == '静定':
            matrix = [
                [k11, 0, 0, -k11, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [-k11, 0, 0, k11, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0]
            ]
        elif jc == '固固':
            matrix = [
                [k11, 0, 0, -k11, 0, 0],
                [0, 12*ei/l**3, 6*ei/l**2, 0, -12*ei/l**3, 6*ei/l**2],
                [0, 6*ei/l**2, 4*ei/l, 0, -6*ei/l**2, 2*ei/l],
                [-k11, 0, 0, k11, 0, 0],
                [0, -12*ei/l**3, -6*ei/l**2, 0, 12*ei/l**3, -6*ei/l**2],
                [0, 6*ei/l**2, 2*ei/l, 0, -6*ei/l**2, 4*ei/l]
            ]
        elif jc == '固铰':
            matrix = [
                [k11, 0, 0, -k11, 0, 0],
                [0, 3*ei/l**3, 3*ei/l**2, 0, -3*ei/l**3, 0],
                [0, 3*ei/l**2, 3*ei/l, 0, -3*ei/l**2, 0],
                [-k11, 0, 0, k11, 0, 0],
                [0, -3*ei/l**3, -3*ei/l**2, 0, 3*ei/l**3, 0],
                [0, 0, 0, 0, 0, 0]
            ]
        elif jc == '铰固':
            matrix = [
                [k11, 0, 0, -k11, 0, 0],
                [0, 3*ei/l**3, 0, 0, -3*ei/l**3, 3*ei/l**2],
                [0, 0, 0, 0, 0, 0],
                [-k11, 0, 0, k11, 0, 0],
                [0, -3*ei/l**3, 0, 0, 3*ei/l**3, -3*ei/l**2],
                [0, 3*ei/l**2, 0, 0, -3*ei/l**2, 3*ei/l]
            ]
        elif jc == '固定':
            matrix = [
                [k11, 0, 0, -k11, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, ei/l, 0, 0, -ei/l],
                [-k11, 0, 0, k11, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, -ei/l, 0, 0, ei/l]
            ]
        elif jc == '定固':
            matrix = [
                [k11, 0, 0, -k11, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, ei/l, 0, 0, -ei/l],
                [-k11, 0, 0, k11, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, -ei/l, 0, 0, ei/l]
            ]
        return matrix

    @ property
    def matrix_elementL_endsFixed(self):
        '''单元坐标系下两端固接的单元刚度矩阵'''
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
        '''节点等效荷载(局部坐标系)'''
        qx = self.q[0]
        qy = self.q[1]
        l = self.length
        # 假设两端钢结的等效荷载
        nodeEquivalentLoads_0 = [0.5*qx*l, qy*l/2,
                                 qy*l**2/12, 0.5*qx*l, qy*l/2, -qy*l**2/12]
        nodeEquivalentLoads = nodeEquivalentLoads_0.copy()
        #
        matrix_element = self.matrix_elementL_endsFixed
        # 根据约束设对应行列为0
        setZeroMatrix = np.matrix(np.zeros((6, 6)))
        for i, v in enumerate(self.junctions):
            if not v:
                setZeroMatrix[i, i] = 1
        matrix_element = np.dot(
            np.dot(setZeroMatrix, matrix_element), setZeroMatrix)
        # 根据约束将对应主对角线约束设1, 荷载归0
        for i, v in enumerate(self.junctions):
            if v:
                matrix_element[i, i] = 1
                nodeEquivalentLoads[i] = 0
        # 解出实际约束条件下的位移
        nodeDeformation = np.dot(np.linalg.inv(
            matrix_element), np.array(nodeEquivalentLoads))
        nodeDeformation = nodeDeformation.reshape((6, 1))
        # 转换为实际的等效荷载
        nodeEquivalentLoads = np.dot(
            self.matrix_elementL_endsFixed, nodeDeformation)
        nodeEquivalentLoads = nodeEquivalentLoads.reshape((1, 6))
        nodeEquivalentLoads = np.add(-nodeEquivalentLoads,
                                     np.array(nodeEquivalentLoads_0))
        return nodeEquivalentLoads

    def misc(self):
        '''杂项代码堆放堆放'''
        # 1
        k11 = 0
        k22 = 0
        k23 = 0
        k26 = 0
        k32 = 0
        k62 = 0
        k33 = 0
        k63 = 0
        k36 = 0
        k66 = 0
        matrix = [
            [k11, 0, 0, -k11, 0, 0],
            [0, k22, k23, 0, -k22, k26],
            [0, k32, k33, 0, -k32, k36],
            [-k11, 0, 0, k11, 0, 0],
            [0, -k22, -k23, 0, k22, -k26],
            [0, k62, k63, 0, -k62, k66]
        ]

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
            for i, v in enumerate(node.calculationLoad):
                self.loadArray[ind+i] = v
        self.loadArray = np.dot(matrix_setConstraintsToZero, self.loadArray)
        # 解出位移矩阵
        self.matrix_deformation = np.dot(np.linalg.inv(
            self.matrix_totalStiffness), self.loadArray)
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
                a-b for a, b in zip(elementForceInLocal, element.nodeEquivalentLoads.tolist()[0])]
            element.solution[self.id] = {'fullForce': elementFullForceInLocal}
        self.isCalcultated = True
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
                    f[0], -f[3], element.length, enp[0], enp[1], element.ang)
            if printForce[1]:
                '''绘制剪力图'''
                cplot.plotShearingForce(
                    f[1], -f[4], element.length, enp[0], enp[1], element.ang)
            if printForce[2]:
                '''绘制弯矩图'''
                cplot.plotBendingMoment(
                    f[2], -f[5], element.q[1], element.length, enp[0], enp[1], element.ang)
        cplot.show()
        return self
