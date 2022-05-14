import random


class Utils:
    @staticmethod
    def geneId(len_=32, dict_='0123456789abcdef'):
        return str().join(random.choices(dict_, k=len_))


class Node:
    def __init__(self, locked=False, position=(0, 0)):
        # 基本参数
        Struction.nodeList.append(self)
        self.id = Utils.geneId()
        self.locked = locked
        self.position = position
        # 其他参数
        self.bendingForceLoad = 0
        self.bendingForce = None

    def updateConstructonForce(self):
        conList = Struction.getConstruction(self.id)
        conList: list[Construction]
        # 更新节点弯矩
        if self.bendingForce == None:
            self.bendingForce = self.bendingForceLoad
        else:
            self.bendingForce = self.bendingForceLoad
            for con in conList:
                nodeOpposite = con.getOppositeNode(self.id)
                if nodeOpposite.locked or nodeOpposite.bendingForce == None:
                    continue
                self.bendingForce += con.bendingForce[con.getNodePosition(
                    nodeOpposite.id)]
        # 更新构件弯矩
        stiffnessSum = sum(i.stiffness for i in conList)
        for con in conList:
            # 更新相邻构件弯矩
            stiffWeight = con.stiffness/stiffnessSum*-0.5
            footBendingForce = -1*self.bendingForce * stiffWeight
            if con.node[0].position == (0,0) and con.node[1].position == (1,0):
                print(footBendingForce)
            con.bendingForce[con.getNodePosition(self.id)] = footBendingForce



class Construction:
    def __init__(self, stiffness: float, leftNode: Node, rightNode: Node, q=0):
        # 基本参数
        Struction.constructionList.append(self)
        self.id = (leftNode.id, rightNode.id)
        self.stiffness = stiffness
        self.node = (leftNode, rightNode)
        self.bendingForce = [0, 0]

    def getNodePosition(self, id=None):
        '''
        获取节点在构件中的位置
        :return int[0|1]: 节点位置
        '''
        if id:
            if id == self.id[0]:
                return 0
            elif id == self.id[1]:
                return 1
        raise Exception('构件中没有该节点')

    def getOppositeNode(self, id=None):
        '''
        获取构件中的另一个节点
        '''
        if id:
            if id == self.id[0]:
                return self.node[1]
            elif id == self.id[1]:
                return self.node[0]
        raise Exception('构件中没有该节点')

    def __str__(self):
        np0 = self.node[0].position
        np1 = self.node[1].position
        fb0 = self.node[0].bendingForceLoad + \
            self.bendingForce[0]*2  # + self.bendingForce[1]
        fb1 = self.node[1].bendingForceLoad + \
            self.bendingForce[1]*2  # + self.bendingForce[0]
        fb0, fb1 = map(lambda x: round(x, 2), [fb0, fb1])
        return f'[{np0}#{fb0}#{fb1}#{np1}]'


class Struction:
    nodeList = []
    nodeList: list[Node]
    constructionList = []
    constructionList: list[Construction]

    @staticmethod
    def getNode(id=None):
        if id:
            for n in Struction.nodeList:
                if n.id == id:
                    return n
        raise Exception('没有找到节点')

    @staticmethod
    def getConstruction(id=None):
        if id:
            conList = []
            for c in Struction.constructionList:
                if id in c.id:
                    conList.append(c)
            if conList:
                return conList
        raise Exception('没用找到构件')


if 1:
    '''生成节点矩阵'''
    bMatri = [
        [Node(), Node(), Node(), Node()],
        [Node(), Node(), Node(), Node()],
        [Node(), Node(), Node(), Node()],
        [Node(), Node(), Node(), Node()],
        [Node(), Node(), Node(), Node()],
        [Node(), Node(), Node(), Node()],
        [Node(True), Node(True), Node(True), Node(True)]
    ]

    for r in range(len(bMatri)):
        for c in range(len(bMatri[0])):
            bMatri[r][c].position = (r, c)

if 1:
    '''生成节点间构件'''
    # 柱线刚度
    stiffnessBetweenRowMatrix = [
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [1.3, 1.3, 1.3, 1.3]
    ]
    # 梁线刚度
    stiffnessBetweenColumnMatrix = [
        [0.32, 0.41, 0.32],
        [0.32, 0.41, 0.32],
        [0.32, 0.41, 0.32],
        [0.32, 0.41, 0.32],
        [0.32, 0.41, 0.32],
        [0.39, 0.5, 0.39],
        [None, None, None]
    ]

    for r in range(len(stiffnessBetweenRowMatrix)):
        for c in range(len(stiffnessBetweenRowMatrix[0])):
            i = stiffnessBetweenRowMatrix[r][c]
            if i == None:
                continue
            Construction(i, bMatri[r][c], bMatri[r+1][c])

    for r in range(len(stiffnessBetweenColumnMatrix)):
        for c in range(len(stiffnessBetweenColumnMatrix[0])):
            i = stiffnessBetweenColumnMatrix[r][c]
            if i == None:
                continue
            Construction(i, bMatri[r][c], bMatri[r][c+1])

if 1:
    '''节点所受荷载'''
    forceMatri = [
        [-62.55, 57.99, -57.99, 62.55],
        [-102.39, 98.1, -98.1, 102.39],
        [-102.39, 98.1, -98.1, 102.39],
        [-102.39, 98.1, -98.1, 102.39],
        [-102.39, 98.1, -98.1, 102.39],
        [-122.09, 125.39, -125.39, 122.09],
        [0, 0, 0, 0]
    ]

    for r in range(len(forceMatri)):
        for c in range(len(forceMatri[0])):
            i = forceMatri[r][c]
            bMatri[r][c].bendingForceLoad = i

if 1:
    '''更新节点'''
    for _ in range(4):
        for r in range(len(bMatri)):
            for c in range(len(bMatri[0])):
                bMatri[r][c].updateConstructonForce()

Struction.constructionList.sort(key=lambda x: x.node[0].position)

print('\n'.join([str(i) for i in Struction.constructionList]))

print('|'.join([str(i.bendingForce) for i in Struction.nodeList]))

print(len(Struction.constructionList))
