import matplotlib.pyplot as plt

from liteTools import MathTools


class FuncLib:
    @staticmethod
    def quadraticFunctionBy3Points(x1, y1, x2, y2, x3, y3):
        '''
        输入二次函数三个点, 输出该二次函数
        :return: function: 返回函数f(x)
        '''
        b = ((x2**2-x3**2)*(y1-y2)-(x1**2-x2**2)*(y2-y3)) / \
            ((x2**2-x3**2)*(x1-x2)-(x1**2-x2**2)*(x2-x3))
        a = (y1-y2-b*(x1-x2))/(x1**2-x2**2)
        c = y1 - a*x1**2 - b*x1
        return lambda x: a*x**2 + b*x + c

    @staticmethod
    def linearFunctionBy2Points(x1, y1, x2, y2):
        '''
        输一次函数三个点, 输出该一次函数
        :return: function: 返回函数f(x)
        '''
        return lambda x: (y2-y1)/(x2-x1)*(x-x1) + y1


class FunctionCurve:
    def __init__(self, func=lambda x: 0, func_a=-100, func_b=-100, func_step=0.1, localX=0, localY=0, localA=0, yScale=1):
        '''
        :param func: 曲线的函数
        :param func_a: 函数曲线范围起始
        :param func_b: 函数曲线范围终点
        :param func_step: 函数曲线绘制步长
        :param localX: 局部坐标系原点在整体坐标系的横坐标
        :param localY: 局部坐标系原点在整体坐标系的纵坐标
        :param localA: 局部坐标系的转角(角度)
        :param yScale: 局部坐标系y轴比例
        '''
        # 基本参数
        self.func = func
        self.func_a = func_a
        self.func_b = func_b
        self.func_step = func_step
        self.localX = localX
        self.localY = localY
        self.localA = localA
        self.yScale = yScale
        # 生成函数曲线
        self.funcX = MathTools.geneSequence(
            self.func_a, self.func_b, self.func_step)
        self.funcY = [self.func(i) for i in self.funcX]
        # 寻找极值点和函数边缘
        self.findExtremePoints()
        # 寻找最值点
        ep = self.extremePoints[:]
        ep.sort(key=lambda x: x[1])
        self.funcMax = ep[-1]
        self.funcMin = ep[0]

    def setLocal(self, localX=0, localY=0, localA=0, yScale=1):
        '''
        设置局部坐标系性质
        :param localX: 局部坐标系原点在整体坐标系的横坐标
        :param localY: 局部坐标系原点在整体坐标系的纵坐标
        :param localA: 局部坐标系的转角(角度)
        :param yScale: 局部坐标系y轴比例
        '''
        self.localX = localX
        self.localY = localY
        self.localA = localA
        self.yScale = yScale
        return self

    def genelocalCurve(self):
        '''
        生成函数曲线
        :return (curveX, curveY): 函数曲线各点的x,y值
        '''
        self.localCurveX, self.localCurveY = [], []
        for x, y in zip(self.funcX, self.funcY):
            localPoint = self.localToGlobal(x, y)
            self.localCurveX.append(localPoint[0])
            self.localCurveY.append(localPoint[1])
        return self.localCurveX, self.localCurveY

    def genelocalXAxis(self):
        '''
        生成局部坐标系X轴
        '''
        line = [self.func_a, self.func_b]
        lineX = [self.localToGlobal(x, 0)[0] for x in line]
        lineY = [self.localToGlobal(x, 0)[1] for x in line]
        return lineX, lineY

    def genelocalYAxis(self):
        '''
        生成局部坐标系Y轴
        '''
        line = [self.funcMin[1], self.funcMax[1]]
        lineX = [self.localToGlobal(0, y)[0] for y in line]
        lineY = [self.localToGlobal(0, y)[1] for y in line]
        return lineX, lineY

    def findExtremePoints(self, sensitivity=2):
        '''
        寻找函数极值点和函数边缘
        '''
        self.extremePoints = []
        for i, x in enumerate(self.funcX):
            y1 = self.func(x)
            y0 = y1 if i == 0 else self.func(self.funcX[i-1])
            y2 = y1 if i == len(self.funcX)-1 else self.func(self.funcX[i+1])
            dy1 = y1 - y0
            dy2 = y2 - y1
            if i == 0 or i == len(self.funcX)-1:
                self.extremePoints.append((x, y1))
            elif dy1*dy2 < 0:
                self.extremePoints.append((x, y1))
            elif dy1*dy2 == 0 and abs(dy1+dy2) > 10**(-sensitivity):
                self.extremePoints.append((x, y1))
        return self.extremePoints

    def localToGlobal(self, x, y):
        '''
        局部坐标转整体坐标
        :param x: 局部x坐标
        :param y: 局部y坐标
        :return (x,y): 整体坐标
        '''
        xArr = MathTools.uniVector(self.localA)
        yArr = MathTools.uniVector(self.localA+90)
        gx = self.localX + x*xArr[0] + y*yArr[0]*self.yScale
        gy = self.localY + x*xArr[1] + y*yArr[1]*self.yScale
        return (gx, gy)


class StructionPlot:
    '''
    结构内力图
    '''

    def __init__(self, outputType=0):
        self.outputType = outputType
        if outputType == 0:
            self.figN = plt.figure(figsize=(12, 7))
            self.axN = self.figN.add_subplot()
            self.axN.set_aspect(1)
            self.figV = plt.figure(figsize=(12, 7))
            self.axV = self.figV.add_subplot()
            self.axV.set_aspect(1)
            self.figM = plt.figure(figsize=(12, 7))
            self.axM = self.figM.add_subplot()
            self.axM.set_aspect(1)
        elif outputType == 1:
            self.fig = plt.figure()
            self.fig.tight_layout()
            self.axN = self.fig.add_subplot(131)
            self.axN.set_aspect(1)
            self.axV = self.fig.add_subplot(132)
            self.axV.set_aspect(1)
            self.axM = self.fig.add_subplot(133)
            self.axM.set_aspect(1)

    def plotBendingMoment(self, m1, m2, q, l, x, y, a, scale=0.1):
        '''
        输入若干参数, 绘制弯矩图
        :param m1: 左侧支座弯矩
        :param m2: 右侧支座弯矩
        :param q: 均布荷载
        :param l: 杆件长度
        :param x: 构件起点x坐标
        :param y: 构件起点y坐标
        :param a: 构件转角(角度)
        :param scale: 弯矩放大系数
        :return (x,y): 整体坐标
        '''
        if q == 0:
            '''没有均布荷载'''
            f = FuncLib.linearFunctionBy2Points(0, m1, l, m2)
        else:
            '''有均布荷载'''
            f = FuncLib.quadraticFunctionBy3Points(
                0, m1, l, m2, l/2, (1/8)*q*l**2+(m1+m2)/2)
        fc = FunctionCurve(f, 0, l, l/1000).setLocal(x, y, a, scale)
        self.axPlot(self.axM, fc, ('#000000', '#FF0000',
                    '#E08389'), 'Bending Moment Diagram')

    def plotShearingForce(self, v1, v2, l, x, y, a, scale=0.1):
        '''
        输入若干参数, 绘制剪力图
        :param v1: 左侧支座剪力
        :param v2: 右侧支座剪力
        :param l: 杆件长度
        :param x: 构件起点x坐标
        :param y: 构件起点y坐标
        :param a: 构件转角(角度)
        :param scale: 放大系数
        :return (x,y): 整体坐标
        '''
        f = FuncLib.linearFunctionBy2Points(0, v1, l, v2)
        fc = FunctionCurve(f, 0, l, l/1000).setLocal(x, y, a, scale)
        self.axPlot(self.axV, fc, ('#000000', '#0070C0',
                    '#54C1F0'), 'Shearing Force Diagram')

    def plotAxialForce(self, n, l, x, y, a, scale=0.1):
        '''
        输入若干参数, 绘制轴力图
        :param n: 杆件轴力
        :param l: 杆件长度
        :param x: 构件起点x坐标
        :param y: 构件起点y坐标
        :param a: 构件转角(角度)
        :param scale: 放大系数
        :return (x,y): 整体坐标
        '''
        f = FuncLib.linearFunctionBy2Points(0, n, l, n)
        fc = FunctionCurve(f, 0, l, l/1000).setLocal(x, y, a, scale)
        self.axPlot(self.axN, fc, ('#000000', '#085820',
                    '#6BD089'), 'Axial Force Diagram')

    @staticmethod
    def axPlot(ax: plt.Axes, fc: FunctionCurve, color=('#000000', '#085820', '#6BD089'), title=''):
        '''
        在指定ax上画内力图
        :params ax: 坐标系对象
        :params fc: 函数曲线对象
        :params color: 颜色样式(x轴, 内力, 标注)
        :params title: 坐标系标题
        '''
        # 增加坐标系标题
        ax.set_title(title)
        # 绘制函数曲线
        x, y = fc.genelocalCurve()
        ax.plot(x, y, color=color[1])
        # 绘制局部坐标x轴
        x, y = fc.genelocalXAxis()
        ax.plot(x, y, color=color[0])
        # 标记极值点
        for x, y in fc.extremePoints:
            lx, ly = fc.localToGlobal(x, y)
            lx0, ly0 = fc.localToGlobal(x, 0)
            ax.plot([lx, lx0], [ly, ly0], color=color[2], alpha=0.5)
            ax.text(lx, ly, str(round(y, 3)), color=color[2])

    def show(self):
        '''
        显示图形
        '''
        plt.show()
