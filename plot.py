import matplotlib.pyplot as plt

from .utils import MathTools, MiscTools


class FuncLib:
    @staticmethod
    def quadratic_function_by_3points(x_1, y_1, x_2, y_2, x_3, y_3):
        """
        输入二次函数三个点, 输出该二次函数
        :return: function: 返回函数f(x)
        """
        b = (
            (x_2**2 - x_3**2) * (y_1 - y_2) - (x_1**2 - x_2**2) * (y_2 - y_3)
        ) / ((x_2**2 - x_3**2) * (x_1 - x_2) - (x_1**2 - x_2**2) * (x_2 - x_3))
        a = (y_1 - y_2 - b * (x_1 - x_2)) / (x_1**2 - x_2**2)
        c = y_1 - a * x_1**2 - b * x_1
        return lambda x: a * x**2 + b * x + c

    @staticmethod
    def linear_function_by_2points(x_1, y_1, x_2, y_2):
        """
        输一次函数三个点, 输出该一次函数
        :return: function: 返回函数f(x)
        """
        return lambda x: (y_2 - y_1) / (x_2 - x_1) * (x - x_1) + y_1


class FunctionCurve:
    def __init__(
        self,
        func=lambda x: 0,
        func_a=-100,
        func_b=-100,
        func_step=0.1,
        local_x=0,
        local_y=0,
        local_a=0,
        y_scale=1,
    ):
        """
        :param func: 曲线的函数
        :param func_a: 函数曲线范围起始
        :param func_b: 函数曲线范围终点
        :param func_step: 函数曲线绘制步长
        :param local_x: 局部坐标系原点在整体坐标系的横坐标
        :param local_y: 局部坐标系原点在整体坐标系的纵坐标
        :param local_a: 局部坐标系的转角(角度)
        :param y_scale: 局部坐标系y轴比例
        """
        # 基本参数
        self.func = func
        self.func_a = func_a
        self.func_b = func_b
        self.func_step = func_step
        self.local_x = local_x
        self.local_y = local_y
        self.local_a = local_a
        self.y_scale = y_scale
        # 生成函数曲线
        self.func_x = MathTools.gene_sequence(self.func_a, self.func_b, self.func_step)
        self.func_y = [self.func(i) for i in self.func_x]
        # 寻找极值点和函数边缘
        self.find_extreme_points()
        # 寻找最值点
        extreme_points = self.extreme_points[:]
        extreme_points.sort(key=lambda x: x[1])
        self.func_max = extreme_points[-1]
        self.func_min = extreme_points[0]

    def set_local_coordinate_system(self, local_x=0, local_y=0, local_a=0, y_scale=1):
        """
        设置局部坐标系性质
        :param local_x: 局部坐标系原点在整体坐标系的横坐标
        :param local_y: 局部坐标系原点在整体坐标系的纵坐标
        :param local_a: 局部坐标系的转角(角度)
        :param y_scale: 局部坐标系y轴比例
        """
        self.local_x = local_x
        self.local_y = local_y
        self.local_a = local_a
        self.y_scale = y_scale
        return self

    def gene_local_curve(self):
        """
        生成函数曲线
        :return (curve_x, curve_y): 函数曲线各点的x,y值
        """
        self.local_curve_x, self.local_curve_y = [], []
        for x, y in zip(self.func_x, self.func_y):
            local_point = self.coord_local_to_global(x, y)
            self.local_curve_x.append(local_point[0])
            self.local_curve_y.append(local_point[1])
        return self.local_curve_x, self.local_curve_y

    def gene_local_x_axis(self):
        """
        生成局部坐标系X轴
        """
        line = [self.func_a, self.func_b]
        linx_x = [self.coord_local_to_global(x, 0)[0] for x in line]
        line_y = [self.coord_local_to_global(x, 0)[1] for x in line]
        return linx_x, line_y

    def gene_local_y_axis(self):
        """
        生成局部坐标系Y轴
        """
        line = [self.func_min[1], self.func_max[1]]
        line_x = [self.coord_local_to_global(0, y)[0] for y in line]
        line_y = [self.coord_local_to_global(0, y)[1] for y in line]
        return line_x, line_y

    def find_extreme_points(self, sensitivity=2):
        """
        寻找函数极值点和函数边缘
        """
        self.extreme_points = []
        for i, x in enumerate(self.func_x):
            y_1 = self.func(x)
            y_0 = y_1 if i == 0 else self.func(self.func_x[i - 1])
            y_2 = y_1 if i == len(self.func_x) - 1 else self.func(self.func_x[i + 1])
            delta_y1 = y_1 - y_0
            delta_y2 = y_2 - y_1
            if i == 0 or i == len(self.func_x) - 1:
                self.extreme_points.append((x, y_1))
            elif delta_y1 * delta_y2 < 0:
                self.extreme_points.append((x, y_1))
            elif delta_y1 * delta_y2 == 0 and abs(delta_y1 + delta_y2) > 10 ** (
                -sensitivity
            ):
                self.extreme_points.append((x, y_1))
        return self.extreme_points

    def coord_local_to_global(self, x, y):
        """
        局部坐标转整体坐标
        :param x: 局部x坐标
        :param y: 局部y坐标
        :return (x,y): 整体坐标
        """
        x_arr = MathTools.unit_vector(self.local_a)
        y_arr = MathTools.unit_vector(self.local_a + 90)
        gx = self.local_x + x * x_arr[0] + y * y_arr[0] * self.y_scale
        gy = self.local_y + x * x_arr[1] + y * y_arr[1] * self.y_scale
        return (gx, gy)


class StructionPlot:
    """
    结构内力图
    """

    def __init__(
        self,
        output_type=0,
        fig_size=(10, 10),
        img_output=("./pic", "pdf"),
        scale=(0.1, 0.1, 0.1),
        decimal=(2, 2, 2),
        show_digit=True,
    ):
        """
        :params output_type:
            0: 各内力图独占画布
            1: 各内力图共用画布
            2: 各内力图共用坐标系
        :params fig_size: 画布大小(单位: 英寸)
        :params img_output: 图片保存的位置和类型(如果为None则不保存)
        :param scale: 放大系数(轴力|剪力|弯矩)
        :params decimal: 数据标记精度(轴力|剪力|弯矩)
        :params show_digit: 是否显示数字
        """
        self.output_type = output_type
        self.img_output = img_output
        self.scale = scale
        self.decimal = decimal
        self.show_digit = show_digit
        self.figs = []
        self.axes = []
        if output_type == 0:
            for _ in range(3):
                fig = plt.figure(figsize=fig_size)
                fig.tight_layout()
                ax = fig.add_subplot()
                ax.set_aspect(1)
                self.figs.append(fig)
                self.axes.append(ax)
        elif output_type == 1:
            fig = plt.figure(figsize=fig_size)
            fig.tight_layout()
            self.figs.append(fig)
            for i in range(131, 134):
                ax = fig.add_subplot(i)
                ax.set_aspect(1)
                self.axes.append(ax)
        elif output_type == 2:
            fig = plt.figure(figsize=fig_size)
            fig.tight_layout()
            self.figs.append(fig)
            ax = fig.add_subplot()
            ax.set_aspect(1)
            self.axes = [ax, ax, ax]

    def plot_bending_moment(self, m1, m2, q, l, x, y, a):
        """
        输入若干参数, 绘制弯矩图
        :param m1: 左侧支座弯矩
        :param m2: 右侧支座弯矩
        :param q: 均布荷载
        :param l: 杆件长度
        :param x: 构件起点x坐标
        :param y: 构件起点y坐标
        :param a: 构件转角(角度)
        :return (x,y): 整体坐标
        """
        if q == 0:
            """没有均布荷载"""
            f = FuncLib.linear_function_by_2points(0, m1, l, m2)
        else:
            """有均布荷载"""
            f = FuncLib.quadratic_function_by_3points(
                0, m1, l, m2, l / 2, (1 / 8) * q * l**2 + (m1 + m2) / 2
            )
        fc = FunctionCurve(f, 0, l, l / 1000).set_local_coordinate_system(
            x, y, a, self.scale[2]
        )
        self.plot_on_ax(
            self.axes[2],
            fc,
            ("#000000", "#FF0000", "#E08389"),
            "Bending Moment Diagram",
            self.decimal[2],
        )

    def plot_shearing_force(self, v1, v2, l, x, y, a):
        """
        输入若干参数, 绘制剪力图
        :param v1: 左侧支座剪力
        :param v2: 右侧支座剪力
        :param l: 杆件长度
        :param x: 构件起点x坐标
        :param y: 构件起点y坐标
        :param a: 构件转角(角度)
        :return (x,y): 整体坐标
        """
        f = FuncLib.linear_function_by_2points(0, v1, l, v2)
        fc = FunctionCurve(f, 0, l, l / 1000).set_local_coordinate_system(
            x, y, a, self.scale[1]
        )
        self.plot_on_ax(
            self.axes[1],
            fc,
            ("#000000", "#0070C0", "#54C1F0"),
            "Shearing Force Diagram",
            self.decimal[1],
        )

    def plot_axial_force(self, n1, n2, l, x, y, a):
        """
        输入若干参数, 绘制轴力图
        :param n1: 杆件轴力
        :param n2: 杆件轴力
        :param l: 杆件长度
        :param x: 构件起点x坐标
        :param y: 构件起点y坐标
        :param a: 构件转角(角度)
        :return (x,y): 整体坐标
        """
        f = FuncLib.linear_function_by_2points(0, n1, l, n2)
        fc = FunctionCurve(f, 0, l, l / 1000).set_local_coordinate_system(
            x, y, a, self.scale[0]
        )
        self.plot_on_ax(
            self.axes[0],
            fc,
            ("#000000", "#138535", "#6BD089"),
            "Axial Force Diagram",
            self.decimal[0],
        )

    def plot_on_ax(
        self,
        ax: plt.Axes,
        fc: FunctionCurve,
        color=("#000000", "#085820", "#6BD089"),
        title="",
        decimal=2,
    ):
        """
        在指定ax上画内力图
        :params ax: 坐标系对象
        :params fc: 函数曲线对象
        :params color: 颜色样式(x轴, 内力, 标注)
        :params title: 坐标系标题
        :params decimal: 数据标记精度
        """
        # 增加坐标系标题
        ax.set_title(title)
        # 绘制函数曲线
        x, y = fc.gene_local_curve()
        ax.plot(x, y, color=color[1])
        # 绘制局部坐标x轴
        x, y = fc.gene_local_x_axis()
        ax.plot(x, y, color=color[0])
        # 标记极值点
        for x, y in fc.extreme_points:
            lx, ly = fc.coord_local_to_global(x, y)
            lx0, ly0 = fc.coord_local_to_global(x, 0)
            ax.plot([lx, lx0], [ly, ly0], color=color[2], alpha=0.5)
            if self.show_digit:
                ax.text(lx, ly, str(round(y, decimal)), color=color[2])

    def show(self):
        """
        显示(保存)图形
        """
        if self.img_output:
            if self.output_type == 0:
                filenames = (
                    "Axial Force Diagram",
                    "Shearing Force Diagram",
                    "Bending Moment Diagram",
                )
            elif self.output_type in (1, 2):
                filenames = ("Internal Force Diagram",)
            filenames = MiscTools.gene_file_folder(
                self.img_output[0], filenames, self.img_output[1]
            )
            for name, fig in zip(filenames, self.figs):
                fig.savefig(name)
        plt.show()
