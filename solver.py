import numpy as np
import math

from . import plot
from .utils import MiscTools


class Node:
    """
    节点: 节点是结构中突变点(截面性质改变/存在集中荷载/分布荷载突变处)
    """

    def __init__(self, position=(0, 0), constraints=(0, 0, 0), load=(0, 0, 0)):
        """
        :params position: Tuple[float, float]: 节点坐标(x, y)
        :params constraints: Tuple[int, int, int]: 节点被支座约束的分量
        :params load: Tuple[float, float, float]: 节点荷载(μ, ν, θ)
        """
        # 基本性质
        self.position: tuple[float, float] = position
        self.constraints = constraints
        self.load = load
        self.id = MiscTools.gene_int_uuid()
        # 接口
        self.element: list[Element] = []
        # 解
        self.solution = {}

    def set_properties(self, x=None, y=None, load=None, constraints=None):
        """
        重设节点参数
        :params x: 节点x坐标
        :params y: 节点y坐标
        :params load: 节点荷载
        :params constraints: 节点支座约束
        """
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
    def unconstraint_component(self):
        """无支座约束的节点分量"""
        return tuple((0 if i else 1) for i in self.constraints)

    @property
    def binding_component(self):
        """受支座或单元约束的非自由分量"""
        junction_constraints = [0, 0, 0]
        direction_vectors = []
        for element in self.element:
            """遍历节点连接的单元"""
            node_index = element.node.index(self)
            element_constraints = element.junctions[
                node_index * 3 : (node_index + 1) * 3
            ]
            element_unit_vector = element.unit_vector
            if element_constraints[0]:
                """如果有轴向约束"""
                direction_vectors.append(element_unit_vector)
            if element_constraints[1]:
                """如果有法向约束"""
                direction_vectors.append(
                    (-element_unit_vector[1], element_unit_vector[0])
                )
            if element_constraints[2]:
                """如果有转角约束"""
                junction_constraints[2] = 1
        # 考虑支座约束
        if self.constraints[0]:
            direction_vectors.append((1, 0))
        if self.constraints[1]:
            direction_vectors.append((0, 1))
        if len(direction_vectors) < 2:
            pass  # f'节点『{self}』缺少来自支座或单元的约束'
        else:
            direction_vector_1 = direction_vectors[0]
            for v in direction_vectors[1:]:
                if np.cross(direction_vector_1, v) != 0:
                    """有向量与向量一不平行, 则节点的x, y位移被约束"""
                    junction_constraints[0] = 1
                    junction_constraints[1] = 1
                    break
            else:
                pass  # f'节点『{self}』缺少来自支座或单元的约束'
        return tuple(junction_constraints)

    @property
    def meaningful_component(self):
        """有意义的变量(没有被支座约束且与构件有连接)"""
        return tuple(
            1 if a and b else 0
            for a, b in zip(self.unconstraint_component, self.binding_component)
        )

    @property
    def calculation_load(self):
        """
        统计单元传来的等效荷载, 以及生成总节点荷载
        """
        load = self.load
        load_equivalent = [0, 0, 0]
        for element in self.element:
            # 遍历单元
            position = element.node.index(self)
            node_equivalent_loads = element.node_equivalent_loads
            # 坐标转换
            node_equivalent_loads = np.dot(
                element.matrix_coord_trans.T, node_equivalent_loads.T
            )
            # 将单元对应本节点部分的节点荷载加上
            node_equivalent_loads = node_equivalent_loads.T.tolist()
            node_equivalent_loads = node_equivalent_loads[0][
                position * 3 : position * 3 + 3
            ]
            load_equivalent = [
                a + b for a, b in zip(load_equivalent, node_equivalent_loads)
            ]
        load_calculate = [a + b for a, b in zip(load_equivalent, load)]
        return load_calculate

    def __str__(self):
        return f"Node {self.id}"


class Element:
    """
    单元: 单元是两个节点之间的联系
    """

    def __init__(
        self,
        node: tuple[Node, Node],
        element_EA=10**10,
        element_EI=1,
        q=(0, 0),
        junctions=(1, 1, 1, 1, 1, 1),
    ):
        """
        :params node: tuple[Node, Node]: 单元连接的节点
        :params element_EA: float: 杆件弹性模量和截面面积的乘积
        :params element_EI: float: 杆件弹性模量和轴惯性矩的乘积
        :params q: tuple[float, float]: 均布荷载(轴向, 法向)
        :params junctions: 单元与节点之间的连接方式,
            分别代表单元端部与节点(x方向, y方向, 转角)是否绑定, 两个节点总共6个项目
        """
        # 基本性质
        self.node = node
        self.junctions = junctions
        self.element_EA = element_EA
        self.element_EI = element_EI
        self.q = q
        self.id = MiscTools.gene_int_uuid()
        # 节点单元联系
        for n in self.node:
            n.element.append(self)
        # 解
        self.solution = {}

    def set_properties(self, element_EA=None, element_EI=None, q=None, junctions=None):
        """
        重设单元参数
        :params element_EA: 单元截面EA
        :params element_EI: 单元截面EI
        :params q: tuple[float, float]: 均布荷载(轴向, 法向)
        :params junctions: 单元与节点之间的连接方式,
            分别代表单元端部与节点(x方向, y方向, 转角)是否绑定, 两个节点总共6个项目
        """
        if element_EA != None:
            self.element_EA = element_EA
        if element_EI != None:
            self.element_EI = element_EI
        if q != None:
            self.q = q
        if junctions != None:
            self.junctions = junctions

    @property
    def rad(self):
        """杆件转角(弧度)"""
        udx, udy = self.unit_vector
        return math.acos(udx) if udy >= 0 else 2 * math.pi - math.acos(udx)

    @property
    def ang(self):
        """杆件转角(弧度)"""
        return math.degrees(self.rad)

    @property
    def unit_vector(self):
        """单位向量"""
        p1 = self.node[0].position
        p2 = self.node[1].position
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        l = (dx**2 + dy**2) ** (0.5)
        udx = dx / l
        udy = dy / l
        unitVector = (udx, udy)
        return unitVector

    @property
    def length(self):
        """杆件长度"""
        p1 = self.node[0].position
        p2 = self.node[1].position
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        l = (dx**2 + dy**2) ** (0.5)
        length = l
        return length

    @property
    def matrix_coord_trans(self):
        """整体转局部坐标转换矩阵"""
        udx, udy = self.unit_vector
        matrix = [
            [udx, udy, 0, 0, 0, 0],
            [-udy, udx, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, udx, udy, 0],
            [0, 0, 0, -udy, udx, 0],
            [0, 0, 0, 0, 0, 1],
        ]
        return np.matrix(matrix)

    @property
    def matrix_element_local(self):
        """单元坐标系下的单元刚度矩阵"""
        element_EA = self.element_EA
        ei = self.element_EI
        l = self.length
        junctions = self.junctions
        # 轴向
        if junctions[0] and junctions[3]:
            k11 = element_EA / l
        elif junctions[0] or junctions[3]:
            k11 = 0
        else:
            raise Exception(f"{self}在轴向缺乏约束")
        # 法向和转角
        junctions_Y_T = tuple(junctions[1:3] + junctions[4:6])
        junction_condition_Y_T_table = {
            (0, 0, 0, 0): "可变",
            (0, 0, 0, 1): "可变",
            (0, 0, 1, 0): "可变",
            (0, 0, 1, 1): "静定",
            (0, 1, 0, 0): "可变",
            (0, 1, 0, 1): "可变",
            (0, 1, 1, 0): "静定",
            (0, 1, 1, 1): "定固",
            (1, 0, 0, 0): "可变",
            (1, 0, 0, 1): "静定",
            (1, 0, 1, 0): "静定",
            (1, 0, 1, 1): "铰固",
            (1, 1, 0, 0): "静定",
            (1, 1, 0, 1): "固定",
            (1, 1, 1, 0): "固铰",
            (1, 1, 1, 1): "固固",
        }
        junction_type = junction_condition_Y_T_table[junctions_Y_T]
        if junction_type == "可变":
            raise Exception("{self}缺乏约束")
        elif junction_type == "静定":
            matrix = [
                [k11, 0, 0, -k11, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [-k11, 0, 0, k11, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
            ]
        elif junction_type == "固固":
            matrix = [
                [k11, 0, 0, -k11, 0, 0],
                [
                    0,
                    12 * ei / l**3,
                    6 * ei / l**2,
                    0,
                    -12 * ei / l**3,
                    6 * ei / l**2,
                ],
                [0, 6 * ei / l**2, 4 * ei / l, 0, -6 * ei / l**2, 2 * ei / l],
                [-k11, 0, 0, k11, 0, 0],
                [
                    0,
                    -12 * ei / l**3,
                    -6 * ei / l**2,
                    0,
                    12 * ei / l**3,
                    -6 * ei / l**2,
                ],
                [0, 6 * ei / l**2, 2 * ei / l, 0, -6 * ei / l**2, 4 * ei / l],
            ]
        elif junction_type == "固铰":
            matrix = [
                [k11, 0, 0, -k11, 0, 0],
                [0, 3 * ei / l**3, 3 * ei / l**2, 0, -3 * ei / l**3, 0],
                [0, 3 * ei / l**2, 3 * ei / l, 0, -3 * ei / l**2, 0],
                [-k11, 0, 0, k11, 0, 0],
                [0, -3 * ei / l**3, -3 * ei / l**2, 0, 3 * ei / l**3, 0],
                [0, 0, 0, 0, 0, 0],
            ]
        elif junction_type == "铰固":
            matrix = [
                [k11, 0, 0, -k11, 0, 0],
                [0, 3 * ei / l**3, 0, 0, -3 * ei / l**3, 3 * ei / l**2],
                [0, 0, 0, 0, 0, 0],
                [-k11, 0, 0, k11, 0, 0],
                [0, -3 * ei / l**3, 0, 0, 3 * ei / l**3, -3 * ei / l**2],
                [0, 3 * ei / l**2, 0, 0, -3 * ei / l**2, 3 * ei / l],
            ]
        elif junction_type == "固定":
            matrix = [
                [k11, 0, 0, -k11, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, ei / l, 0, 0, -ei / l],
                [-k11, 0, 0, k11, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, -ei / l, 0, 0, ei / l],
            ]
        elif junction_type == "定固":
            matrix = [
                [k11, 0, 0, -k11, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, ei / l, 0, 0, -ei / l],
                [-k11, 0, 0, k11, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, -ei / l, 0, 0, ei / l],
            ]
        return matrix

    @property
    def matrix_element_local_bothfixed(self):
        """单元坐标系下两端固接的单元刚度矩阵"""
        ea = self.element_EA
        ei = self.element_EI
        l = self.length
        matrix = [
            [ea / l, 0, 0, -ea / l, 0, 0],
            [
                0,
                12 * ei / l**3,
                6 * ei / l**2,
                0,
                -12 * ei / l**3,
                6 * ei / l**2,
            ],
            [0, 6 * ei / l**2, 4 * ei / l, 0, -6 * ei / l**2, 2 * ei / l],
            [-ea / l, 0, 0, ea / l, 0, 0],
            [
                0,
                -12 * ei / l**3,
                -6 * ei / l**2,
                0,
                12 * ei / l**3,
                -6 * ei / l**2,
            ],
            [0, 6 * ei / l**2, 2 * ei / l, 0, -6 * ei / l**2, 4 * ei / l],
        ]
        matrix = np.matrix(matrix)
        return matrix

    @property
    def matrix_element_global(self):
        """结构坐标系下的单元刚度矩阵"""
        matrix = np.dot(
            np.dot(self.matrix_coord_trans.T, self.matrix_element_local),
            self.matrix_coord_trans,
        )
        return matrix

    @property
    def node_equivalent_loads(self):
        """节点等效荷载(局部坐标系)"""
        qx, qy = self.q
        l = self.length
        # 假设两端钢结的等效荷载
        node_equivalent_loads_0 = [
            qx * l / 2,
            qy * l / 2,
            qy * l**2 / 12,
            qx * l / 2,
            qy * l / 2,
            -qy * l**2 / 12,
        ]
        node_equivalent_loads = node_equivalent_loads_0.copy()
        #
        matrix_element = self.matrix_element_local_bothfixed
        # 根据约束设对应行列为0
        set_zero_matrix = np.matrix(np.zeros((6, 6)))
        for i, v in enumerate(self.junctions):
            if not v:
                set_zero_matrix[i, i] = 1
        matrix_element = np.dot(
            np.dot(set_zero_matrix, matrix_element), set_zero_matrix
        )
        # 根据约束将对应主对角线约束设1, 荷载归0
        for i, v in enumerate(self.junctions):
            if v:
                matrix_element[i, i] = 1
                node_equivalent_loads[i] = 0
        # 解出实际约束条件下的位移
        node_deformation = np.dot(
            np.linalg.inv(matrix_element), np.array(node_equivalent_loads)
        )
        node_deformation = node_deformation.reshape((6, 1))
        # 转换为实际的等效荷载
        node_equivalent_loads = np.dot(
            self.matrix_element_local_bothfixed, node_deformation
        )
        node_equivalent_loads = node_equivalent_loads.reshape((1, 6))
        node_equivalent_loads = np.add(
            -node_equivalent_loads, np.array(node_equivalent_loads_0)
        )
        return node_equivalent_loads

    def misc(self):
        """杂项代码堆放堆放"""
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
            [0, k62, k63, 0, -k62, k66],
        ]

    def __str__(self):
        return f"Element {self.id}"


class Struction:
    """
    结构: 由若干有相互联系的节点和单元组成
    """

    def __init__(self, first_node: Node):
        """
        输入结构中的任意节点, 自动寻找与之有联系的所有单元和节点
        """
        # 生成结构元素列表
        self.id = MiscTools.gene_int_uuid()
        self.first_node = first_node
        self.is_calcultated = False

    def __str__(self):
        return f"Struction {self.id}"

    @property
    def node_list(self):
        """
        与初始节点有联系的节点列表
        """
        node_list = [self.first_node]
        element_list = []
        for n in node_list:
            """遍历节点列表"""
            for e in n.element:
                """遍历该节点连接的单元列表"""
                if e not in element_list:
                    element_list.append(e)
                for n2 in e.node:
                    """遍历该单元连接的节点列表"""
                    if n2 not in node_list:
                        node_list.append(n2)
        node_list.sort(key=lambda x: x.id)
        node_list: tuple[Node]
        return tuple(node_list)

    @property
    def element_list(self):
        """
        与初始节点有联系的单元列表
        """
        node_list = [self.first_node]
        element_list = []
        for n in node_list:
            """遍历节点列表"""
            for e in n.element:
                """遍历该节点连接的单元列表"""
                if e not in element_list:
                    element_list.append(e)
                for n2 in e.node:
                    """遍历该单元连接的节点列表"""
                    if n2 not in node_list:
                        node_list.append(n2)
        element_list.sort(key=lambda x: x.id)
        element_list: tuple[Element]
        return tuple(element_list)

    @property
    def all_meaningful_components(self):
        """
        整体刚度矩阵中需要计算的分量
        """
        unconstraints = ()
        for node in self.node_list:
            unconstraints += node.meaningful_component
        return unconstraints

    @property
    def matrix_set_constraints_to_zero(self):
        """
        将受约束的结构节点分量设为0的矩阵
        (对向量需要前乘本矩阵)
        (对二维矩阵需要前乘本矩阵再后乘本矩阵)
        (总钢矩阵中, 受约束的节点分量行列被设为0(除了主对角线), 本矩阵用于设0)
        """
        all_meaningful_components = self.all_meaningful_components
        # 将无约束分量向量组装到矩阵主对角线上
        matrix_size = len(all_meaningful_components)
        matrix_unconstraint = np.zeros((matrix_size, matrix_size))
        for i in range(matrix_size):
            matrix_unconstraint[i, i] = all_meaningful_components[i]
        return matrix_unconstraint

    def calculate(self):
        """计算结构内力"""
        node_list = self.node_list
        element_list = self.element_list
        # 生成总钢矩阵
        matrix_size = len(node_list) * 3
        self.matrix_total_stiffness = np.zeros((matrix_size, matrix_size))
        # 组合单钢矩阵
        for element in element_list:
            i1 = node_list.index(element.node[0]) * 3
            i2 = node_list.index(element.node[1]) * 3
            self.matrix_total_stiffness[
                i1 : i1 + 3, i1 : i1 + 3
            ] += element.matrix_element_global[0:3, 0:3]
            self.matrix_total_stiffness[
                i1 : i1 + 3, i2 : i2 + 3
            ] += element.matrix_element_global[0:3, 3:6]
            self.matrix_total_stiffness[
                i2 : i2 + 3, i1 : i1 + 3
            ] += element.matrix_element_global[3:6, 0:3]
            self.matrix_total_stiffness[
                i2 : i2 + 3, i2 : i2 + 3
            ] += element.matrix_element_global[3:6, 3:6]
        # 处理不需要计算的分量: 对应行列设0
        matrix_set_constraints_to_zero = self.matrix_set_constraints_to_zero
        self.matrix_total_stiffness = np.dot(
            np.dot(matrix_set_constraints_to_zero, self.matrix_total_stiffness),
            matrix_set_constraints_to_zero,
        )
        # 处理不需要计算的分量: 对应行列的主对角线元素设1
        for i, v in enumerate(self.all_meaningful_components):
            if not v:
                self.matrix_total_stiffness[i][i] = 1
        # 节点荷载向量 (含等效节点荷载)
        self.load_array = np.zeros(matrix_size)
        for node in node_list:
            ind = node_list.index(node) * 3
            for i, v in enumerate(node.calculation_load):
                self.load_array[ind + i] = v
        self.load_array = np.dot(matrix_set_constraints_to_zero, self.load_array)
        # 解出位移矩阵
        self.matrix_deformation = np.dot(
            np.linalg.inv(self.matrix_total_stiffness), self.load_array
        )
        # 将变形保存到节点中
        for node in node_list:
            ind = node_list.index(node) * 3
            node.solution[self.id] = {
                "deformation": self.matrix_deformation[ind : ind + 3]
            }
        # 计算单元杆端应力
        for element in element_list:
            n1s = element.node[0].solution[self.id]
            n2s = element.node[1].solution[self.id]
            element_deformation = np.concatenate(
                (n1s["deformation"], n2s["deformation"])
            )
            element_force = np.dot(element.matrix_element_global, element_deformation)
            element_force_local = np.dot(
                element.matrix_coord_trans, np.transpose(element_force)
            )
            element_force_local = [float(element_force_local[i][0]) for i in range(6)]
            element.solution[self.id] = {"force": element_force_local}
            element_full_force_in_local = [
                a - b
                for a, b in zip(
                    element_force_local, element.node_equivalent_loads.tolist()[0]
                )
            ]
            element.solution[self.id]["fullForce"] = element_full_force_in_local
            element.solution[self.id]["q"] = element.q

        self.is_calcultated = True
        return self

    def print_image(
        self,
        scale=(0.1, 0.1, 0.1),
        print_force=(True, True, True),
        output_type=0,
        fig_size=(10, 10),
        img_output=("./pic", "pdf"),
        decimal=(2, 2, 2),
    ):
        """
        根据计算结果绘制内力图
        :params scale: tuple[float, float, float]: 轴力 剪力 弯矩图的放大系数
        :params print_force: tuple[bool, bool, bool]: 是否绘制(轴力 剪力 弯矩)图
        :params output_type:
            0: 各内力图独占画布
            1: 各内力图共用画布
        :params fig_size: 画布大小(单位: 英寸)
        :params img_output: 图片保存的位置和类型(如果为None则不保存)
        :params decimal: 数据标记精度(轴力|剪力|弯矩)
        """
        if not self.is_calcultated:
            self.calculate()
        # 画图
        cplot = plot.StructionPlot(
            output_type,
            fig_size=fig_size,
            img_output=img_output,
            scale=scale,
            decimal=decimal,
        )
        for element in self.element_list:
            """逐个单元进行绘图"""
            f = element.solution[self.id]["fullForce"]
            enp = element.node[0].position
            if print_force[0]:
                """绘制轴力图"""
                cplot.plot_axial_force(
                    f[0], -f[3], element.length, enp[0], enp[1], element.ang
                )
            if print_force[1]:
                """绘制剪力图"""
                cplot.plot_shearing_force(
                    f[1], -f[4], element.length, enp[0], enp[1], element.ang
                )
            if print_force[2]:
                """绘制弯矩图"""
                cplot.plot_bending_moment(
                    f[2],
                    -f[5],
                    element.q[1],
                    element.length,
                    enp[0],
                    enp[1],
                    element.ang,
                )
        cplot.show()
        return self
