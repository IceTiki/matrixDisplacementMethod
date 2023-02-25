# matrixDisplacementMethod

**——基于矩阵位移法的二维结构内力图绘制器**

![image-20220515190240478](readme.assets/image-20220515190240478.png)

---

[<img src="https://github-readme-stats.vercel.app/api?username=IceTiki&show_icons=true" alt="Readme Card" style="zoom: 80%;" />](https://github.com/IceTiki)[<img src="https://github-readme-stats.vercel.app/api/pin/?username=IceTiki&repo=matrixDisplacementMethod" alt="Readme Card" style="zoom: 104%;" />](https://github.com/IceTiki/matrixDisplacementMethod)

### 功能

- [ ] matplotlib绘图
  - [x] 轴力
  - [x] 剪力
  - [x] 弯矩
  - [ ] 绘制荷载和约束
- [x] 结构定义
  - [x] 自由定义节点位置
  - [x] 自由定义构件的EA和EI
  - [x] 自由定义节点约束条件
  - [x] 自由定义节点与构件的间约束条件

    > 假设两端钢结点的情况写出等效节点荷载
    >
    > 再设置节点约束，求解矩阵，得到真正的等效节点荷载
  
- [x] 荷载定义
  - [x] 节点上的集中荷载
  - [x] 在构件上布置均布荷载
    - [x] 构件法向均布荷载
    - [x] 构件轴向均布荷载
  - [ ] 在单元上任意布置荷载
    - [ ] 函数曲线荷载
    - [ ] 集中荷载列表


### 特性

* 坐标系是数学上常用的右手坐标系（结构力学书上常用左手坐标系）
  * 力在x, y轴正方向为正
  * 力矩逆时针方向为正

### 问题

* 内力图的正负在特殊情况可能相反

## 教程

### 示例

> 1. 新建文件夹
> 2. 将代码文件夹`matrixDisplacementMethod`（注意，不要在解压时套了两层文件夹）移动到新建的文件夹内
> 3. 创建一个`*.py`文件，输入以下示例代码。并运行。

1. 定义**节点**(输入<u>坐标, 约束条件, 荷载......</u>)
2. 定义**构件**(输入<u>两个节点组成一个构件</u>)
3. 定义**结构**(输入<u>一个节点</u>，自动寻找所有与之相关的节点和构件)
4. **计算**和**绘图**

```python
from matrixDisplacementMethod import Node, Element, Struction

n1 = Node(position=(0, 0), constraints=(1, 1, 1)) # 定义节点位置、约束情况(x,y,转角)、荷载(见n3节点)
n2 = Node(position=(2, 0))
n3 = Node(position=(2, 1), load=(0, 0, 10))
n4 = Node(position=(3, 1), constraints=(1, 1, 1))
e1 = Element(node=(n1, n2), q=(0, -100)) # 定义杆件由哪两个节点相连、均布荷载的方向与大小(皆使用笛卡尔坐标系)
e2 = Element(node=(n2, n3))
e3 = Element(node=(n3, n4))
c1 = Struction(first_node=n1) # 选择结构内的一个节点(通过BFS建立整个结构)
c1.print_image(output_type=1) # 绘制内力图
```

## 开发与原理

### 开发

* 『相比弯矩分配法』<u>弯矩分配法不考虑水平和垂直位移的分配</u>，在非对称结构/荷载下与矩阵位移法算出的结构误差较大(可能达到10%)。对所有节点都增加水平和垂直位移的约束(考虑到对轴力图的影响以及杆件轴向变形比较小，垂直约束可以不用加)，与弯矩分配法的数值差可以减少两个左右的数量级。

* 『相比D值法』

  * D值法的抗侧移刚度是通过<u>本层梁柱</u>进行计算的。而矩阵位移法中，柱子受到<u>上层框架和地面的(转角)约束</u>越大，抗侧移刚度更大。越接近顶层抗侧移刚度越小。

  * 如果只在某两层之间施加剪力，在D值法中只有这两层之间存在层间水平位移。但在矩阵位移法中，因为<u>转角的变形协调</u>，会使得<u>其他层也产生层间水平位移</u>（可以达到存在剪力层的30%）。

    在风荷载计算中，用矩阵位移法计算时，因为其他层的影响，层间水平位移会比仅本层有剪力的情况大100%左右。但是最终算出的水平位移和D值法大多偏差只有5%。(换而言之，D值法的抗侧移刚度会考虑其他层剪力对本层的影响。但这个考虑是静态的，不会随着实际其他楼层情况而改变(不会实际计算其他层对本层抗侧移刚度的影响)。)(但首层(D值法偏大)和顶层(D值法偏小)会差30%\~40%)。
    
  * 虽然抗侧移刚度没有考虑其他层对本层端部转角约束的影响，但是在内力计算时，反弯点的位置会随着层数/相邻层高比而变化。所以内力图比侧移量精确得多。

* 『坐标手性』将坐标系从左手坐标系转到数学常用的右手坐标系不需要修改坐标转换矩阵/单元刚度矩阵等...因为力和位移的坐标系同时被转换了。

* 『通用性增强』一般后处理法只能处理两端刚节点的单元。但是只需要根据两端约束情况写适配的单元刚度矩阵，即可适配所有情况。(虽然两端6个分量总共有64种情况，但是轴向与法向和转角无关，分出来为4*16种情况。这16种情况中，只有5种是超静定。(静定结构不因位移产生内力))

* 『边界条件』边界条件使用置零法，将无意义分量对应行列(总钢矩阵)设为零，将主对角线值设为1，对应的荷载分量设为0。

### 原理

本项目基于<u>矩阵位移法(后处理法)</u>进行内力计算。

> 矩阵位移法是求弹性构件内力精确解的方法。
>
> 本项目使用后处理法，与先处理法相比
>
> * 优点
>   * 编码方便
>   * 稍微更好理解一些
> * 缺点
>   * 资料少，同济的《结构力学》没教
>   * 通用性差（但我用了一些特殊的方法修复了通用性）
>
> 矩阵位移法的核心是总刚度矩阵，在后处理法中，总钢矩阵的长和宽都是节点数的三倍。每一行每一列各代表节点的一个分量（$\nu, \mu, \theta$）。
>
> 单元是附加在节点上的联系，这种联系用单元刚度矩阵描述（两个节点之间相对位移产生的力），代表着两个节点之间位移与力之间的关系。在总钢矩阵中，单元刚度矩阵会叠加在单元对应的两个节点上(因为不同单元的对节点的作用是叠加的)。
>
> 组装完成的总钢矩阵代表着所有节点位移与力的关系，求解该矩阵与荷载，其实就是求解一个状态——所有节点都在合适的位置，在这个位置上拉扯单元产生的力恰好与荷载相同。

## 参考

* [第9章-弯矩分配法 - 豆丁网](https://www.docin.com/p-2175792518.html)

  > 原本想写弯矩分配法，但没成功
  
* [Python与Ansys apdl有限元系列二：矩阵位移法计算桁架结构_王.伟的博客-CSDN博客](https://blog.csdn.net/weixin_43717845/article/details/105515372)

  > 大概看了一眼，知道实现矩阵位移法的代码量不大

* [01 -- 矩阵位移法 - 知乎](https://zhuanlan.zhihu.com/p/57871511)

  > 这篇很有用，特别是关于边界条件的处理
  
* [结构力学.NET](http://www.jglx.net/)

  > 一个基于矩阵位移法的结构力学计算器，我编写的过程中，数据就是与这个计算器做对比。
  > (我当时不知道这个结构力学计算器可以精确编辑EA/EI和杆件长度，所以自己写了一个)
  >
  > (这个结构计算器在计算定向链杆连接时似乎会出错)
