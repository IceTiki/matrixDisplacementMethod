# matrixDisplacementMethod

**——基于矩阵位移法的二维结构内力图绘制器**

[![Readme Card](https://github-readme-stats.vercel.app/api/pin/?username=IceTiki&repo=matrixDisplacementMethod)](https://github.com/IceTiki/matrixDisplacementMethod)

![image-20220515190240478](readme.assets/image-20220515190240478.png)

### 功能

- [x] matplotlib绘图
  - [x] 轴力
  - [x] 剪力
  - [x] 弯矩
  - [ ] 绘制荷载和约束
- [ ] 结构定义
  - [x] 自由定义节点位置
  - [x] 自由定义构件的EA和EI
  - [x] 自由定义节点约束条件
  - [ ] 自由定义节点与构件的间约束条件

    > 假设两端钢结点的情况写出等效节点荷载
    >
    > 再设置节点约束，求解矩阵，得到真正的等效节点荷载
  
- [ ] 荷载定义
  - [x] 节点上的集中荷载
  - [x] 在构件上布置均布荷载


### 特性

* 坐标系是数学上常用的右手坐标系（结构力学书上常用左手坐标系）
  * 力在x, y轴正方向为正
  * 力矩逆时针方向为正

### 问题

* 内力图的正负在特殊情况可能相反

## 教程

### 示例

1. 定义**节点**(输入<u>坐标, 约束条件, 荷载......</u>)
2. 定义**构件**(输入<u>两个节点组成一个构件</u>)
3. 定义**结构**(输入<u>一个节点</u>，自动寻找所有与之相关的节点和构件)
4. **计算**和**绘图**

```python
from matrixDisplacementMethod import Node, Element, Struction

n1 = Node((0, 0), (1, 1, 1))
n2 = Node((2, 0), (0, 0, 0))
n3 = Node((2, 1))
n4 = Node((3, 1), (1, 1, 1))
e1 = Element((n1, n2), q=-100)
e2 = Element((n2, n3))
e3 = Element((n3, n4))
c1 = Struction(n1)
c1.calculate()
c1.printImage(outputType=1)
```

## 参考文献

* [第9章-弯矩分配法 - 豆丁网](https://www.docin.com/p-2175792518.html)

  > 原本想写弯矩分配法，但没成功
* [Python与Ansys apdl有限元系列二：矩阵位移法计算桁架结构_王.伟的博客-CSDN博客](https://blog.csdn.net/weixin_43717845/article/details/105515372)

  > 大概看了一眼，知道实现矩阵位移法的代码量不大

* [01 -- 矩阵位移法 - 知乎](https://zhuanlan.zhihu.com/p/57871511)

  > 这篇很有用，特别是关于边界条件的处理

