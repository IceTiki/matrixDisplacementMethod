## 说明

基于矩阵位移法的二维结构内力图绘制器。

![image-20220515190240478](readme.assets/image-20220515190240478.png)

### 功能

- [x] matplotlib内力图绘制
- [x] 精确的节点位置
- [x] 自定义构件的EA和EI
- [x] 在构件上布置均布荷载

### 特性

* 坐标系是数学上常用的右手坐标系（结构力学书上常用左手坐标系）
  * 力在x, y轴正方向为正
  * 力矩逆时针方向为正

### 尚未修复

* 内力图的正负和绘图方向在某些情况下会绘图错误

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
* [Python与Ansys apdl有限元系列二：矩阵位移法计算桁架结构_王.伟的博客-CSDN博客](https://blog.csdn.net/weixin_43717845/article/details/105515372)
* [01 -- 矩阵位移法 - 知乎](https://zhuanlan.zhihu.com/p/57871511)

