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
c1.printImage()
