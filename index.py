from matrixDisplacementMethod import Node, Element, Struction


n1 = Node((0, 0), (1, 1, 0))
n2 = Node((1, 0))
n3 = Node((1, -1), (1, 1, 1))
n4 = Node((2, 0))
n5 = Node((3, 0))
n6 = Node((4, 1), (1, 1, 0))
n7 = Node((5, 0))
n8 = Node((6, 0))
n9 = Node((7, 2), load=(0, -10, 0))
e10 = Element((n1, n2), q=(0, -5))
e11 = Element((n2, n3), junctions=(1, 1, 0, 1, 1, 1))
e12 = Element((n2, n4), junctions=(1, 1, 1, 1, 0, 1))
e13 = Element((n4, n5), junctions=(1, 0, 1, 1, 1, 1))
e14 = Element((n5, n6))
e15 = Element((n5, n7))
e16 = Element((n7, n8), junctions=(1, 1, 1, 1, 1, 0))
e17 = Element((n7, n9), junctions=(1, 1, 0, 1, 1, 0))
e18 = Element((n8, n9), junctions=(1, 1, 0, 1, 1, 0))
c19 = Struction(n1)
c19.printImage(outputType=0)