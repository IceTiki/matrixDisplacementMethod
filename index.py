from matrixDisplacementMethod import Node, Element, Struction

n1 = Node(position=(0, 0), constraints=(1, 1, 1))
n2 = Node(position=(2, 0))
n3 = Node(position=(2, 1), load=(0, 0, 10))
n4 = Node(position=(3, 1), constraints=(1, 1, 1))
e1 = Element(node=(n1, n2), q=-100)
e2 = Element(node=(n2, n3))
e3 = Element(node=(n3, n4))
c1 = Struction(firstNode=n1)
c1.printImage(outputType=1)
