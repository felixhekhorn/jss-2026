from feynman import Diagram

gluon_style = dict(style="linear loopy", xamp=0.025, yamp=0.035, nloops=7)

diagram = Diagram(figsize=(10.0, 10.0))
inB = diagram.vertex(xy=(0, 1), marker="")
inG = diagram.vertex(xy=(0, 0), marker="")
vB = diagram.vertex(xy=(0.5, 0.7))
vG = diagram.vertex(xy=(0.5, 0.3))
outQ = diagram.vertex(xy=(1, 0.7), marker="")
outQb = diagram.vertex(xy=(1, 0.3), marker="")

b = diagram.line(inB, vB, style="wiggly")
g = diagram.line(inG, vG, **gluon_style)
outQ = diagram.line(outQ, vB)
t = diagram.line(vB, vG)
outQb = diagram.line(vG, outQb)

diagram.plot()
diagram.savefig("nlo-g.pdf")
