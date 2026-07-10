from feynman import Diagram

gluon_style = dict(style="linear loopy", xamp=0.025, yamp=0.035, nloops=7)

diagram = Diagram(figsize=(10.0, 10.0))
inB = diagram.vertex(xy=(0, 1), marker="")
inQ = diagram.vertex(xy=(0, 0), marker="")
vB = diagram.vertex(xy=(0.5, 0.7))
vg = diagram.vertex(xy=(0.5, 0.3))
outQ = diagram.vertex(xy=(1, 0.7), marker="")
outG = diagram.vertex(xy=(1, 0.3), marker="")

b = diagram.line(inB, vB, style="wiggly")
inQ = diagram.line(inQ, vg)
diagram.line(vg, vB)
g = diagram.line(vg, outG, **gluon_style)
outQ = diagram.line(vB, outQ)

diagram.plot()
diagram.savefig("nlo-q-rt.pdf")
