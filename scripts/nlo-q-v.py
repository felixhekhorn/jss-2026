from feynman import Diagram

gluon_style = dict(
    style="elliptic loopy", ellipse_spread=-0.55, xamp=0.04, yamp=0.05, nloops=7
)

diagram = Diagram(figsize=(10.0, 10.0))
inB = diagram.vertex(xy=(0, 1), marker="")
inQ = diagram.vertex(xy=(0, 0), marker="")
vB = diagram.vertex(xy=(0.4, 0.5))
vGin = diagram.vertex(xy=(3.0 / 5.0 * 0.4, 0.3))
vGout = diagram.vertex(xy=(0.7, 0.5))
outQ = diagram.vertex(xy=(1, 0.5), marker="")

b = diagram.line(inB, vB, style="wiggly")
inQ = diagram.line(inQ, vGin)
diagram.line(vGin, vB)
diagram.line(vB, vGout)
diagram.line(vGin, vGout, **gluon_style)
outQ = diagram.line(vGout, outQ)

diagram.plot()
diagram.savefig("nlo-q-v.pdf")
