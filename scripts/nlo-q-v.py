import matplotlib.pyplot as plt
from feynman import Diagram

fig = plt.figure(figsize=(10.0, 10.0))
ax = fig.add_axes([0, 0, 1, 1], frameon=False)
gluon_style = dict(
    style="elliptic loopy", ellipse_spread=-0.55, xamp=0.04, yamp=0.05, nloops=7
)

diagram = Diagram(ax)
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

# b.text(r"$b(q)$",t=.2,y=.07,fontsize=40)
# inQ.text(r"$q(p_1)$",t=.3,fontsize=40)
# outQ.text(r"$q'(p_2)$",t=.4,y=.03,fontsize=40)

diagram.plot()
plt.savefig("nlo-q-v.pdf")
