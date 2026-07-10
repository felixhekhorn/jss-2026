import matplotlib.pyplot as plt
from feynman import Diagram

fig = plt.figure(figsize=(10.0, 10.0))
ax = fig.add_axes([0, 0, 1, 1], frameon=False)
gluon_style = dict(style="linear loopy", xamp=0.025, yamp=0.035, nloops=7)

diagram = Diagram(ax)
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

# b.text(r"$b(q)$", t=0.35, y=0.05, fontsize=40)
# outQ.text(r"$\bar q(p_1)$", t=0.5, fontsize=40)
# outQb.text(r"$q(p_2)$", t=0.5, fontsize=40)
# g.text(r"$g(k)$", t=0.35, y=-0.07, fontsize=40)

diagram.plot()
# plt.show()
plt.savefig("nlo-g.pdf")
