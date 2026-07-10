from feynman import Diagram

diagram = Diagram(figsize=(10.0, 10.0))
inB = diagram.vertex(xy=(0, 1), marker="")
inQ = diagram.vertex(xy=(0, 0), marker="")
v = diagram.vertex(xy=(0.5, 0.5))
outQ = diagram.vertex(xy=(1, 0.5), marker="")

b = diagram.line(inB, v, style="wiggly")
inQ = diagram.line(inQ, v)
outQ = diagram.line(v, outQ)

diagram.plot()
diagram.savefig("lo-q.pdf")
