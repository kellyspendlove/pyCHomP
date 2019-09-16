### LapGradedComplex.py
### MIT LICENSE 2019 Kelly Spendlove

from pychomp._chomp import *
from pychomp.CondensationGraph import *
from pychomp.StronglyConnectedComponents import *
from pychomp.DirectedAcyclicGraph import *
from pychomp.Poset import *

def LapGradedComplex(complex, laps):
  """
  Overview:
    Given a complex and a function on its top dimensional cells,
    produce a GradedComplex such that the preimage of a down set
    is the collection of cells in the closure of all the 
    associated top cells
  Inputs:
    complex    : a complex
    laps : a function from vertices to the integers
  """
  grading = construct_grading(complex, laps)
  return GradedComplex(complex, grading) # lambda x : mapping[x])

  #return poset, chompy.GradedComplex(complex, lambda x : mapping[x])

#Produce total order given lap numbers (used for visualization)
def TotalOrder(complex, laps):
  MI = min([laps(cell) for cell in complex(complex.dimension())])
  MA = max([laps(cell) for cell in complex(complex.dimension())]) 
  dag = DirectedAcyclicGraph()
  for cell in complex(complex.dimension()):
      v = laps(cell)
      #Don't add maximum value, this handles the right fringe
      if laps(cell) < MA:
          dag.add_vertex(v)
  for v in dag.vertices():
      if v>MI:
          dag.add_edge(v,v-1)
  return dag

