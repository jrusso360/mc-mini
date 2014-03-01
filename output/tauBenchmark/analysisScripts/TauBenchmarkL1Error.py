from math import sqrt
class TauBenchmarkL1ErrorExpression(SimplePythonExpression):
  def __init__(self):
    SimplePythonExpression.__init__(self)
    self.name = "tauBenchmarkUVelocity"
    self.description = "Analytic solution U Velocity from Tau Benchmark problem."
    self.output_is_point_var = False
    self.output_dimension = 1
  def modify_contract(self,contract):
    pass
  def derive_variable(self,ds_in,domain_id):
    ds_bounds = ds_in.GetBounds()
    x_ext = ds_bounds[1] - ds_bounds[0]
    y_ext = ds_bounds[3] - ds_bounds[2]
    uVals = ds_in.GetCellData().GetArray(self.input_var_names[0])
    vVals = ds_in.GetCellData().GetArray(self.input_var_names[1])
    ncells = ds_in.GetNumberOfCells()
    res = vtk.vtkFloatArray()
    res.SetNumberOfComponents(1)
    res.SetNumberOfTuples(ncells)
    l1Norm = 0
    for i in xrange(ncells):
      cell = ds_in.GetCell(i)
      bounds = cell.GetBounds()
      uVel = uVals.GetTuple1(i)
      vVel = vVals.GetTuple1(i)
      l1Norm += sqrt(uVel * uVel + vVel * vVel) 
    for i in xrange(ncells):
      res.SetTuple1(i, abs(l1Norm / ncells))
    return res
py_filter = TauBenchmarkL1ErrorExpression
