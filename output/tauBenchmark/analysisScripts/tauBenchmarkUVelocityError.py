from math import sin, cos, pi
class TauBenchmarkUVelocityErrorExpression(SimplePythonExpression):
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
    cell_vals = ds_in.GetCellData().GetArray(self.input_var_names[0])
    ncells = ds_in.GetNumberOfCells()
    res = vtk.vtkFloatArray()
    res.SetNumberOfComponents(1)
    res.SetNumberOfTuples(ncells)
    for i in xrange(ncells):
      cell = ds_in.GetCell(i)
      bounds = cell.GetBounds()
      xCenter = (bounds[0] + bounds[1]) / 2
      yCenter = (bounds[2] + bounds[3]) / 2
      val = cos (xCenter) * sin (yCenter)
      val -= cell_vals.GetTuple1(i)
      res.SetTuple1(i, abs(val))
    return res
py_filter = TauBenchmarkUVelocityErrorExpression
