import aesara
from aesara.graph.op import Op
from aesara.graph.basic import Apply
from aesara.tensor import as_tensor_variable

class Impacts(Op):
    
    __props__ = () 
    
    def make_node(self, *inargs):
        in_args = map(as_tensor_variable, args)
        return Apply(self, [x], [x.type()] * 15)
        
    def perform(#...)
        # this can include the derivatives, we'll get them in the gradient function by calling self
        


class DoubleOp1(Op):
    __props__ = ()

    def make_node(self, x):
        x = aesara.tensor.as_tensor_variable(x)
        # Note: using x_.type() is dangerous, as it copies x's broadcasting
        # behaviour
        return Apply(self, [x], [x.type()])

    def perform(self, node, inputs, output_storage):
        x = inputs[0]
        z = output_storage[0]
        z[0] = x * 2

    def infer_shape(self, fgraph, node, i0_shapes):
        return i0_shapes

    def grad(self, inputs, output_grads):
        return [output_grads[0] * 2]

    def R_op(self, inputs, eval_points):
        # R_op can receive None as eval_points.
        # That mean there is no diferientiable path through that input
        # If this imply that you cannot compute some outputs,
        # return None for those.
        if eval_points[0] is None:
            return eval_points
        return self.grad(inputs, eval_points)

doubleOp1 = DoubleOp1()