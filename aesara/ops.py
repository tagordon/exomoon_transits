import aesara
from aesara.graph.op import Op
from aesara.graph.basic import Apply
from aesara.tensor import as_tensor_variable
from aesara import tensor as tt
import numpy as np

import ctypes
from ctypes import byref
clib = ctypes.CDLL("../fortran/wrapper.so")

#aesara.config.optimizer='None'

clib.grad_impacts.restype = None
clib.flux.restype = None

class FluxOp(Op):
    
    def make_node(self, *args):
        inputs = list(map(as_tensor_variable, args))
        outputs = [
            tt.tensor(
                broadcastable=tuple(inputs[0].broadcastable) + (False,),
                dtype=inputs[0].dtype,
            )
            for _ in range(8)
        ]
        return Apply(self, inputs, outputs)
    
    def perform(self, node, inputs, outputs):
        u1, u2, rp, rm, bp, bpm, theta = inputs
        j = len(bp)
        
        u1, u2, rp, rm = list(
            map(
                lambda a: byref(ctypes.c_double(a)), 
                (u1, u2, rp, rm)
            )
        )
        cth = np.cos(theta)
        sth = np.sin(theta)
        bp, bpm, cth, sth = list(
            map(
                lambda a: byref((ctypes.c_double * j).from_buffer(a)), 
                (bp, bpm, cth, sth)
            )
        )
        
        out = ((ctypes.c_double * 8) * j).from_buffer(np.zeros((8, j)))
        
        clib.flux(u1, u2, rp, rm, bp, bpm, cth, sth, out, byref(ctypes.c_int(j)))
        out = np.array(out).T
        
        for i in range(8):
            outputs[i][0] = out[i]
        
    def grad(self, inputs, gradients):
        f, dfdu1, dfdu2, dfdrp, dfdrm, dfdbp, dfdbpm, dfdtheta = self(*inputs)
        dcdf = gradients
        
        if not isinstance(dcdf[0].type, aesara.gradient.DisconnectedType):
            return [tt.dot(dcdf, dfdx) for dfdx in (dfdu1, dfdu2, dfdrp, dfdrm, dfdbp, dfdbpm, dfdtheta)]
        
    def R_op(self, inputs, eval_points):
        if eval_points[0] is None:
            return eval_points
        return self.grad(inputs, eval_points)
    
    def do_constant_folding(self):
        return False

class ImpactsOp(Op):
        
    def __init__(self, t):
        self.t = t
    
    def make_node(self, *args):
        inputs = list(map(as_tensor_variable, args))
        outputs = [
            tt.tensor(
                broadcastable=tuple(inputs[0].broadcastable) + (False,),
                dtype=inputs[0].dtype,
            )
            for _ in range(45)
        ]
        return Apply(self, inputs, outputs)
        
    def perform(self, node, inputs, outputs):        
        j = len(self.t)
        
        bp = (ctypes.c_double * j).from_buffer(np.zeros(j))
        bpm = (ctypes.c_double * j).from_buffer(np.zeros(j))
        theta = (ctypes.c_double * j).from_buffer(np.zeros(j))
        
        dbp = ((ctypes.c_double * j) * 14).from_buffer(np.zeros((j, 14)))
        dbpm = ((ctypes.c_double * j) * 14).from_buffer(np.zeros((j, 14)))
        dtheta = ((ctypes.c_double * j) * 14).from_buffer(np.zeros((j, 14)))

        args = list(map(lambda a: byref(ctypes.c_double(a)), inputs))
        t = byref((ctypes.c_double * j).from_buffer(self.t))
        clib.grad_impacts(t, *args, byref(ctypes.c_int(j)), 
                        bp, bpm, theta, dbp, dbpm, dtheta)
        
        outputs[0][0] = np.array(bp)        
        outputs[1][0] = np.array(bpm)
        outputs[2][0] = np.array(theta)
        
        for i in range(14):
            outputs[3 + i][0] = np.array(dbp[i])
            outputs[17 + i][0] = np.array(dbpm[i])
            outputs[31 + i][0] = np.array(dtheta[i])
                    
    def grad(self, inputs, gradients):
        outs = self(*inputs)
        
        bp = outs[0]
        bpm = outs[1]
        theta = outs[2]
        
        dbp = [outs[3 + i] for i in range(14)]
        dbpm = [outs[17 + i] for i in range(14)]
        dtheta = [outs[31 + i] for i in range(14)]
        
        dcdf = gradients
        g = [0] * 14
        
        for i in range(14):
            if not isinstance(dcdf[0].type, aesara.gradient.DisconnectedType):
                g[i] += tt.dot(dcdf[0], dbp[i]) 
            if not isinstance(dcdf[1].type, aesara.gradient.DisconnectedType):
                g[i] += tt.dot(dcdf[1], dbpm[i])
            if not isinstance(dcdf[2].type, aesara.gradient.DisconnectedType):
                g[i] += tt.dot(dcdf[2], dtheta[i])
                    
        return g #[tt.sum(x) for x in g]
        
    def R_op(self, inputs, eval_points):
        if eval_points[0] is None:
            return eval_points
        return self.grad(inputs, eval_points)
    
    def do_constant_folding(self):
        return False