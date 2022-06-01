#!/usr/bin/env python

import numpy as np
import optimization_module

layer = optimization_module.generate_multilayer(0.076135, 0.054135, 22, 0.00119)
optimization_module.plot_points(layer['x'], layer['y'])
optimization_module.generate_input_file('test.txt', layer['x'], layer['y'], 10.)
