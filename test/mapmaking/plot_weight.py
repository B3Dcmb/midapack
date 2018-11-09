#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

lambdareduce = 10
lambdaBlock = 2**4
ntimesamples = 20

W = np.zeros(ntimesamples)
W[0] = 10
state = np.random.RandomState(lambdaBlock)
for i in range(1,lambdareduce):
    W[i] = -1./i
W[lambdareduce:lambdaBlock] = 0.01 * state.rand(lambdaBlock-lambdareduce)
print(W[:lambdaBlock])
plt.plot(np.arange(ntimesamples),W,'r.-')
plt.title("Weight profile")
plt.grid(True)
plt.show()
