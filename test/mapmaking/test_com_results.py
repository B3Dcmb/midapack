#!/usr/bin/python

#Import important modules
import matplotlib.pyplot as plt
import numpy as np


#Define data
procs0 = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
procs1 = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
procs2 = [32, 64, 128, 256, 512, 1024, 2048, 4096]

pre_ring =[15.572764, 7.866721, 4.147319, 2.375836, 1.577882, 1.523011, 1.517357, 2.079866, 4.005907, 9.617817]
pre_butterfly = [15.487716, 7.743700, 3.957187, 1.994612, 1.029207, 0.491761, 0.262946, 0.160273, 0.082454, 0.054311]
pre_nonblocking = [15.461293, 7.808733, 4.096966, 2.307679, 1.536028, 1.270955, 1.453602, 1.985907, 3.232442, 6.658669]
pre_noempty = [15.449770, 7.830982, 4.093966, 2.302852, 1.674734, 1.287849, 1.498801, 2.037193, 3.296579]
pre_alltoallv = [15.474844, 7.810185, 4.094200, 2.302213, 1.599996, 2.087711, 1.466961, 2.543838]
pre_allreduce = [15.390989, 7.745258, 3.849943, 1.937037, 0.985112, 0.450683, 0.224913, 0.118297]
pre_butterfly_blocking1 = [15.459431, 7.783160, 3.931129, 2.002513, 1.035985, 0.471262, 0.245123, 0.134910]
pre_butterfly_blocking2 = [15.458532, 7.796658, 3.954938, 2.025759, 1.084015, 0.493348, 0.260885, 0.143447]
pre_noemptystepring = [15.441426, 7.870011, 4.127127, 2.409728, 1.639390, 1.267909, 1.582869, 2.033446]

pre_ring_2 = [15.512140, 7.865609, 4.167007, 2.699866, 4.857002, 4.778619, 1.466484, 5.398875, 7.984991, 18.962984]
pre_butterfly_2 = [15.426553, 7.767025, 3.974832, 2.242556, 1.173114, 0.633624, 0.253199, 0.139997, 0.086791, 0.170484]
pre_nonblocking_2 = [15.425877, 7.845265, 4.060812, 3.395557, 3.879409, 5.119724, 1.316132, 3.938339, 8.088806, 10.011191]
pre_noempty_2 = [15.425674, 7.834665, 4.060883, 4.023938, 2.836675, 1.145492, 1.964005, 1.843994, 5.680362]
pre_butterfly_blocking1_2 = [15.432156, 7.791317, 3.929135, 2.006978, 1.154006, 0.492901, 0.257471, 0.255673, 0.277881]
pre_butterfly_blocking2_2 = [15.427023, 7.812590, 3.942093, 1.997413, 1.180721, 0.493134, 0.390924, 0.283915, 0.252247]
pre_noemptystepring_2 = [15.424956, 7.828379, 4.095896, 2.235222, 3.042000, 1.879656, 3.505643, 5.525759, 16.269996]

ring = [0.454288, 0.303669, 0.278023, 0.276576, 0.309547, 0.556726, 0.453026, 0.698802, 0.779377, 0.796029]
butterfly = [0.421683, 0.213914, 0.129819, 0.071936, 0.042685, 0.029103, 0.019902, 0.016711, 0.016845, 0.01874]
nonblocking = [0.537412, 0.518241, 0.790133, 1.590417, 2.116983, 2.247586, 3.01164, 5.464721, 9.896763, 20.632889]
noempty = [0.523406, 0.520298, 0.82462, 1.488553, 1.876628, 2.455192, 2.956803, 5.556456, 9.987831]
alltoallv = [0.546789, 0.482164, 0.696876, 0.937968, 1.296889, 1.854625, 3.147458, 5.223848]
allreduce = [0.435936, 0.238148, 0.169611, 0.340098, 0.152843, 0.168728, 0.249569, 0.314614]
butterfly_blocking1 = [0.423522, 0.217502, 0.119248, 0.068527, 0.040777, 0.027492, 0.020293, 0.018633]
butterfly_blocking2 = [0.422853, 0.217992, 0.121326, 0.072131, 0.040063, 0.028567, 0.020408, 0.018821]
noemptystepring = [0.454598, 0.288653, 0.277056, 0.398685, 0.46821, 0.644056, 1.403764, 2.684668]

ring_2 = [0.434830, 0.285190, 0.301304, 1.605559, 2.001266, 3.343376, 0.495055, 0.621897, 0.725988, 0.703469]
butterfly_2 = [0.405223, 0.218798, 0.119759, 0.078155, 0.106737, 0.074618, 0.019405, 0.015921, 0.016246, 0.078614]
nonblocking_2 = [0.507661, 0.501352, 0.957200, 1.426564, 1.983871, 2.369013, 2.926455, 5.046992, 9.104677, 19.564788]
noempty_2 = [0.510111, 0.512390, 0.754336, 1.305708, 2.198547, 2.741100, 2.979942, 5.082322, 9.185337]
butterfly_blocking1_2 = [0.405630, 0.225253, 0.120319, 0.066938, 0.097662, 0.027414, 0.019845, 0.092445, 0.081835]
butterfly_blocking2_2 = [0.405095, 0.224558, 0.120417, 0.118269, 0.114112, 0.026952, 0.087121, 0.048479, 0.077827]
noemptystepring_2 = [0.435016, 0.291484, 0.277280, 0.298172, 1.119666, 0.640839, 5.782456, 9.994717, 13.853718]

#Plot results
import matplotlib
fig1, ax1 = plt.subplots()
# plt.xscale('log', basex=2)
ax1.plot(procs0, pre_ring,'k+--', label= "ring")
ax1.plot(procs0, pre_butterfly,'r+--', label="butterfly")
ax1.plot(procs0, pre_nonblocking,'b+--', label="nonblocking")
ax1.plot(procs1, pre_noempty,'c+--', label="noempty")
ax1.plot(procs2, pre_alltoallv,'m+--', label="alltoallv")
ax1.plot(procs2, pre_allreduce,'y+--', label="allreduce")
ax1.plot(procs2, pre_butterfly_blocking1,'rv--', label="butterfly_blocking_1")
ax1.plot(procs2, pre_butterfly_blocking2,'r^--', label="butterfly_blocking_2")
ax1.plot(procs2, pre_noemptystepring,'g+--', label="no_empty_step_ring")

ax1.set_xscale('log', basex=2)
ax1.set_xticks(procs0)
ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.grid(True)
plt.xlabel("MPI Tasks")
plt.ylabel("exec time (s)")
plt.title("Strong scaling: MatInit. Data volume = 2 x 10^9 nts, nside = 512.")
plt.legend()
plt.show()

fig2, ax2 = plt.subplots()
# plt.xscale('log', basex=2)
ax2.plot(procs0, ring,'k+--', label= "ring")
ax2.plot(procs0, butterfly,'r+--', label="butterfly")
ax2.plot(procs0, nonblocking,'b+--', label="nonblocking")
ax2.plot(procs1, noempty,'c+--', label="noempty")
ax2.plot(procs2, alltoallv,'m+--', label="alltoallv")
ax2.plot(procs2, allreduce,'y+--', label="allreduce")
ax2.plot(procs2, butterfly_blocking1,'rv--', label="butterfly_blocking_1")
ax2.plot(procs2, butterfly_blocking2,'r^--', label="butterfly_blocking_2")
ax2.plot(procs2, noemptystepring,'g+--', label="no_empty_step_ring")

ax2.set_xscale('log', basex=2)
ax2.set_xticks(procs0)
ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax2.set_yscale('log', basey=2)
ax2.set_yticks([0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16])
ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

plt.grid(True)
plt.xlabel("MPI Tasks")
plt.ylabel("exec time (s)")
plt.title("Strong scaling: V = At x d. Data volume = 2 x 10^9 nts, nside = 512.")
plt.legend()
plt.show()

fig3, ax3 = plt.subplots()
ax3.plot(procs0, pre_ring_2,'k+--', label= "ring")
ax3.plot(procs0, pre_butterfly_2,'r+--', label="butterfly")
ax3.plot(procs0, pre_nonblocking_2,'b+--', label="nonblocking")
ax3.plot(procs1, pre_noempty_2,'c+--', label="noempty")
ax3.plot(procs1, pre_butterfly_blocking1_2,'rv--', label="butterfly_blocking_1")
ax3.plot(procs1, pre_butterfly_blocking2_2,'r^--', label="butterfly_blocking_2")
ax3.plot(procs1, pre_noemptystepring_2,'g+--', label="no_empty_step_ring")

ax3.set_xscale('log', basex=2)
ax3.set_xticks(procs0)
ax3.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.grid(True)
plt.xlabel("MPI Tasks")
plt.ylabel("exec time (s)")
plt.title("Second Trial Strong scaling: MatInit. Data volume = 2 x 10^9 nts, nside = 512.")
plt.legend()
plt.show()

fig4, ax4 = plt.subplots()
ax4.plot(procs0, ring_2,'k+--', label= "ring")
ax4.plot(procs0, butterfly_2,'r+--', label="butterfly")
ax4.plot(procs0, nonblocking_2,'b+--', label="nonblocking")
ax4.plot(procs1, noempty_2,'c+--', label="noempty")
ax4.plot(procs1, butterfly_blocking1_2,'rv--', label="butterfly_blocking_1")
ax4.plot(procs1, butterfly_blocking2_2,'r^--', label="butterfly_blocking_2")
ax4.plot(procs1, noemptystepring_2,'g+--', label="no_empty_step_ring")

ax4.set_xscale('log', basex=2)
ax4.set_xticks(procs0)
ax4.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax4.set_yscale('log', basey=2)
ax4.set_yticks([0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16])
ax4.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

plt.grid(True)
plt.xlabel("MPI Tasks")
plt.ylabel("exec time (s)")
plt.title("Second Trial Strong scaling: V = At x d. Data volume = 2 x 10^9 nts, nside = 512.")
plt.legend()
plt.show()
