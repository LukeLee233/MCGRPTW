import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt('/home/luke/CLionProjects/MCGRP_debug/instance/debug_test/log/2020-Dec-15 17:06:15/TWe10_t.dat.log',
                     delimiter=',')

print(data.shape)

# cost_len = data[:, 0]
# vio_load = data[:, 1]
# beta = data[:, 2]
# fitness = data[:, 3]
#
# print(fitness.min())
#
# # x = np.ones(600)*3000
# fig = plt.figure()
# # plt.title(f'')
#
# plt.ylabel('total cost')
# plt.xlabel('searching step')
# plt.plot(beta, label='beta')
# # plt.plot(x)
#
# plt.legend(loc='upper right')
# plt.show()
#
# x = np.ones(600)*3000
# fig = plt.figure()
# plt.title(f'')

plt.ylabel('total cost')
plt.xlabel('searching step')
plt.plot(data, label='local search')
# plt.axhline(y=12000, color='r', label='milestone')
# plt.axhline(y=11794, color='g', label='best')

# plt.plot(x)

plt.legend(loc='upper right')
plt.show()

# fig.savefig('/home/luke/CLionProjects/MCGRP_debug/instance/bhw_test/local threshold comparison.pdf')
