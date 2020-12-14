import json

filename = '/home/luke/CLionProjects/MCGRP_debug/Visualize/result/sol_TWe10.json'


class Solution:
    node_seq = []


with open(filename, 'r') as f:
    data = json.load(f)

results = data['solution'][0]

# output node level info
for route in results.values():
    seq = route[0]['route']
    print('{', end='')
    for i in seq:
        print(f'{i + 1}', end=',')
    print('\b',end='')
    print('},')


print('\n'* 3)

for route in results.values():
    seq = route[1]['service']
    print('{', end='')
    for pair in seq:
        print('{' + f'{pair[0] + 1},{pair[1] + 1}' + '}', end=',')
    print('\b',end='')
    print('},')