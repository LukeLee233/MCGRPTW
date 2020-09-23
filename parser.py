from dataclasses import dataclass
import os


@dataclass
class Edge:
    source: int
    sink: int
    cost: int
    demand: int
    begin_window: int
    end_window: int


def transfer(input_file: str, output_file: str):
    instance_name = os.path.split(input_file)[-1].split('.')[0]
    with open(input_file, 'r') as FILE:
        lines = FILE.readlines()
        nodes_num = int(lines[0].strip())
        edges_num = int(lines[1].strip())

        req_edges = []
        non_req_edges = []
        for idx in range(2, 2 + edges_num):
            edge = list(map(int, lines[idx].split()))
            if len(edge) == 6:
                req_edges.append(Edge(*edge))
            else:
                non_req_edges.append(Edge(*edge, -1, -1))

        capacity = int(lines[-5])
        phi = int(lines[-4])
        psi = int(lines[-3])
        lower_bound = int(lines[-2])
        upper_bound = int(lines[-1])

    with open(output_file, 'w') as fp:
        fp.write(f"Name:\t\t{instance_name}\n")

        if lower_bound == upper_bound:
            fp.write(f"Optimal value:\t{lower_bound}\n", )
        else:
            fp.write(f"Optimal value:\n")

        fp.write("#Vehicles:\t-1\n")
        fp.write(f'Capacity:\t{capacity}\n')
        fp.write(f'Depot Node:\t{1}\n')
        fp.write(f'#Nodes:\t{nodes_num}\n')
        fp.write(f'#Edges:\t{edges_num}\n')
        fp.write(f'#Arcs:\t0\n')
        fp.write(f'#Required N:\t0\n')
        fp.write(f'#Required E:\t{len(req_edges)}\n')
        fp.write(f'#Required A:\t0\n')

        fp.write('\n')
        fp.write('ReN.	DEMAND	S. COST\n')

        fp.write('\n')
        fp.write('ReE.\tFrom N.\tTo N.\tT. COST\tDEMAND\tS. COST\tT. TIME\tS. TIME\tB. TIME\tE. TIME\n')
        for idx, edge in enumerate(req_edges, start=1):
            fp.write(f'E{idx}\t'
                     f'{edge.source + 1}\t'
                     f'{edge.sink + 1}\t'
                     f'{edge.cost}\t'
                     f'{edge.demand}\t'
                     f'{edge.cost}\t'
                     f'{phi * edge.cost}\t'
                     f'{psi * edge.cost}\t'
                     f'{edge.begin_window}\t'
                     f'{edge.end_window}\n')

        fp.write('\n')
        fp.write('EDGE\tFROM N.\tTO N.\tT. COST\tT. TIME\n')
        for idx, edge in enumerate(non_req_edges, start=1):
            fp.write(f'NrE{idx}\t'
                     f'{edge.source + 1}\t'
                     f'{edge.sink + 1}\t'
                     f'{edge.cost}\t'
                     f'{phi * edge.cost}\n')

        fp.write('\n')
        fp.write('ReA.	FROM N.	TO N.	T. COST	DEMAND	S. COST\n')

        fp.write('\n')
        fp.write('ARC	FROM N.	TO N.	T. COST\n')


if __name__ == '__main__':
    input_dir = 'instance/carp_tw_a'
    output_dir = 'instance/carp_tw_a_trans'

    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass

    files = list(filter(lambda x: x.endswith(".dat"), os.listdir(input_dir)))

    for file in files:
        file = os.path.join(input_dir, file)
        name = os.path.split(file)[-1].split('.')[0]
        output = output_dir + '/' + name + "_t.dat"
        transfer(file, output)
