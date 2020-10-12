from dataclasses import dataclass
import os
import openpyxl as xlsx


@dataclass
class Edge:
    source: int
    sink: int
    cost: int
    demand: int
    begin_window: int
    end_window: int


@dataclass
class Task:
    source: int
    sink: int
    travel_time: int
    serve_time: int
    cost: int
    demand: int
    arc: str
    begin_window: int
    end_window: int


def find_files(filename, search_path):
    for root, _, files in os.walk(search_path):
        if filename in files:
            return os.path.join(root,filename)
    return ''


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
            fp.write(f"Optimal value:\t{-1}\n")

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


def transfer_xlsx(filename: str):
    wb = xlsx.load_workbook(filename)
    output_dir = 'instance/istanze/mcgrp_tw'

    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass

    for sheet in wb:
        if sheet.title == 'Istanze MCGRPTW':
            continue

        metafile = find_files(sheet.title[0:2] + '-' + sheet.title[2:] + '.dat', 'instance/istanze')
        with open(metafile, 'r') as FILE:
            lines = FILE.readlines()
            nodes_num = int(lines[0].strip())

            capacity = int(lines[-5])

        data = []
        for row in sheet.iter_rows(min_row=4):
            data_piece = []
            for cell in row:
                data_piece.append(cell.value)
            data.append(Task(*data_piece))

        req_edges = []
        non_req_edges = []
        req_arcs = []
        non_req_arcs = []
        req_nodes = []
        for task in data:
            if task.source == task.sink:
                req_nodes.append(task)
            elif task.demand != 0:
                if task.arc == 'true':
                    req_arcs.append(task)
                else:
                    req_edges.append(task)
            else:
                if task.arc == 'true':
                    non_req_arcs.append(task)
                else:
                    non_req_edges.append(task)

        output_file = output_dir + '/' + sheet.title + "_t.dat"

        with open(output_file, 'w') as fp:
            fp.write(f"Name:\t\t{sheet.title}\n")

            fp.write(f"Optimal value:\t{-1}\n")
            fp.write("#Vehicles:\t-1\n")
            fp.write(f'Capacity:\t{capacity}\n')
            fp.write(f'Depot Node:\t{1}\n')
            fp.write(f'#Nodes:\t{nodes_num}\n')
            fp.write(f'#Edges:\t{len(non_req_edges) + len(req_edges)}\n')
            fp.write(f'#Arcs:\t{len(non_req_arcs) + len(req_arcs)}\n')
            fp.write(f'#Required N:\t{len(req_nodes)}\n')
            fp.write(f'#Required E:\t{len(req_edges)}\n')
            fp.write(f'#Required A:\t{len(req_arcs)}\n')

            fp.write('\n')
            fp.write('ReN.\tDEMAND\tS. COST\tS. TIME\tB. TIME\tE. TIME\n')
            for node in req_nodes:
                fp.write(f'N{node.source + 1}\t'
                         f'{node.demand}\t'
                         f'{node.cost}\t'
                         f'{node.serve_time}\t'
                         f'{node.begin_window}\t'
                         f'{node.end_window}\n')

            fp.write('\n')
            fp.write('ReE.\tFrom N.\tTo N.\tT. COST\tDEMAND\tS. COST\tT. TIME\tS. TIME\tB. TIME\tE. TIME\n')
            for idx, edge in enumerate(req_edges, start=1):
                fp.write(f'E{idx}\t'
                         f'{edge.source + 1}\t'
                         f'{edge.sink + 1}\t'
                         f'{edge.cost}\t'
                         f'{edge.demand}\t'
                         f'{edge.cost}\t'
                         f'{edge.travel_time}\t'
                         f'{edge.serve_time}\t'
                         f'{edge.begin_window}\t'
                         f'{edge.end_window}\n')

            fp.write('\n')
            fp.write('EDGE\tFROM N.\tTO N.\tT. COST\tT. TIME\n')
            for idx, edge in enumerate(non_req_edges, start=1):
                fp.write(f'NrE{idx}\t'
                         f'{edge.source + 1}\t'
                         f'{edge.sink + 1}\t'
                         f'{edge.cost}\t'
                         f'{edge.travel_time}\n')

            fp.write('\n')
            fp.write('ReA.\tFrom N.\tTo N.\tT. COST\tDEMAND\tS. COST\tT. TIME\tS. TIME\tB. TIME\tE. TIME\n')
            for idx, arc in enumerate(req_arcs, start=1):
                fp.write(f'A{idx}\t'
                         f'{arc.source + 1}\t'
                         f'{arc.sink + 1}\t'
                         f'{arc.cost}\t'
                         f'{arc.demand}\t'
                         f'{arc.travel_time}\t'
                         f'{arc.serve_time}\t'
                         f'{arc.begin_window}\t'
                         f'{arc.end_window}\n')

            fp.write('\n')
            fp.write('ARC\tFROM N.\tTO N.\tT. COST\tT. TIME\n')
            for idx, arc in enumerate(non_req_arcs, start=1):
                fp.write(f'NrA{idx}\t'
                         f'{arc.source + 1}\t'
                         f'{arc.sink + 1}\t'
                         f'{arc.cost}\t'
                         f'{arc.travel_time}\n')


if __name__ == '__main__':
    transfer_xlsx('instance/istanze/Istanze ABCes.xlsx')

    # input_dir = 'instance/carp_tw_a'
    # output_dir = 'instance/carp_tw_a_trans'

    # try:
    #     os.mkdir(output_dir)
    # except FileExistsError:
    #     pass
    #
    # files = list(filter(lambda x: x.endswith(".dat"), os.listdir(input_dir)))
    #
    # for file in files:
    #     file = os.path.join(input_dir, file)
    #     name = os.path.split(file)[-1].split('.')[0]
    #     output = output_dir + '/' + name + "_t.dat"
    #     transfer(file, output)
