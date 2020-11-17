import json
import os

if __name__ == '__main__':
    filename = "/home/luke/CLionProjects/MCGRP_debug/instance/gurobi_solution/sol_TWe12.json"
    output = "/home/luke/CLionProjects/MCGRP_debug/instance/debug_test/e12.dat"

    with open(output,"w") as wrt:
        with open(filename) as f:
            data = json.load(f)
            routes = data["solution"][0]
            for vehicle in routes:
                tasks = routes[vehicle][1]['service']
                for pair in tasks:
                    wrt.write(f'{pair[0] + 1}->{pair[1] + 1},')
                wrt.seek(wrt.tell() - 1, os.SEEK_SET)
                wrt.write("\n")
