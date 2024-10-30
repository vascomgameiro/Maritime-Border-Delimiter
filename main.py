"""
This module contains classes and methods for processing geographical points and calculating distances.
"""

import time

from src.approx_delimitation import ApproxDelimitation
from src.convex_hull import ConvexHull
from src.filehandler import FileHandler
from src.optimal_delimitation import OptimalDelimitation


def main():
    """
    Main function of the program, which will continuously prompt the user for input commands.
    Commands available are:
    - import_points <file_name> : Load points from a file.
    - draw_hull : Draw the convex hull of the imported points.
    - draw_optimal_delimitation <max_distance> : Draw optimal delimitation for the specified max distance.
    - draw_approx_delimitation <max_distance> : Draw approximate delimitation for the specified max distance.
    - make_map : Generate a map of the points with the delimitation.
    - quit : Exit the program.
    """
    file_handler = None
    delimitation = None

    while True:
        command = input().strip().split()
        if len(command) == 0:
            continue
        cmd = command[0]

        if cmd == "import_points":
            if len(command) < 2:
                print("Error: missing file name.")
                continue
            file_name = command[1]

            start_time = time.time()
            file_handler = FileHandler(file_name)
            num_points = file_handler.read_points().get_size()
            elapsed_time = (time.time() - start_time) * 1000

            print(f"Size: {num_points}")
            print(f"Time: {elapsed_time:.2f}\n")
        elif cmd == "help":
            print("Commands available are:")
            print("- import_points <file_name> : Load points from a file.")
            print("- draw_hull : Draw the convex hull of the imported points.")
            print(
                "- draw_optimal_delimitation <max_distance> : Draw optimal delimitation for the specified max distance."
            )
            print(
                "- draw_approx_delimitation <max_distance> : Draw approximate delimitation for the specified max distance."
            )
            print("- make_map : Generate a map of the points with the delimitation.")
            print("- quit : Exit the program.\n")
        elif cmd == "draw_hull":
            if file_handler is None:
                print("Error: No points imported. Please use 'import_points' first.")
                continue
            start_time = time.time()
            delimitation = ConvexHull(file_handler.read_points()).find_delimitation()
            try:
                area = delimitation.get_area()
            except ValueError:
                area = 0
            elapsed_time = (time.time() - start_time) * 1000

            print("Type: Convex Hull")
            print(f"Line: {delimitation}")
            print(f"Area: {area:.2f}")
            print(f"Time: {elapsed_time:.2f}\n")

        elif cmd == "draw_optimal_delimitation":
            if len(command) < 2:
                print("Error: missing max distance argument.")
                continue
            if file_handler is None:
                print("Error: No points imported. Please use 'import_points' first.")
                continue

            max_distance = int(command[1])
            start_time = time.time()
            delimitation = OptimalDelimitation(file_handler.read_points(), max_distance).find_delimitation()
            elapsed_time = (time.time() - start_time) * 1000
            try:
                area = delimitation.get_area()
            except ValueError:
                area = 0

            print("Type: Optimal Delimitation")
            print(f"Line: {delimitation}")
            print(f"Area: {area:.2f}")
            print(f"Time: {elapsed_time:.2f}\n")

        elif cmd == "draw_approx_delimitation":
            if len(command) < 2:
                print("Error: missing max distance argument.")
                continue
            if file_handler is None:
                print("Error: No points imported. Please use 'import_points' first.")
                continue

            max_distance = int(command[1])
            start_time = time.time()
            delimitation = ApproxDelimitation(file_handler.read_points(), max_distance).find_delimitation()
            elapsed_time = (time.time() - start_time) * 1000
            try:
                area = delimitation.get_area()
            except ValueError:
                area = 0

            print("Type: Approximate Delimitation")
            print(f"Line: {delimitation}")
            print(f"Area: {area:.2f}")
            print(f"Time: {elapsed_time:.2f}\n")

        elif cmd == "make_map":
            if file_handler is None:
                print("Error: No points available.")
                continue
            start_time = time.time()
            file_handler.make_map(delimitation)
            elapsed_time = (time.time() - start_time) * 1000
            print(f"Map: {file_handler.get_map_name()}")
            print(f"Time: {elapsed_time:.2f}\n")

        elif cmd == "quit":
            print("Done!\n")
            break

        else:
            print(
                f"Unknown command: {cmd}. Valid Commands: import_points, draw_hull, draw_optimal_delimitation <max_distance>, draw_approx_delimitation <max_distance>, make_map, quit"
            )


if __name__ == "__main__":
    main()
