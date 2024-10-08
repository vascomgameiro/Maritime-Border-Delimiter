import time

import folium

from src.delimitation import Delimitation
from src.point import Point
from src.valid_points import ValidPoints


class FileHandler:
    """
    The FileHandler ADT is responsible for managing file interactions within the program.
    It handles reading a specified file containing data about points and generating HTML maps
    using those points.

    Attributes:
        __filename (str): The name of the file from which point data will be read.

    Methods:
        __init__(str): Initializes a new FileHandler object with the specified file name.
        The file handler is automatically returned at instantiation.

        read_points(): Reads data about points from the specified file and
        creates a ValidPoints object to store the Point objects. Returns the ValidPoints object.

        make_map(Delimitation): Creates an HTML file containing a map
        that represents the points stored in the ValidPoints object as circles,
        and the provided Delimitation object as a line. Uses the folium module to generate the map.

        __repr__(): Returns a string representation of the FileHandler object, including
        the file name and the points that have been imported.
    """

    def __init__(self, filename: str):
        """
        This method initializes a FileHandler object.

        Arguments: filename(str)

        Complexity: O(1)
        """
        self.__filename = filename
        self.__points = self.read_points()

    def read_points(self) -> ValidPoints:
        """
        This method reads point data from a file and stores it in a ValidPoints object.

        Returns: ValidPoints

        Complexity: O(p), p being the number of points stored in the file
        """
        valid_points = []
        try:
            with open(self.__filename, "r") as md:
                for line in md:
                    items = line.strip().split()
                    if len(items) == 3:
                        valid_points.append(
                            Point(items[0], float(items[1]), float(items[2])),
                        )
        except FileNotFoundError:
            raise FileNotFoundError(f"The file '{self.__filename}' was not found.")
        except ValueError as ve:
            raise ValueError(f"The file '{self.__filename}' has the following value error: {ve}.")
        except Exception as e:
            print(f"Error reading points from {self.__filename}: {e}")
        return ValidPoints(valid_points)

    def make_map(self, delimitation: Delimitation | None):
        """
        This method creates a HTML file containing a map that represents the points stored
        in the file handler's ValidPoints object as circles, and the delimitation as
        a line.

        Arguments: delimitation(Delimitation)

        Complexity: O(p+d), p being the number of points of the ValidPoints object and d
        being the number of points of the Delimitation object
        """
        map = folium.Map()
        points = self.__points.get_all_points()
        coords = [(point.get_longitude(), point.get_latitude()) for point in points]

        try:
            for longitude, latitude in coords:
                folium.CircleMarker(location=[latitude, longitude], radius=2, weight=5).add_to(map)

            if delimitation:
                del_points = delimitation.get_points()
                del_coordinates = [(del_point.get_latitude(), del_point.get_longitude()) for del_point in del_points]
                folium.PolyLine(locations=del_coordinates, color="blue", weight=2.5).add_to(map)

            self.__filename = f"map_{time.time()}.html"
            map.save(self.__filename)
        except Exception as e:
            print(f"Error creating map: {e}")

    def __repr__(self) -> str:
        """
        This method returns a representation of the FileHandler object.

        Returns: str

        Complexity: O(p), p being the number of points in the FileHandler object
        """
        return f"Filename: {self.__filename}\n Points: {self.read_points()}"
