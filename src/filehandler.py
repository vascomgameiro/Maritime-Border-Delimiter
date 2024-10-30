"""
This module contains the FileHandler class, which is responsible for managing file interactions within the program.
"""

import time

import folium

from src.point import Point
from src.valid_points import ValidPoints


class FileHandler:
    """
    The FileHandler ADT is responsible for managing file interactions within the program.
    It handles reading a specified file containing data about points and generating HTML maps
    using those points.

    Attributes:
        __filename (str): The name of the file from which point data will be read.
        __points (ValidPoints): The ValidPoints object made of the file points
        __map_name (None/ str): The name of the map file created

    Methods:
        __init__(str): Initializes a new FileHandler object with the specified file name.
        The file handler is automatically returned at instantiation.

        read_points(): Reads data about points from the specified file and
        creates a ValidPoints object to store the Point objects. Returns the ValidPoints object.

        make_map(Delimitation): Creates an HTML file containing a map
        that represents the points stored in the ValidPoints object as circles,
        and the provided Delimitation object as a line. Uses the folium module to generate the map.

        get_map_name(): Return the name of the map file create

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
        self.__map_name = None

    def read_points(self) -> ValidPoints:
        """
        This method reads point data from a file and stores it in a ValidPoints object.

        Returns: ValidPoints

        Complexity: O(p), p being the number of points stored in the file
        """
        valid_points = []
        open_file = open(self.__filename, "r", encoding="utf-8")
        if open_file:
            for line in open_file:
                items = line.split()
                valid_points.append(Point(items[0], float(items[1]), float(items[2])))
        else:
            print(f"Error opening file {self.__filename}")

        open_file.close()
        return ValidPoints(valid_points)

    def make_map(self, delimitation=None):
        """
        This method creates a HTML file containing a map that represents the points stored
        in the file handler's ValidPoints object as circles, and the delimitation as
        a line.

        Arguments: delimitation(Delimitation)

        Complexity: O(plog(p)), p being the number of points of the ValidPoints object
        """
        folium_map = folium.Map()
        points = self.__points.get_all_points()
        coords = [(point.get_latitude(), point.get_longitude()) for point in points]

        min_lat = self.__points.get_min_lat()
        min_lon = self.__points.get_min_lon()

        max_lat = self.__points.get_max_lat()
        max_lon = self.__points.get_max_lon()

        folium_map.fit_bounds([[min_lat, min_lon], [max_lat, max_lon]])

        for latitude, longitude in coords:
            folium.CircleMarker(location=[latitude, longitude], radius=2, weight=5).add_to(folium_map)

        if delimitation:
            if delimitation.size() > 0:
                del_points = delimitation.get_points()
                del_coordinates = [(del_point.get_latitude(), del_point.get_longitude()) for del_point in del_points]
                del_coordinates.append(
                    (
                        delimitation.get_first().get_latitude(),
                        delimitation.get_first().get_longitude(),
                    )
                )
                folium.PolyLine(locations=del_coordinates, color="red", weight=2.5).add_to(folium_map)

        map_name = f"map_{time.time()}.html"
        self.__map_name = map_name
        folium_map.save(map_name)

    def get_map_name(self) -> str:
        """
        This method returns the name of the map file created.

        Returns: str

        Complexity: O(1)
        """
        return self.__map_name  # type: ignore

    def __repr__(self) -> str:
        """
        This method returns a representation of the FileHandler object.

        Returns: str

        Complexity: O(p), p being the number of points in the FileHandler object
        """
        return f"Filename: {self.__filename}; {self.__points}"
