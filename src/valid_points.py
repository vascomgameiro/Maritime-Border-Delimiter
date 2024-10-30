"""
This module contains the ValidPoints class, which represents a collection of valid points.
"""

from src.point import Point

class ValidPoints:
    """
    Represents a collection of valid points.

    Internal Representation:
    - __points (List[Point]): A list of Point objects.

    Methods:
        - get_size: Returns the number of points stored.
        - get_min_lat, get_max_lat, get_min_lon, get_max_lon: Return the min/max latitude/longitude.
        - get_all_points: Returns all the points stored.
        - get_points_vicinity: Returns the set of points within a given distance from a specified point.
    """

    def __init__(self, points: list):
        """
        Initializes the ValidPoints object with a list of Point objects.

        Args:
            points (List[Point]): A list of valid Point objects.

        Complexity:
            O(P) where P is the number of points in the list.
        """
        self.__points = points

    def get_size(self) -> int:
        """
        Returns the number of points stored.

        Returns:
            int: The number of points.

        Complexity:
            O(1)
        """
        return len(self.__points)

    def get_min_lat(self) -> float:
        """
        Returns the minimum latitude among the points stored.

        Returns:
            float: The minimum latitude.

        Complexity:
            O(P), where P is the number of points.
        """
        return min(point.get_latitude() for point in self.__points)

    def get_max_lat(self) -> float:
        """
        Returns the maximum latitude among the points stored.

        Returns:
            float: The maximum latitude.

        Complexity:
            O(P), where P is the number of points.
        """
        return max(point.get_latitude() for point in self.__points)

    def get_min_lon(self) -> float:
        """
        Returns the minimum longitude among the points stored.

        Returns:
            float: The minimum longitude.

        Complexity:
            O(P), where P is the number of points.
        """
        return min(point.get_longitude() for point in self.__points)

    def get_max_lon(self) -> float:
        """
        Returns the maximum longitude among the points stored.

        Returns:
            float: The maximum longitude.

        Complexity:
            O(P), where P is the number of points.
        """
        return max(point.get_longitude() for point in self.__points)

    def get_all_points(self) -> list:
        """
        Returns the list of all points stored, sorted by their identifier.

        Returns:
            List[Point]: A list of all Point objects.

        Complexity:
            O(P log P) where P is the number of points (due to sorting).
        """
        return sorted(self.__points, key=lambda point: point.get_id())

    def get_points_vicinity(self, ref_point: "Point", max_distance: int) -> set:
        """
        Returns a set of points that are within a given maximum distance from the reference point.

        Args:
            ref_point (Point): The reference point.
            max_distance (int): The maximum distance in nautical miles.

        Returns:
            Set[Point]: A set of points that are within the max_distance from the reference point.

        Complexity:
            O(P), where P is the number of points (since we must check each point).
        """
        return {point for point in self.__points if point.in_vicinity_of(ref_point, max_distance)}

    def __repr__(self) -> str:
        """
        Returns the string representation of the ValidPoints object, showing all points.

        Returns:
            str: A string representation of the points stored.

        Complexity:
            O(P), where P is the number of points.
        """
        return f"ValidPoints({self.__points})"
        return f"ValidPoints({self.__points})"
