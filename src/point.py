import math
from typing import Tuple

from geopy.distance import geodesic


class Point:
    """
    Represents a point on the globe with an identifier, latitude, and longitude.

    Internal Representation:
    - __id (str): Unique identifier of the point.
    - __latitude (float): Latitude of the point.
    - __longitude (float): Longitude of the point.
    """

    def __init__(self, point_id: str, latitude: float, longitude: float) -> None:
        """
        Initializes a new Point object.

        Args:
            point_id (str): Identifier of the point.
            latitude (float): Latitude of the point.
            longitude (float): Longitude of the point.

        Complexity:
            O(1)
        """
        self.__id = point_id
        self.__latitude = latitude
        self.__longitude = longitude

    def get_id(self) -> str:
        """
        Returns the identifier of the point.

        Returns:
            str: Identifier of the point.

        Complexity:
            O(1)
        """
        return self.__id

    def get_latitude(self) -> float:
        """
        Returns the latitude of the point.

        Returns:
            float: Latitude of the point.

        Complexity:
            O(1)
        """
        return self.__latitude

    def get_longitude(self) -> float:
        """
        Returns the longitude of the point.

        Returns:
            float: Longitude of the point.

        Complexity:
            O(1)
        """
        return self.__longitude

    def distance(self, other: "Point") -> float:
        """
        Calculates the geodesic distance between this point and another point.

        Args:
            other (Point): Another point to calculate the distance to.

        Returns:
            float: Distance in nautical miles between the two points.

        Complexity:
            O(1) - Uses geodesic calculation which is constant time.

        Note:
            The distance is returned in nautical miles.
        """
        coords_1 = (self.__latitude, self.__longitude)
        coords_2 = (other.__latitude, other.__longitude)
        distance = geodesic(coords_1, coords_2).nm

        return distance

    def in_vicinity_of(self, other: "Point", max_distance: int) -> bool:
        """
        Determines if this point is within a certain distance of another point.

        Args:
            other (Point): Another point to compare with.
            max_distance (int): Maximum distance in nautical miles.

        Returns:
            bool: True if within max_distance, False otherwise.

        Complexity:
            O(1)
        """
        return self.distance(other) <= max_distance

    def get_forward_angle(self, point_b: "Point", point_c: "Point") -> float:
        """
        Calculates the forward azimuth angle from this point to point_b,
        then to point_c.

        Args:
            point_b (Point): The second point in the sequence.
            point_c (Point): The third point in the sequence.

        Returns:
            float: The angle in degrees between the lines self->point_b and self->point_c.

        Complexity:
            O(1)

        Note:
            The angle is measured clockwise from the line point_b->self to self->point_c.
        """

        vector_ba = self.__calculate_semi_line(point_b, self)
        vector_ac = self.__calculate_semi_line(self, point_c)

        magnitude_ba = point_b.distance(self)
        magnitude_ac = self.distance(point_c)

        dot_product = self.__calculate_dot_product(vector_ba, vector_ac)
        cos_theta = dot_product / (magnitude_ba * magnitude_ac)

        forward_angle = 360 - math.degrees(math.acos(cos_theta))

        return forward_angle

    @staticmethod
    def __calculate_semi_line(point_A: "Point", point_B: "Point") -> Tuple[float, float]:
        """
        Calculates the direction vector from point_A to point_B.

        Args:
            point_A (Point): The first point.
            point_B (Point): The second point.

        Returns:
            Tuple[float, float]: A tuple representing the direction vector (delta_latitude, delta_longitude).
                                This is the difference between the coordinates of the two points.

        Complexity:
            O(1)
        """
        return (point_A.get_latitude() - point_B.get_latitude(), point_A.get_longitude() - point_B.get_longitude())

    @staticmethod
    def __calculate_dot_product(vector_A: Tuple[float, float], vector_B: Tuple[float, float]) -> float:
        """
        Calculates the dot product of two vectors.

        Args:
            vector_A (Tuple[float, float]): The first vector.
            vector_B (Tuple[float, float]): The second vector.

        Returns:
            float: The dot product of the two vectors.

        Complexity:
            O(1)
        """
        return vector_A[0] * vector_B[0] + vector_A[1] * vector_B[1]

    def __repr__(self) -> str:
        """
        Returns the string representation of the Point object.

        Returns:
            str: String representation of the point.

        Complexity:
            O(1)
        """
        return f"Point {self.__id} ({self.__latitude}, {self.__longitude})"

    def __eq__(self, other: object) -> bool:
        """
        Checks if this point is equal to another point.

        Args:
            other (object): The object to compare with.

        Returns:
            bool: True if equal, False otherwise.

        Complexity:
            O(1)
        """
        if not isinstance(other, Point):
            return False
        return self.__id == other.__id and self.__latitude == other.__latitude and self.__longitude == other.__longitude

    def __hash__(self):
        """
        Returns the hash value of the Point object.

        The hash is based on the identifier, latitude and longitude of the point.

        Returns:
            int: The hash value of the point.

        Complexity:
            O(1)
        """
        return hash((self.__id, self.__latitude, self.__longitude))
