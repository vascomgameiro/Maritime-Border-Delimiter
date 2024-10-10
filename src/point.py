import math

from geopy.distance import geodesic
from geopy.units import nautical


class Vector:
    def __init__(self, point_a: "Point", point_b: "Point") -> None:
        """
        Initializes a Vector object with two points.

        Args:
            point_a (Point): The starting point.
            point_b (Point): The ending point.

        Attributes:
            x (float): The x component of the vector.
            y (float): The y component of the vector.
        Complexity:
            O(1)
        """

        self.__x = point_b.get_proj_x() - point_a.get_proj_x()
        self.__y = point_b.get_proj_y() - point_a.get_proj_y()

    def get_X(self):
        """
        Returns the x component of the vector.

        Returns:
            float: The x component of the vector.
        Complexity:
            O(1)
        """

        return self.__x

    def get_Y(self):
        """
        Returns the y component of the vector.

        Returns:
            float: The y component of the vector.

        Complexity:
            O(1)
        """
        return self.__y

    def __repr__(self) -> str:
        """
        Returns the string representation of the Point object.

        Returns:
            str: String representation of the point.

        Complexity:
            O(1)
        """
        return f"Vector ({self.__x}, {self.__y})"

    def __eq__(self, other: "Vector") -> bool:
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
        return self.__x == other.__x and self.__y == other.__y


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
        earth_radius = 6371009  # in meters
        lat_dist = math.pi * earth_radius / 180.0
        self.__y = nautical(meters=self.get_latitude() * lat_dist)  # type: ignore
        self.__x = nautical(meters=self.get_longitude() * lat_dist * math.cos(math.radians(self.get_latitude())))  # type: ignore

    def get_proj_x(self):
        return self.__x

    def get_proj_y(self):
        return self.__y

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

        vector1 = Vector(self, point_b)
        vector2 = Vector(self, point_c)

        angle = math.degrees(
            math.atan2(vector1.get_Y(), vector1.get_X()) - math.atan2(vector2.get_Y(), vector2.get_X())
        )
        return (360 + angle) % 360

    def __repr__(self) -> str:
        """
        Returns the string representation of the Point object.

        Returns:
            str: String representation of the point.

        Complexity:
            O(1)
        """
        return f"Point {self.__id} ({self.__latitude}, {self.__longitude})"

    def __eq__(self, other: "Point") -> bool:
        """
        Checks if this point is equal to another point.

        Args:
            other (Point): The Point to compare with.

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
