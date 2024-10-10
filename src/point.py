import math

from geopy.distance import geodesic
from geopy.units import nautical

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

        vector1=(point_b.get_proj_x() - self.get_proj_x(), point_b.get_proj_y() - self.get_proj_y())
        
        vector2=(point_c.get_proj_x() - self.get_proj_x(), point_c.get_proj_y() - self.get_proj_y())

        angle = math.degrees(
            math.atan2(vector1[1], vector1[0]) - math.atan2(vector2[1], vector2[0])
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
